from datetime import timedelta, datetime
from os import path, makedirs, sep
import imp
import argparse

from SMARTinout import \
    get_dict_rain_series_simu, get_dict_peva_series_simu, get_dict_discharge_series, write_flow_file_from_dict
from SMARTparameters import get_parameters_from_file
from SMARTstructure import run
from SMARTtime import TimeFrame


class SMART(object):
    def __init__(self, catchment, c_area_m2, g_area_m2, start, end,
                 time_delta_simu, time_delta_report, warm_up_days, in_fmt, root):
        # general information
        self.catchment = catchment
        self.area = c_area_m2
        # directory information
        self.root_f = root
        self.in_f = sep.join([self.root_f, 'in', self.catchment, sep])
        self.out_f = sep.join([self.root_f, 'out', self.catchment, sep])
        if not path.exists(self.out_f):
            makedirs(self.out_f)
        # temporal information
        self.start = start
        self.end = end
        self.delta_simu = time_delta_simu
        self.delta_report = time_delta_report
        self.timeframe = TimeFrame(self.start, self.end, self.delta_simu, self.delta_report)
        self.timeseries = self.timeframe.get_series_simu()
        self.timeseries_report = self.timeframe.get_series_report()
        self.warm_up = warm_up_days
        # physical information
        extra_ext = '.nc' if in_fmt == 'netcdf' else ''
        self.rain = get_dict_rain_series_simu(''.join([self.in_f, self.catchment, '.rain' + extra_ext]), in_fmt,
                                              self.timeseries[1], self.timeseries[-1], self.delta_simu)
        self.peva = get_dict_peva_series_simu(''.join([self.in_f, self.catchment, '.peva' + extra_ext]), in_fmt,
                                              self.timeseries[1], self.timeseries[-1], self.delta_simu)
        self.flow = get_dict_discharge_series(''.join([self.in_f, self.catchment, '.flow']),
                                              self.timeframe.get_series_report()[1],
                                              self.timeframe.get_series_report()[-1],
                                              c_area_m2, g_area_m2)

    def simulate(self, parameters):
        db = dict()
        return run(self.area, self.delta_simu, self.rain, self.peva,
                   parameters, db, self.timeseries, self.timeseries_report,
                   warm_up=self.warm_up)


def simulate(catchment_id, c_area_m2, g_area_m2,
             start, end, time_delta_simu, time_delta_report, warm_up_days,
             in_format, root):

    smart_model = SMART(catchment_id, c_area_m2, g_area_m2, start, end,
                        time_delta_simu, time_delta_report, warm_up_days, in_format, root)

    # get sets of parameters
    parameters = get_parameters_from_file(''.join([smart_model.in_f, smart_model.catchment, '.parameters']))

    # run simulations for each set of parameters
    db = dict()

    # simulate
    dict_discharge, gw_ratio = smart_model.simulate(parameters)

    write_flow_file_from_dict(smart_model.timeframe, dict_discharge,
                              ''.join([smart_model.out_f, catchment_id, '.mod.flow']),
                              method='raw')

    write_flow_file_from_dict(smart_model.timeframe, smart_model.flow,
                              ''.join([smart_model.out_f, catchment_id, '.obs.flow']),
                              method='raw')

    # explicit garbage collection
    del db


def valid_date(s):
    try:
        return datetime.strptime(s, "%d/%m/%Y_%H:%M:%S")
    except ValueError:
        raise argparse.ArgumentTypeError("Not a valid date: '{0}'.".format(s))


def valid_delta(n):
    try:
        return timedelta(minutes=int(n))
    except ValueError:
        raise argparse.ArgumentTypeError("Not a valid time delta: '{0}'.".format(n))


def valid_file_format(fmt):
    if fmt.lower() == "netcdf":
        try:
            imp.find_module('netCDF4')
            return "netcdf"
        except ImportError:
            raise argparse.ArgumentTypeError("NetCDF4 module is not installed, please choose another file format.")
    elif fmt.lower() == "csv":
        try:
            imp.find_module('csv')
            return "csv"
        except ImportError:
            raise argparse.ArgumentTypeError("CSV module is not installed, please choose another file format.")
    else:
        raise argparse.ArgumentTypeError("File format not recognised: '{0}'.".format(fmt))


if __name__ == '__main__':
    # Define the root of the SMARTpy package
    smart_root = path.realpath('..')  # move to parent directory of this current python file

    # Collect the arguments of the program call
    parser = argparse.ArgumentParser(description="simulate lumped catchment hydrology"
                                                 "for one catchment and one time period")
    parser.add_argument('catchment', type=str,
                        help="name of the catchment")
    parser.add_argument('start', type=valid_date,
                        help="simulation start date and time [format DD/MM/YYYY_HH:MM:SS]")
    parser.add_argument('end', type=valid_date,
                        help="simulation end date and time [format DD/MM/YYYY_HH:MM:SS]")
    parser.add_argument('delta_simu', type=valid_delta,
                        help="time delta for simulation in minutes")
    parser.add_argument('delta_report', type=valid_delta,
                        help="time delta for reporting in minutes")
    parser.add_argument('catchment_area', type=float,
                        help="catchment area in square kilometers")
    parser.add_argument('gauged_area', type=float,
                        help="area upstream of the hydrometric gauge in square kilometers")
    parser.add_argument('-w', '--warm_up', type=int, default=365,
                        help="warm-up duration in days")
    parser.add_argument('-i', '--in_format', type=valid_file_format, default='csv',
                        help="format of input data files [csv or netcdf]")

    args = parser.parse_args()

    # Run the main() function
    simulate(args.catchment,
             args.catchment_area * 1E6, args.gauged_area * 1E6,
             args.start, args.end, args.delta_simu, args.delta_report, args.warm_up,
             args.in_format, smart_root)
