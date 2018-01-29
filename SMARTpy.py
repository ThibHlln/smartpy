from datetime import timedelta, datetime
from os import path, makedirs, sep
import argparse

from SMARTfiles import \
    get_dict_rain_series_simu, get_dict_peva_series_simu, get_dict_discharge_series, write_flow_file_from_dict
from SMARTparameters import get_parameters_from_file
from SMARTstructure import run
from SMARTtime import TimeFrame
from SMARTobjective import calculate_obj_fn


class SMART(object):
    def __init__(self, catchment, area_m2, start, end, time_delta_simu, time_delta_report, warm_up_days, root):
        # general information
        self.catchment = catchment
        self.area = area_m2
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
        self.rain = get_dict_rain_series_simu(''.join([self.in_f, self.catchment, '.rain']),
                                              self.timeseries[1], self.timeseries[-1], self.delta_simu)
        self.peva = get_dict_peva_series_simu(''.join([self.in_f, self.catchment, '.peva']),
                                              self.timeseries[1], self.timeseries[-1], self.delta_simu)
        self.flow = get_dict_discharge_series(''.join([self.in_f, self.catchment, '.flow']),
                                              self.timeframe.get_series_report()[1],
                                              self.timeframe.get_series_report()[-1])

    def simulate(self, parameters):
        db = dict()
        return run(self.area, self.delta_simu, self.rain, self.peva,
                   parameters, db, self.timeseries, self.timeseries_report,
                   warm_up=self.warm_up)


def simulate(catchment_id, area_m2,
             start, end, time_delta_simu, time_delta_report, warm_up_days,
             root):

    smart_model = SMART(catchment_id, area_m2, start, end, time_delta_simu, time_delta_report, warm_up_days, root)

    # get sets of parameters
    parameters = get_parameters_from_file(''.join([smart_model.in_f, smart_model.catchment, '.parameters']))

    # run simulations for each set of parameters
    db = dict()

    # simulate
    dict_discharge, gw_ratio = smart_model.simulate(parameters)

    write_flow_file_from_dict(smart_model.timeframe, dict_discharge,
                              ''.join([smart_model.out_f, catchment_id, '.flow']),
                              method='raw')

    # calculate objective function
    obs = [observation for observation in smart_model.flow.itervalues()]
    mod = [dict_discharge[dt] for dt in smart_model.flow.iterkeys()]

    nse = calculate_obj_fn(mod, obs, method='nse')
    print nse, gw_ratio

    # explicit garbage collection
    del db

    # save results in file


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


if __name__ == '__main__':
    # Define the root of the CSF
    smart_root = path.realpath('..')  # move to parent directory of this current python file

    # Collect the arguments of the program call
    parser = argparse.ArgumentParser(description="simulate lumped catchment hydrology"
                                                 "for one catchment and one time period")
    parser.add_argument('catchment', type=str,
                        help="name of the catchment")
    parser.add_argument('outlet', type=str,
                        help="european code of the catchment outlet [format IE_XX_##X######]")
    parser.add_argument('start', type=valid_date,
                        help="simulation start date and time [format DD/MM/YYYY_HH:MM:SS]")
    parser.add_argument('end', type=valid_date,
                        help="simulation end date and time [format DD/MM/YYYY_HH:MM:SS]")
    parser.add_argument('delta_simu', type=valid_delta,
                        help="time delta for simulation in minutes")
    parser.add_argument('delta_report', type=valid_delta,
                        help="time delta for reporting in minutes")
    parser.add_argument('area', type=float,
                        help="catchment area in square kilometers")
    parser.add_argument('-w', '--warm_up', type=int, default=365,
                        help="warm-up duration in days")

    args = parser.parse_args()

    # Run the main() function
    simulate('_'.join([args.catchment.capitalize(), args.outlet.upper()]), args.area * 1E6,
             args.start, args.end, args.delta_simu, args.delta_report, args.warm_up,
             smart_root)
