# This file is part of SMARTpy - An open-source rainfall-runoff model in Python
# Copyright (C) 2018-2022  Thibault Hallouin (1)
#
# (1) Dooge Centre for Water Resources Research, University College Dublin, Ireland
#
# SMARTpy is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SMARTpy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with SMARTpy. If not, see <http://www.gnu.org/licenses/>.

from builtins import dict, zip
from csv import DictReader, writer
from datetime import datetime, timedelta
import numpy as np
from collections import OrderedDict
import argparse
try:
    from netCDF4 import Dataset
except ImportError:
    Dataset = None

from .timeframe import get_required_resolution, check_interval_in_list, \
    rescale_time_resolution_of_irregular_mean_data, rescale_time_resolution_of_regular_cumulative_data
from .version import __version__


def get_dict_rain_series_simu(file_location, file_format, start_simu, end_simu, time_delta_simu):
    dict_rain, start_data, end_data, time_delta_data = read_rain_file(file_location, file_format)

    if (start_data - time_delta_data + time_delta_simu <= start_simu) and (end_simu <= end_data):
        time_delta_res = get_required_resolution(start_data, start_simu, time_delta_data, time_delta_simu)
        return rescale_time_resolution_of_regular_cumulative_data(dict_rain,
                                                                  start_data, end_data, time_delta_data,
                                                                  time_delta_res,
                                                                  start_simu, end_simu, time_delta_simu)
    else:
        raise Exception('Rain data not sufficient for simulation.')


def get_dict_peva_series_simu(file_location, file_format, start_simu, end_simu, time_delta_simu):
    dict_peva, start_data, end_data, time_delta_data = read_peva_file(file_location, file_format)

    if (start_data - time_delta_data + time_delta_simu <= start_simu) and (end_simu <= end_data):
        time_delta_res = get_required_resolution(start_data, start_simu, time_delta_data, time_delta_simu)
        return rescale_time_resolution_of_regular_cumulative_data(dict_peva,
                                                                  start_data, end_data, time_delta_data,
                                                                  time_delta_res,
                                                                  start_simu, end_simu, time_delta_simu)
    else:
        raise Exception('PEva data not sufficient for simulation.')


def get_dict_discharge_series(file_location, file_format, start_report, end_report, catchment_area, gauged_area):
    data_flow = read_flow_file(file_location, file_format)

    scaling_factor = catchment_area / gauged_area

    start_date = (start_report - timedelta(days=2)).date()
    end_date = (end_report + timedelta(days=1)).date()

    # select subset of observations in simulation period + apply scaling factor
    dict_flow = OrderedDict()
    for dt in data_flow:
        d = dt.date()
        if (start_date <= d) and (d <= end_date):
            dict_flow[dt] = data_flow[dt] * scaling_factor

    # return rescaled observations time series
    return rescale_time_resolution_of_irregular_mean_data(dict_flow, start_report, end_report,
                                                          timedelta(days=1), timedelta(hours=1))


def get_dict_simulation_settings(file_location):
    my_dict_args = read_simulation_settings_file(file_location)
    # CATCHMENT AREA [float, m2]
    try:
        c_area = float(my_dict_args["catchment_area_km2"]) * 1e6
    except KeyError:
        raise Exception('Setting CATCHMENT AREA is missing from simulation file.')
    except ValueError:
        raise Exception('Setting CATCHMENT AREA could not be converted to a float.')
    # GAUGED AREA [float, m2]
    try:
        g_area = float(my_dict_args["gauged_area_km2"]) * 1e6
    except KeyError:
        g_area = c_area  # i.e. there is no gauge or area is the same
    except ValueError:
        raise Exception('Setting GAUGED AREA could not be converted to a float.')
    # START SIMULATION [datetime]
    try:
        start = datetime.strptime(my_dict_args["start_datetime"], '%d/%m/%Y %H:%M:%S')
    except KeyError:
        raise Exception('Setting START is missing from simulation file.')
    except ValueError:
        raise Exception('Setting START could not be converted to a datetime [format required: DD/MM/YYYY HH:MM:SS].')
    # END SIMULATION [datetime]
    try:
        end = datetime.strptime(my_dict_args["end_datetime"], '%d/%m/%Y %H:%M:%S')
    except KeyError:
        raise Exception('Setting END is missing from simulation file.')
    except ValueError:
        raise Exception('Setting END could not be converted to a datetime [format required: DD/MM/YYYY HH:MM:SS].')
    # SIMULATION TIME STEP [timedelta]
    try:
        delta_simu = timedelta(minutes=int(my_dict_args["simu_timedelta_min"]))
    except KeyError:
        raise Exception('Setting DELTA SIMU is missing from simulation file.')
    except ValueError:
        raise Exception('Setting DELTA SIMU could not be converted to an integer/timedelta.')
    # REPORTING TIME STEP [timedelta]
    try:
        delta_report = timedelta(minutes=int(my_dict_args["report_timedelta_min"]))
    except KeyError:
        raise Exception('Setting DELTA REPORT is missing from simulation file.')
    except ValueError:
        raise Exception('Setting DELTA REPORT could not be converted to an integer/timedelta.')
    # WARM-UP DURATION [int, days]
    try:
        warm_up = int(my_dict_args["warm_up_days"])
    except KeyError:
        raise Exception('Setting WARM UP DURATION is missing from simulation file.')
    except ValueError:
        raise Exception('Setting WARM UP DURATION could not be converted to an integer.')
    # GROUNDWATER CONSTRAINT [float]
    try:
        gw_constraint = float(my_dict_args["gw_constraint"])
    except KeyError:
        gw_constraint = None  # i.e. no simulation constraint required
    except ValueError:
        raise Exception('Setting GROUNDWATER CONSTRAINT could not be converted to a float.')

    return c_area, g_area, start, end, delta_simu, delta_report, warm_up, gw_constraint


def read_rain_file(file_location, file_format):
    if file_format == 'netcdf':
        if Dataset:
            return read_netcdf_time_series_with_delta_check(file_location, key_variable='DateTime', val_variable='rain')
        else:
            raise Exception("The use of 'netcdf' as the input file format requires the package 'netCDF4', "
                            "please install it and retry, or choose another file format.")
    else:
        return read_csv_time_series_with_delta_check(file_location, key_header='DateTime', val_header='rain')


def read_peva_file(file_location, file_format):
    if file_format == 'netcdf':
        if Dataset:
            return read_netcdf_time_series_with_delta_check(file_location, key_variable='DateTime', val_variable='peva')
        else:
            raise Exception("The use of 'netcdf' as the input file format requires the package 'netCDF4', "
                            "please install it and retry, or choose another file format.")
    else:
        return read_csv_time_series_with_delta_check(file_location, key_header='DateTime', val_header='peva')


def read_flow_file(file_location, file_format):
    if file_format == 'netcdf':
        if Dataset:
            return read_netcdf_time_series_with_missing_check(file_location,
                                                              key_variable='DateTime', val_variable='flow')
        else:
            raise Exception("The use of 'netcdf' as the input file format requires the package 'netCDF4', "
                            "please install it and retry, or choose another file format.")
    else:
        return read_csv_time_series_with_missing_check(file_location, key_header='DateTime', val_header='flow')


def read_simulation_settings_file(file_location):
    my_dict_args = dict()
    try:
        with open(file_location, 'r', encoding='utf8') as my_file:
            my_reader = DictReader(my_file)
            for row in my_reader:
                my_dict_args[row['ARGUMENT']] = row['VALUE']
    except KeyError:
        raise Exception("There is no 'ARGUMENT' or 'VALUE' column in {}.".format(file_location))
    except IOError:
        raise Exception("There is no simulation file at {}.".format(file_location))

    return my_dict_args


def read_csv_time_series_with_delta_check(csv_file, key_header, val_header):
    try:
        with open(csv_file, 'r', encoding='utf8') as my_file:
            my_dict_data = dict()
            my_list_dt = list()
            my_reader = DictReader(my_file)
            try:
                for row in my_reader:
                    my_dict_data[datetime.strptime(row[key_header], "%Y-%m-%d %H:%M:%S")] = np.float64(row[val_header])
                    my_list_dt.append(datetime.strptime(row[key_header], "%Y-%m-%d %H:%M:%S"))
            except KeyError:
                raise Exception('Field {} or {} does not exist in {}.'.format(key_header, val_header, csv_file))

        start_data, end_data, time_delta = check_interval_in_list(my_list_dt, csv_file)

        return my_dict_data, start_data, end_data, time_delta
    except IOError:
        raise Exception('File {} could not be found.'.format(csv_file))


def read_netcdf_time_series_with_delta_check(netcdf_file, key_variable, val_variable):
    try:
        with Dataset(netcdf_file, "r") as my_file:
            my_file.set_auto_mask(False)
            my_dict_data = dict()
            my_list_dt = list()
            try:
                my_list_dt += \
                    [datetime(1970, 1, 1) + timedelta(seconds=tstamp) for tstamp in my_file.variables[key_variable][:]]
                for idx, dt in enumerate(my_list_dt):
                    my_dict_data[dt] = my_file.variables[val_variable][idx]
            except KeyError:
                raise Exception('Variable {} or {} does not exist in {}.'.format(key_variable, val_variable,
                                                                                 netcdf_file))

        start_data, end_data, time_delta = check_interval_in_list(my_list_dt, netcdf_file)

        return my_dict_data, start_data, end_data, time_delta
    except IOError:
        raise Exception('File {} could not be found.'.format(netcdf_file))


def read_csv_time_series_with_missing_check(csv_file, key_header, val_header):
    try:
        with open(csv_file, 'r', encoding='utf8') as my_file:
            my_dict_data = OrderedDict()
            my_reader = DictReader(my_file)
            try:
                for row in my_reader:
                    try:
                        if row[val_header] != '':  # not an empty string (that would mean missing data)
                            if np.float64(row[val_header]) != -99.0:  # flag for missing data
                                my_dict_data[datetime.strptime(row[key_header], "%Y-%m-%d %H:%M:%S")] = \
                                    np.float64(row[val_header])
                    except ValueError:
                        raise Exception('Field {} in {} cannot be converted to float '
                                        'at {}.'.format(val_header, csv_file, row[key_header]))
            except KeyError:
                raise Exception('Field {} or {} does not exist in {}.'.format(key_header, val_header, csv_file))
        return my_dict_data
    except IOError:
        raise Exception('File {} could not be found.'.format(csv_file))


def read_netcdf_time_series_with_missing_check(netcdf_file, key_variable, val_variable):
    try:
        with Dataset(netcdf_file, 'r') as my_file:
            my_dict_data = OrderedDict()
            try:
                my_dts = \
                    [datetime(1970, 1, 1) + timedelta(seconds=tstamp) for tstamp in my_file.variables[key_variable][:]]
                my_flows = my_file.variables[val_variable][:]

                for dt, flow in zip(my_dts, my_flows):
                    if not np.isnan(flow):  # flag for missing data
                        my_dict_data[dt] = flow

            except KeyError:
                raise Exception('Variable {} or {} does not exist in {}.'.format(key_variable, val_variable,
                                                                                 netcdf_file))
        return my_dict_data
    except IOError:
        raise Exception('File {} could not be found.'.format(netcdf_file))


def write_flow_file_from_nds(series_report, discharge, the_file, out_file_format, parallel=False):
    if out_file_format == 'netcdf':
        if Dataset:
            write_flow_netcdf_file_from_nds(series_report, discharge, the_file, parallel=parallel)
        else:
            raise Exception("The use of 'netcdf' as the output file format requires the package 'netCDF4', "
                            "please install it and retry, or choose another file format.")
    elif out_file_format == 'csv':
        write_flow_csv_file_from_nds(series_report, discharge, the_file)
    else:
        raise Exception("The output format type \'{}\' cannot be written by SMARTpy, "
                        "choose from: \'csv\', \'netcdf\'.".format(out_file_format))


def write_flow_csv_file_from_nds(series_report, discharge, csv_file):
    with open(csv_file, 'w', newline='', encoding='utf8') as my_file:
        my_writer = writer(my_file, delimiter=',')
        my_writer.writerow(['DateTime', 'flow'])
        for dt, val in zip(series_report, discharge):
            my_writer.writerow([dt, '%e' % val])


def write_flow_netcdf_file_from_nds(series_report, discharge, netcdf_file, parallel):
    with Dataset(netcdf_file + '.nc', 'w', format='NETCDF4', parallel=parallel) as my_file:
        my_file.description = "Discharge file generated with SMARTpy v{}.".format(__version__)
        my_file.createDimension('DateTime', len(series_report))
        t = my_file.createVariable("DateTime", np.float64, ('DateTime',))
        t.units = 'seconds since 1970-01-01 00:00:00.0'
        my_file.createVariable('flow', np.float32, ('DateTime',))

        my_file.variables['DateTime'][0:len(series_report)] = \
            (np.asarray(series_report, dtype='datetime64[us]') - np.datetime64('1970-01-01T00:00:00')) / \
            np.timedelta64(1, 's')
        my_file.variables['flow'][0:len(series_report)] = discharge


def valid_file_format(fmt):
    if fmt.lower() == "netcdf":
        if Dataset:  # i.e. netCDF4 import was successful
            return "netcdf"
        else:
            raise argparse.ArgumentTypeError("NetCDF4 module is not installed, please choose another file format.")
    elif fmt.lower() == "csv":
        return "csv"
    else:
        raise argparse.ArgumentTypeError("File format not recognised: '{0}'.".format(fmt))
