# -*- coding: utf-8 -*-

# This file is part of SMARTpy - An open-source rainfall-runoff model in Python
# Copyright (C) 2018  Thibault Hallouin (1)
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

from builtins import range, dict
from csv import DictReader, writer
from datetime import datetime, timedelta
from numpy import float64
from collections import OrderedDict
import sys
import io
import argparse
import imp
try:
    from netCDF4 import Dataset
except ImportError:
    Dataset = None

from .timeframe import get_required_resolution, check_interval_in_list, \
    rescale_time_resolution_of_irregular_mean_data, rescale_time_resolution_of_regular_cumulative_data


def open_csv_rb(my_file):
    if sys.version_info[0] < 3:
        return io.open(my_file, 'rb')
    else:
        return io.open(my_file, 'r', encoding='utf8')


def open_csv_wb(my_file):
    if sys.version_info[0] < 3:
        return io.open(my_file, 'wb')
    else:
        return io.open(my_file, 'w', newline='', encoding='utf8')


def open_csv_ab(my_file):
    if sys.version_info[0] < 3:
        return io.open(my_file, 'ab')
    else:
        return io.open(my_file, 'a', newline='', encoding='utf8')


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
                                                                  time_delta_simu,
                                                                  start_simu, end_simu, time_delta_res)
    else:
        raise Exception('PEva data not sufficient for simulation.')


def get_dict_discharge_series(file_location, start_report, end_report, catchment_area, gauged_area):
    data_flow = read_flow_file(file_location)

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
        raise Exception('Setting GAUGED AREA is missing from simulation file.')
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
        gw_constraint = -999.0  # i.e. no simulation constraint required
    except ValueError:
        raise Exception('Setting GROUNDWATER CONSTRAINT could not be converted to a float.')

    return c_area, g_area, start, end, delta_simu, delta_report, warm_up, gw_constraint


def read_rain_file(file_location, file_format):
    if file_format == 'netcdf':
        if Dataset:
            return read_netcdf_time_series_with_delta_check(file_location, key_variable='DateTime', val_variable='rain')
        else:
            raise Exception("The use of 'netcdf' as the output file format requires the package 'netCDF4', "
                            "please install it and retry, or choose another file format.")
    else:
        return read_csv_time_series_with_delta_check(file_location, key_header='DateTime', val_header='rain')


def read_peva_file(file_location, file_format):
    if file_format == 'netcdf':
        if Dataset:
            return read_netcdf_time_series_with_delta_check(file_location, key_variable='DateTime', val_variable='peva')
        else:
            raise Exception("The use of 'netcdf' as the output file format requires the package 'netCDF4', "
                            "please install it and retry, or choose another file format.")
    else:
        return read_csv_time_series_with_delta_check(file_location, key_header='DateTime', val_header='peva')


def read_flow_file(file_location):
    return read_csv_time_series_with_missing_check(file_location, key_header='DateTime', val_header='flow')


def read_simulation_settings_file(file_location):
    my_dict_args = dict()
    try:
        with open_csv_rb(file_location) as my_file:
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
        with open_csv_rb(csv_file) as my_file:
            my_dict_data = dict()
            my_list_dt = list()
            my_reader = DictReader(my_file)
            try:
                for row in my_reader:
                    my_dict_data[datetime.strptime(row[key_header], "%Y-%m-%d %H:%M:%S")] = float64(row[val_header])
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
                my_list_dt += [datetime.utcfromtimestamp(tstamp) for tstamp in my_file.variables[key_variable][:]]
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
        with open_csv_rb(csv_file) as my_file:
            my_dict_data = OrderedDict()
            my_reader = DictReader(my_file)
            try:
                for row in my_reader:
                    try:
                        if row[val_header] != '':  # not an empty string (that would mean missing data)
                            if float64(row[val_header]) != -99.0:  # flag for missing data
                                my_dict_data[datetime.strptime(row[key_header], "%Y-%m-%d %H:%M:%S")] = \
                                    float64(row[val_header])
                    except ValueError:
                        raise Exception('Field {} in {} cannot be converted to float '
                                        'at {}.'.format(val_header, csv_file, row[key_header]))
            except KeyError:
                raise Exception('Field {} or {} does not exist in {}.'.format(key_header, val_header, csv_file))
        return my_dict_data
    except IOError:
        raise Exception('File {} could not be found.'.format(csv_file))


def write_flow_file_from_list(timeframe, discharge, csv_file, report='save_gap', method='summary'):
    # Select the relevant list of DateTime given the argument used during function call
    if report == 'save_gap':  # standard situation
        my_list_datetime = timeframe.get_series_save()  # list of DateTime to be written in file
        simu_steps_per_reporting_step = \
            int(timeframe.get_gap_report().total_seconds() / timeframe.get_gap_simu().total_seconds())
    elif report == 'simu_gap':  # useful for debugging
        my_list_datetime = timeframe.get_series_simu()  # list of DateTime to be written in file
        simu_steps_per_reporting_step = 1
    else:
        raise Exception('Unknown reporting time gap for updating simulations files.')

    if method == 'summary':
        with open_csv_wb(csv_file) as my_file:
            my_writer = writer(my_file, delimiter=',')
            my_writer.writerow(['DateTime', 'flow'])
            my_index_simu = simu_steps_per_reporting_step   # ignoring first value that is for initial conditions
            my_index_report = 1  # ignoring first value that is for initial conditions
            while my_index_report <= len(my_list_datetime) - 1:
                my_values = list()
                for my_sub_index in range(0, -simu_steps_per_reporting_step, -1):
                    my_values.append(discharge[my_index_simu + my_sub_index])
                my_value = sum(my_values) / len(my_values)
                my_writer.writerow([my_list_datetime[my_index_report], '%e' % my_value])
                my_index_simu += simu_steps_per_reporting_step
                my_index_report += 1
    elif method == 'raw':
        with open_csv_wb(csv_file) as my_file:
            my_writer = writer(my_file, delimiter=',')
            my_writer.writerow(['DateTime', 'flow'])
            my_index_simu = simu_steps_per_reporting_step  # ignoring first value that is for initial conditions
            my_index_report = 1  # ignoring first value that is for initial conditions
            while my_index_report <= len(my_list_datetime):
                my_value = discharge[my_index_simu]
                my_writer.writerow([my_list_datetime[my_index_report], '%e' % my_value])
                my_index_simu += simu_steps_per_reporting_step
                my_index_report += 1
    else:
        raise Exception("Unknown method for updating simulations files.")


def write_flow_file_from_dict(timeframe, discharge, csv_file, report='save_gap', method='summary'):
    # Select the relevant list of DateTime given the argument used during function call
    if report == 'save_gap':  # standard situation
        my_list_datetime = timeframe.get_series_save()  # list of DateTime to be written in file
        simu_steps_per_reporting_step = \
            int(timeframe.get_gap_report().total_seconds() / timeframe.get_gap_simu().total_seconds())
    elif report == 'simu_gap':  # useful for debugging
        my_list_datetime = timeframe.get_series_simu()  # list of DateTime to be written in file
        simu_steps_per_reporting_step = 1
    else:
        raise Exception('Unknown reporting time gap for updating simulations files.')

    if method == 'summary':
        with open_csv_wb(csv_file) as my_file:
            my_writer = writer(my_file, delimiter=',')
            my_writer.writerow(['DateTime', 'flow'])
            for step in my_list_datetime[1:]:
                my_values = list()
                for my_sub_step in range(0, -simu_steps_per_reporting_step, -1):
                    my_values.append(
                        discharge[step + my_sub_step * timeframe.gap_simu])
                my_value = sum(my_values) / len(my_values)
                my_writer.writerow([step, '%e' % my_value])
    elif method == 'raw':
        with open_csv_wb(csv_file) as my_file:
            my_writer = writer(my_file, delimiter=',')
            my_writer.writerow(['DateTime', 'flow'])
            for step in my_list_datetime[1:]:
                try:
                    my_writer.writerow([step, '%e' % discharge[step]])
                except KeyError:
                    my_writer.writerow([step, ''])
    else:
        raise Exception("Unknown method for updating simulations files.")


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
