from csv import DictReader, writer
from datetime import datetime, timedelta
from fractions import gcd
from numpy import float64
from collections import OrderedDict


def get_dict_rain_series_simu(file_location, start_simu, end_simu, time_delta_simu):
    dict_rain, start_data, end_data, time_delta_data = read_rain_file(file_location)

    if (start_data - time_delta_data + time_delta_simu <= start_simu) and (end_simu <= end_data):
        time_delta_res = get_required_resolution(start_data, start_simu, time_delta_data, time_delta_simu)
        return rescale_time_resolution_of_regular_cumulative_data(dict_rain,
                                                                  start_data, end_data, time_delta_data,
                                                                  time_delta_res,
                                                                  start_simu, end_simu, time_delta_simu)
    else:
        raise Exception('Rain data not sufficient for simulation.')


def get_dict_peva_series_simu(file_location, start_simu, end_simu, time_delta_simu):
    dict_peva, start_data, end_data, time_delta_data = read_peva_file(file_location)

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
    for dt in data_flow.iterkeys():
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


def read_rain_file(file_location):
    return read_csv_time_series_with_delta_check(file_location, key_header='DATETIME', val_header='RAIN')


def read_peva_file(file_location):
    return read_csv_time_series_with_delta_check(file_location, key_header='DATETIME', val_header='PEVA')


def read_flow_file(file_location):
    return read_csv_time_series_with_missing_check(file_location, key_header='DATETIME', val_header='FLOW')


def read_simulation_settings_file(file_location):
    my_dict_args = dict()
    try:
        with open(file_location, 'rb') as my_file:
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
        with open(csv_file, 'rb') as my_file:
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


def read_csv_time_series_with_missing_check(csv_file, key_header, val_header):
    try:
        with open(csv_file, 'rb') as my_file:
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


def check_interval_in_list(list_of_dt, csv_file):
    list_intervals = list()
    for i in range(len(list_of_dt) - 1):
        list_intervals.append(list_of_dt[i+1] - list_of_dt[i])
    interval = list(set(list_intervals))
    if len(interval) == 1:
        if list_of_dt[0] + interval[0] * (len(list_of_dt) - 1) == list_of_dt[-1]:
            return list_of_dt[0], list_of_dt[-1], interval[0]
        else:
            raise Exception('Missing Data: {} is missing at least one datetime in period.'.format(csv_file))
    else:
        raise Exception('Inconsistent Interval: {} does not feature a single time interval.'.format(csv_file))


def get_required_resolution(start_data, start_simu, delta_data, delta_simu):
    # GCD(delta_data, delta_simu) gives the maximum time resolution possible to match data and simu
    # shift = start_data - start_simu gives the data shift (e.g. data starting at 8am, simu starting at 9am)
    # GCD(shift, GCD(delta_data, delta_simu)) gives the maximum time resolution to match both the difference in
    # start dates and the difference in data/simu time deltas.
    return timedelta(seconds=gcd((start_data - start_simu).total_seconds(),
                                 gcd(delta_data.total_seconds(), delta_simu.total_seconds())))


def increase_time_resolution_of_regular_cumulative_data(dict_info, start_lo, end_lo,
                                                        time_delta_lo, time_delta_hi):
    """ Use the low resolution to create the high resolution """
    my_dt_lo = start_lo
    (divisor, remainder) = divmod(int(time_delta_lo.total_seconds()), int(time_delta_hi.total_seconds()))
    if remainder != 0:
        raise Exception("Increase Resolution: Time Deltas are not multiples of each other.")
    elif divisor < 1:
        raise Exception("Increase Resolution: Low resolution lower than higher resolution "
                        "{} < {}.".format(time_delta_lo, time_delta_hi))

    new_dict_info = dict()
    while (start_lo <= my_dt_lo) and (my_dt_lo <= end_lo):
        my_value = dict_info[my_dt_lo]
        my_portion = my_value / divisor
        for my_sub_step in xrange(0, -divisor, -1):
            new_dict_info[my_dt_lo + my_sub_step * time_delta_hi] = my_portion
        my_dt_lo += time_delta_lo

    return new_dict_info


def decrease_time_resolution_of_regular_cumulative_data(dict_info, start_lo, end_lo,
                                                        time_delta_lo, time_delta_hi):
    """ Use the high resolution to create the low resolution """
    my_dt_lo = start_lo
    (divisor, remainder) = divmod(int(time_delta_lo.total_seconds()), int(time_delta_hi.total_seconds()))
    if remainder != 0:
        raise Exception("Decrease Resolution: Time Deltas are not multiples of each other.")
    elif divisor < 1:
        raise Exception("Decrease Resolution: Low resolution lower than higher resolution "
                        "{} < {}.".format(time_delta_lo, time_delta_hi))

    new_dict_info = dict()
    while (start_lo <= my_dt_lo) and (my_dt_lo <= end_lo):
        my_portion = 0.0
        for my_sub_step in xrange(0, -divisor, -1):
            my_portion += dict_info[my_dt_lo + my_sub_step * time_delta_hi]
        new_dict_info[my_dt_lo] = my_portion
        my_dt_lo += time_delta_lo

    return new_dict_info


def rescale_time_resolution_of_regular_cumulative_data(dict_data,
                                                       start_data, end_data, time_delta_data,
                                                       time_delta_res,
                                                       start_simu, end_simu, time_delta_simu):

    if time_delta_data > time_delta_res:  # i.e. information resolution too low to generate simu timeseries
        my_tmp_dict = increase_time_resolution_of_regular_cumulative_data(dict_data, start_data, end_data,
                                                                          time_delta_data, time_delta_res)
    else:  # i.e. information resolution suitable to generate simu timeseries
        # i.e. time_delta_data == time_delta_res (time_delta_data < time_delta_res cannot be true because use of GCD)
        my_tmp_dict = dict_data

    if time_delta_simu > time_delta_res:  # i.e. information resolution too high to generate simu timeseries
        my_new_dict = decrease_time_resolution_of_regular_cumulative_data(my_tmp_dict, start_simu, end_simu,
                                                                          time_delta_simu, time_delta_res)
    else:  # i.e. information resolution suitable to generate simu timeseries
        # i.e. time_delta_simu == time_delta_res (time_delta_simu < time_delta_res cannot be true because use of GCD)
        my_new_dict = decrease_time_resolution_of_regular_cumulative_data(my_tmp_dict, start_simu, end_simu,
                                                                          time_delta_simu, time_delta_res)
        # use decrease_data_time_resolution anyway for the only purpose to reduce the size of the dict
        # to the only required DateTimes in the simulation period

    return my_new_dict


def increase_time_resolution_of_irregular_mean_data(dict_info, time_delta_lo, time_delta_hi):
    """
    Create high resolution mean data from lower resolution mean data
    using backwards replication.
    """
    new_dict_info = dict()
    # get the series of DateTime in the data
    my_dts = dict_info.keys()
    # special case for first datetime that does not have an antecedent value
    my_dts.insert(0, my_dts[0] - time_delta_lo)  # add one virtual date before

    for i, my_dt in enumerate(my_dts[1:]):
        # determine the duration of the cumulative data between the two time steps
        my_delta = my_dt - my_dts[i]
        if my_delta >= timedelta(seconds=1.5 * time_delta_lo.total_seconds()):  # i.e. gap is likely > one missing data
            my_delta = time_delta_lo  # so set delta arbitrarily to the delta of the low resolution data ('standard')
        # determine the number of hours the cumulative data has to be spread over
        (divisor, remainder) = divmod(int(my_delta.total_seconds()), int(time_delta_hi.total_seconds()))
        if remainder != 0:
            raise Exception("Increase Resolution: Time Deltas are not multiples of each other.")
        elif divisor < 1:
            raise Exception("Increase Resolution: Low resolution lower than higher resolution "
                            "{} < {}.".format(my_delta, time_delta_hi))
        # determine the high resolution value
        try:
            my_value = float(dict_info[my_dt])
        except ValueError:  # no data for this time step
            my_value = ''
        # spread the hourly value backwards
        for my_sub_step in xrange(0, -divisor, -1):
            if new_dict_info.get(my_dt + my_sub_step * time_delta_hi):  # should not happen
                raise Exception("Increase Resolution: Overwriting already existing data for datetime.")
            else:
                new_dict_info[my_dt + my_sub_step * time_delta_hi] = my_value

    return new_dict_info


def decrease_time_resolution_of_irregular_mean_data(dict_info, dt_start, dt_end, time_delta_hi, time_delta_lo):
    """ Creates low resolution cumulative data from high resolution cumulative data
    using arithmetic mean. """
    my_dt = dt_start
    new_dict_info = OrderedDict()

    (divisor, remainder) = divmod(int(time_delta_lo.total_seconds()), int(time_delta_hi.total_seconds()))
    if remainder != 0:
        raise Exception("Decrease Resolution: Time Deltas are not multiples of each other.")
    elif divisor < 1:
        raise Exception("Decrease Resolution: Low resolution lower than higher resolution "
                        "{} < {}.".format(time_delta_lo, time_delta_hi))

    while (dt_start <= my_dt) and (my_dt <= dt_end):
        try:
            my_values = 0.0
            for my_sub_step in xrange(0, -divisor, -1):
                my_values += dict_info[my_dt + my_sub_step * time_delta_hi]
            new_dict_info[my_dt] = my_values / divisor
            my_dt += time_delta_lo
        except KeyError:  # at least one of the values is not available (i.e. missing value)
            my_dt += time_delta_lo
        except TypeError:  # at least one of the values was not a float [string] (i.e. missing value)
            my_dt += time_delta_lo

    return new_dict_info


def rescale_time_resolution_of_irregular_mean_data(dict_data, start_data, end_data, time_delta_lo, time_delta_hi):
    my_tmp_dict = increase_time_resolution_of_irregular_mean_data(dict_data, time_delta_lo, time_delta_hi)
    my_new_dict = decrease_time_resolution_of_irregular_mean_data(my_tmp_dict, start_data, end_data,
                                                                  time_delta_hi, time_delta_lo)

    return my_new_dict


def write_flow_file_from_list(timeframe, discharge, csv_file, report='gap_report', method='summary'):
    # Select the relevant list of DateTime given the argument used during function call
    if report == 'gap_report':  # standard situation
        my_list_datetime = timeframe.get_series_report()  # list of DateTime to be written in file
        simu_steps_per_reporting_step = \
            int(timeframe.get_gap_report().total_seconds() / timeframe.get_gap_simu().total_seconds())
    elif report == 'gap_simu':  # useful for debugging
        my_list_datetime = timeframe.get_series_simu()  # list of DateTime to be written in file
        simu_steps_per_reporting_step = 1
    else:
        raise Exception('Unknown reporting time gap for updating simulations files.')

    if method == 'summary':
        with open(csv_file, 'wb') as my_file:
            my_writer = writer(my_file, delimiter=',')
            my_writer.writerow(['DATETIME', 'FLOW'])
            my_index_simu = simu_steps_per_reporting_step   # ignoring first value that is for initial conditions
            my_index_report = 1  # ignoring first value that is for initial conditions
            while my_index_report <= len(my_list_datetime) - 1:
                my_values = list()
                for my_sub_index in xrange(0, -simu_steps_per_reporting_step, -1):
                    my_values.append(discharge[my_index_simu + my_sub_index])
                my_value = sum(my_values) / len(my_values)
                my_writer.writerow([my_list_datetime[my_index_report], '%e' % my_value])
                my_index_simu += simu_steps_per_reporting_step
                my_index_report += 1
    elif method == 'raw':
        with open(csv_file, 'wb') as my_file:
            my_writer = writer(my_file, delimiter=',')
            my_writer.writerow(['DATETIME', 'FLOW'])
            my_index_simu = simu_steps_per_reporting_step  # ignoring first value that is for initial conditions
            my_index_report = 1  # ignoring first value that is for initial conditions
            while my_index_report <= len(my_list_datetime):
                my_value = discharge[my_index_simu]
                my_writer.writerow([my_list_datetime[my_index_report], '%e' % my_value])
                my_index_simu += simu_steps_per_reporting_step
                my_index_report += 1
    else:
        raise Exception("Unknown method for updating simulations files.")


def write_flow_file_from_dict(timeframe, discharge, csv_file, report='gap_report', method='summary'):
    # Select the relevant list of DateTime given the argument used during function call
    if report == 'gap_report':  # standard situation
        my_list_datetime = timeframe.get_series_report()  # list of DateTime to be written in file
        simu_steps_per_reporting_step = \
            int(timeframe.get_gap_report().total_seconds() / timeframe.get_gap_simu().total_seconds())
    elif report == 'gap_simu':  # useful for debugging
        my_list_datetime = timeframe.get_series_simu()  # list of DateTime to be written in file
        simu_steps_per_reporting_step = 1
    else:
        raise Exception('Unknown reporting time gap for updating simulations files.')

    if method == 'summary':
        with open(csv_file, 'wb') as my_file:
            my_writer = writer(my_file, delimiter=',')
            my_writer.writerow(['DATETIME', 'FLOW'])
            for step in my_list_datetime[1:]:
                my_values = list()
                for my_sub_step in xrange(0, -simu_steps_per_reporting_step, -1):
                    my_values.append(
                        discharge[step + my_sub_step * timeframe.gap_simu])
                my_value = sum(my_values) / len(my_values)
                my_writer.writerow([step, '%e' % my_value])
    elif method == 'raw':
        with open(csv_file, 'wb') as my_file:
            my_writer = writer(my_file, delimiter=',')
            my_writer.writerow(['DATETIME', 'FLOW'])
            for step in my_list_datetime[1:]:
                try:
                    my_writer.writerow([step, '%e' % discharge[step]])
                except KeyError:
                    my_writer.writerow([step, ''])
    else:
        raise Exception("Unknown method for updating simulations files.")
