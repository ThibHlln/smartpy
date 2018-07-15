from fractions import gcd
from datetime import timedelta
from collections import OrderedDict


class TimeFrame(object):
    """
    This class defines the temporal attributes of the simulation period. It contains the start and the end of the
    simulation as well as the lists of DateTime series for the simulation time steps and the reporting time steps (that
    can be identical or nested).

    N.B. 1: The simulation gap needs to be a multiple if the reporting gap and the simulation gap can ony
    be lower than or equal to the report gap
    N.B. 2: The start and the end of the simulation are defined by the user, the class always adds one data step
    prior to the start date in order to set the initial conditions, one or more simulation steps are added in
    consequence depending if the simulation step in a multiple of the data step or not (i.e. equal)
    """
    def __init__(self, datetime_start, datetime_end, simu_timedelta, report_timedelta):
        assert datetime_start <= datetime_end, "TimeFrame: Start > End"
        assert report_timedelta.total_seconds() % simu_timedelta.total_seconds() == 0, \
            "Reporting TimeDelta is not a multiple of Simulation TimeDelta."
        assert (datetime_end - datetime_start).total_seconds() % simu_timedelta.total_seconds() == 0, \
            "Simulation Period is not a multiple of Simulation TimeDelta."
        # DateTime of the start of the time period simulated
        self.start = datetime_start
        # DateTime of the end of the time period simulated
        self.end = datetime_end
        # TimeDelta of the simulation
        self.gap_simu = simu_timedelta
        # TimeDelta of the reporting
        self.gap_report = report_timedelta
        # List of DateTime for the reporting (i.e. list of time steps)
        self.series_report = TimeFrame._create_list_datetime(self, 'report')
        # List of DateTime for the simulation (i.e. list of time steps)
        self.series_simu = TimeFrame._create_list_datetime(self, 'simu')

    def _create_list_datetime(self, option):
        """
        This function returns a list of DateTime by using the start and the end of the simulation and the time gap
        (either the reporting time gap or the simulation time gap, using the option parameter to specify which one).

        N.B. For the initial conditions, the function always adds:
            - [if 'report' option] one data step prior to the reporting start date
            - [if 'simu' option] one (or more if reporting gap > simulation gap) simulation step(s)
            prior to the simulation start date

        :param option: choice to specify if function should work on reporting or on simulation series
        :type option: str()
        :return: a list of DateTime
        :rtype: list()
        """
        extent = self.end - self.start
        options = {'report': self.gap_report, 'simu': self.gap_simu}

        start_index = int((self.gap_report.total_seconds()) / (options[option].total_seconds()))
        end_index = int(extent.total_seconds() // (options[option].total_seconds())) + 1

        my_list_datetime = list()
        for factor in xrange(-start_index, end_index, 1):  # add one or more datetime before start
            my_datetime = self.start + factor * options[option]
            my_list_datetime.append(my_datetime)

        return my_list_datetime

    def get_gap_simu(self):
        return self.gap_simu

    def get_gap_report(self):
        return self.gap_report

    def get_series_simu(self):
        return self.series_simu

    def get_series_report(self):
        return self.series_report


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
            my_value = float('nan')
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
            new_dict_info[my_dt] = float('nan')
            my_dt += time_delta_lo
        except TypeError:  # at least one of the values was not a float [string] (i.e. missing value)
            new_dict_info[my_dt] = float('nan')
            my_dt += time_delta_lo

    return new_dict_info


def rescale_time_resolution_of_irregular_mean_data(dict_data, start_data, end_data, time_delta_lo, time_delta_hi):
    my_tmp_dict = increase_time_resolution_of_irregular_mean_data(dict_data, time_delta_lo, time_delta_hi)
    my_new_dict = decrease_time_resolution_of_irregular_mean_data(my_tmp_dict, start_data, end_data,
                                                                  time_delta_hi, time_delta_lo)

    return my_new_dict
