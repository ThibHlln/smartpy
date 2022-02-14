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

from builtins import range
from datetime import datetime, timedelta
import argparse
from collections import OrderedDict
from math import gcd


class TimeFrame(object):
    """
    This class gathers the temporal information provided by the user. Then, it defines the temporal attributes
    for the simulation.

    Two type of temporal attributes are considered: 'simu' refers to the internal time used by the model (i.e.
    model temporal discretisation), and 'save' refers to the output from the SMARTpy (i.e. what will be written in the
    output files).

    For each type, three attributes (start, end, gap) are required to build timeseries.

    In terms of starts/ends, the user is only required to provide them for 'save', start/end for 'simu'
    are inferred using 'save' start/end/gap and the 'simu' gap.

    In terms of timeseries, both 'simu' and 'save' timeseries are built using their respective start/end/gap.

    To allow for initial conditions to be set up, one or more time steps are added prior the respective starts of
    'simu' and 'save'. The 'save' time is added one datetime before save_start, while the 'simu' is added one or more
    datetime if it requires several simu_gap to cover one data_gap.

    N.B. The 'save' gap is required to be a multiple of 'simu' gap because this simulator is not intended to
    interpolate on the simulation results, the simulator can only extract or summarise from simulation steps.
    Instead, the user is expected to adapt their simulation gap to match the required reporting gap.
    """
    def __init__(self, dt_save_start, dt_save_end,
                 simu_increment, save_increment):

        # Time Attributed for Output Data (Save/Write in Files)
        self.save_start = dt_save_start  # DateTime
        self.save_gap = save_increment  # TimeDelta
        self.save_end = self._check_save_end(dt_save_end)  # DateTime

        # Time Attributes for Simulation (Internal to the Simulator)
        self.simu_gap = simu_increment  # TimeDelta
        self.simu_start, self.simu_end = \
            TimeFrame._get_simu_start_end_given_save_start_end(self)  # DateTime, DateTime

        # DateTime Series for Data, Save, and Simulation
        self.save_series = TimeFrame._get_list_save_dt_with_initial_conditions(self)
        self.simu_series = TimeFrame._get_list_simu_dt_with_initial_conditions(self)

    def _check_save_end(self, save_end):
        # check whether saving/reporting period makes sense on its own
        if not self.save_start <= save_end:
            raise Exception("Save Start is greater than Save End.")

        # determine whether the save_end is coherent with save_start and save_gap
        (divisor, remainder) = divmod(int((save_end - self.save_start).total_seconds()),
                                      int(self.save_gap.total_seconds()))

        if remainder != 0:
            save_end = self.save_start + timedelta(seconds=self.save_gap.total_seconds()) * divisor
            raise Exception("The combination of (start, end) datetimes and the saving time delta "
                            "are not compatible. For a start at {}, and a time delta of {}, the "
                            "latest end in the period is {}".format(self.save_start, self.save_gap, save_end))
        else:
            return save_end

    def _get_simu_start_end_given_save_start_end(self):
        # check whether the simulation time gap will allow to report on the saving/reporting time gap
        if not self.save_gap.total_seconds() % self.simu_gap.total_seconds() == 0:
            raise Exception("Save Gap is not greater and a multiple of Simulation Gap.")

        # determine the simulation period required to cover the whole saving/reporting period
        simu_start = self.save_start - self.save_gap + self.simu_gap
        simu_end = self.save_end

        return simu_start, simu_end

    def _get_list_save_dt_with_initial_conditions(self):

        # generate a list of DateTime for Saving/Reporting Period with one extra prior step for initial conditions
        my_list_datetime = list()
        my_dt = self.save_start - self.save_gap
        while my_dt <= self.save_end:
            my_list_datetime.append(my_dt)
            my_dt += self.save_gap

        return my_list_datetime

    def _get_list_simu_dt_with_initial_conditions(self):

        # generate a list of DateTime for Simulation Period with one extra prior step for initial conditions
        my_list_datetime = list()
        my_dt = self.simu_start - self.simu_gap
        while my_dt <= self.simu_end:
            my_list_datetime.append(my_dt)
            my_dt += self.simu_gap

        return my_list_datetime

    def get_gap_simu(self):
        return self.simu_gap

    def get_gap_report(self):
        return self.save_gap

    def get_series_simu(self):
        return self.simu_series

    def get_series_save(self):
        return self.save_series


def valid_date(s):
    try:
        return datetime.strptime(s, "%d/%m/%Y_%H:%M:%S")
    except ValueError:
        raise argparse.ArgumentTypeError("Not a valid date: '{0}'.".format(s))


def valid_delta_min(n):
    try:
        return timedelta(minutes=int(n))
    except ValueError:
        raise argparse.ArgumentTypeError("Not a valid time delta: '{0}'.".format(n))


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
    return timedelta(seconds=gcd(int((start_data - start_simu).total_seconds()),
                                 gcd(int(delta_data.total_seconds()), int(delta_simu.total_seconds()))))


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
        for my_sub_step in range(0, -divisor, -1):
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
        for my_sub_step in range(0, -divisor, -1):
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
    my_dts = list(dict_info)
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
        for my_sub_step in range(0, -divisor, -1):
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
            for my_sub_step in range(0, -divisor, -1):
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
