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

from os import path, makedirs, sep
import numpy as np
from builtins import dict

from .timeframe import TimeFrame
from .parameters import Parameters
from .inout import \
    get_dict_rain_series_simu, get_dict_peva_series_simu, get_dict_discharge_series, write_flow_file_from_dict
from .structure import run


class SMART(object):
    def __init__(self, catchment, c_area_m2, g_area_m2, start, end,
                 time_delta_simu, time_delta_save, warm_up_days, in_format, root):
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
        self.delta_save = time_delta_save
        self.timeframe = TimeFrame(self.start, self.end, self.delta_simu, self.delta_save)
        self.timeseries = self.timeframe.get_series_simu()
        self.timeseries_report = self.timeframe.get_series_save()
        self.warm_up = warm_up_days
        # physical information
        extra_ext = '.nc' if in_format == 'netcdf' else ''
        self.rain = get_dict_rain_series_simu(''.join([self.in_f, self.catchment, '.rain' + extra_ext]), in_format,
                                              self.timeseries[1], self.timeseries[-1], self.delta_simu)
        self.peva = get_dict_peva_series_simu(''.join([self.in_f, self.catchment, '.peva' + extra_ext]), in_format,
                                              self.timeseries[1], self.timeseries[-1], self.delta_simu)
        self.flow = get_dict_discharge_series(''.join([self.in_f, self.catchment, '.flow']),
                                              self.timeframe.get_series_save()[1],
                                              self.timeframe.get_series_save()[-1],
                                              c_area_m2, g_area_m2)
        # optional extra information for setting up initial levels in reservoirs
        self.extra = None
        # parameters
        self.parameters = Parameters()
        # model outputs
        self.outputs = None
        self.discharge = None
        self.gw_contribution = None

    def simulate(self, param):
        db = dict()
        self.outputs = run(self.area, self.delta_simu, self.rain, self.peva,
                           param, self.extra, db, self.timeseries, self.timeseries_report,
                           warm_up=self.warm_up)
        self.discharge = self.outputs[0]
        self.gw_contribution = self.outputs[1]
        return self.outputs

    def write_output_files(self):
        if self.outputs:
            write_flow_file_from_dict(self.timeframe, self.discharge,
                                      ''.join([self.out_f, self.catchment, '.mod.flow']),
                                      method='raw')

            write_flow_file_from_dict(self.timeframe, self.flow,
                                      ''.join([self.out_f, self.catchment, '.obs.flow']),
                                      method='raw')
        else:
            raise Exception("The output files cannot be written because the outputs attribute is unassigned."
                            "Please make sure to call the simulate method before writing the output files.")

    def get_simulation_array(self):
        if self.outputs:
            return np.asarray([self.discharge[dt] for dt in self.flow])
        else:
            raise Exception("The simulation array cannot be provided because the outputs attribute is unassigned."
                            "Please make sure to call the simulate method before requesting the simulation array.")

    def get_evaluation_array(self):
        return np.asarray([val for val in self.flow.values()])
