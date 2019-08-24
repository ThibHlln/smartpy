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

from .timeframe import TimeFrame
from .parameters import Parameters
from .inout import \
    get_dict_rain_series_simu, get_dict_peva_series_simu, get_dict_discharge_series, write_flow_file_from_nds
from .structure import run


class SMART(object):
    def __init__(self, catchment, catchment_area_m2, start, end,
                 time_delta_simu, time_delta_save, warm_up_days,
                 in_format, out_format, root,
                 gauged_area_m2=None):
        # general information
        self.catchment = catchment
        self.area = catchment_area_m2
        # directory information
        self.in_fmt = in_format
        self.out_fmt = out_format
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
        # physical information as dictionaries
        extra_ext = '.nc' if self.in_fmt == 'netcdf' else ''
        self.rain = get_dict_rain_series_simu(''.join([self.in_f, self.catchment, '.rain' + extra_ext]), self.in_fmt,
                                              self.timeseries[1], self.timeseries[-1], self.delta_simu)
        self.peva = get_dict_peva_series_simu(''.join([self.in_f, self.catchment, '.peva' + extra_ext]), self.in_fmt,
                                              self.timeseries[1], self.timeseries[-1], self.delta_simu)
        self.flow = get_dict_discharge_series(''.join([self.in_f, self.catchment, '.flow' + extra_ext]),
                                              self.in_fmt,
                                              self.timeframe.get_series_save()[1],
                                              self.timeframe.get_series_save()[-1],
                                              catchment_area_m2, gauged_area_m2) if gauged_area_m2 else None
        # physical information as numpy arrays
        self.nd_rain = np.array([self.rain[dt] for dt in self.timeseries[1:]])
        self.nd_peva = np.array([self.peva[dt] for dt in self.timeseries[1:]])
        self.nd_flow = np.array([self.flow[dt] for dt in self.timeseries_report[1:]]) if gauged_area_m2 else None
        # optional extra information for setting up initial levels in reservoirs
        self.extra = None
        # parameters
        self.parameters = Parameters()
        # model outputs
        self.outputs = None
        self.nd_discharge = None
        self.gw_contribution = None

    def simulate(self, param, report='summary'):
        nd_parameters = np.array([param[name] for name in self.parameters.names])
        self.outputs = run(self.area, self.delta_simu, self.nd_rain, self.nd_peva,
                           nd_parameters, self.extra, self.timeseries, self.timeseries_report,
                           report=report, warm_up=self.warm_up)
        self.nd_discharge = self.outputs[0]
        self.gw_contribution = self.outputs[1]
        return self.outputs

    def write_output_files(self, which='both', parallel=False):
        if (which == 'both') or (which == 'modelled'):
            if self.nd_discharge is not None:
                write_flow_file_from_nds(self.timeseries_report[1:], self.nd_discharge,
                                         ''.join([self.out_f, self.catchment, '.mod.flow']),
                                         out_file_format=self.out_fmt, parallel=parallel)
            else:
                raise Exception("The modelled flow output file cannot be written. Please make sure to call the "
                                "simulate method of your SMART instance before writing this output file.")

        if (which == 'both') or (which == 'observed'):
            if self.nd_flow is not None:
                write_flow_file_from_nds(self.timeseries_report[1:], self.nd_flow,
                                         ''.join([self.out_f, self.catchment, '.obs.flow']),
                                         out_file_format=self.out_fmt, parallel=parallel)
            else:
                raise Exception("The observed flow output file cannot be written. Please make sure that a value is "
                                "assigned to the gauged_area_m2 attribute of the SMART class instance.")

    def get_simulation_array(self):
        if self.nd_discharge is not None:
            return self.nd_discharge
        else:
            raise Exception("The simulation array cannot be retrieved. Please make sure to call the simulate "
                            "method of your SMART instance before requesting this output array.")

    def get_evaluation_array(self):
        if self.nd_flow is not None:
            return self.nd_flow
        else:
            raise Exception("The observation array does not exist. Please make sure that a value is assigned "
                            "to the gauged_area_m2 attribute of your SMART class instance.")
