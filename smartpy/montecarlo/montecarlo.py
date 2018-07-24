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

from builtins import map
from csv import DictReader
import numpy as np
from os import sep, remove
from io import open
import gzip
import shutil

try:
    import spotpy
except ImportError:
    raise Exception('montecarlo requires the package spotpy to be installed.')

try:
    from netCDF4 import Dataset
except ImportError:
    Dataset = None

from ..smart import SMART
from ..inout import get_dict_simulation_settings, open_csv_rb, open_csv_wb
from ..objfunctions import \
    groundwater_constraint, bounded_nash_sutcliffe, sqrt_nash_sutcliffe, spearman_rank_corr, mean_abs_rel_error


class MonteCarlo(object):
    def __init__(self, catchment, root_f, in_format, out_format,
                 parallel, save_sim, func, settings_filename):
        in_f = sep.join([root_f, 'in', catchment, sep])

        # collect the simulation information from the .sttngs file
        if settings_filename:
            c_area, g_area, start, end, delta_simu, delta_report, warm_up, gw_constraint = \
                get_dict_simulation_settings(''.join([in_f, settings_filename]))
        else:
            c_area, g_area, start, end, delta_simu, delta_report, warm_up, gw_constraint = \
                get_dict_simulation_settings(''.join([in_f, catchment, '.sttngs']))

        # generate an instance of the SMART model class
        self.model = SMART(catchment, c_area, start, end, delta_simu, delta_report, warm_up, in_format, root_f, g_area)

        # set the technical aspects of the simulation
        self.parallel = parallel  # using mpi to run on several cores
        self.p = True if parallel == 'mpi' else False
        self.save_sim = save_sim  # saving simulation timeseries alongside objective functions and parameter values

        # set the possible additional constraint for SMART on the portion of base flow in runoff
        self.constraints = {'gw': gw_constraint}

        # collect list of parameters and objective functions names that are common for all Monte Carlo classes
        self.param_names = self.model.parameters.names
        self.obj_fn_names = \
            ['NSE', 'lgNSE', 'rtNSE', 'C2M', 'KGE', 'KGEc', 'KGEa', 'KGEb',
             'Bias', 'PBias', 'RMSE', 'Rho', 'MARE'] \
            if self.constraints['gw'] == -999.0 else \
            ['NSE', 'lgNSE', 'rtNSE', 'C2M', 'KGE', 'KGEc', 'KGEa', 'KGEb',
             'Bias', 'PBias', 'RMSE', 'Rho', 'MARE', 'GW']

        # define attributes to be used for sample of model parameters
        self.p_map = None
        self.params = None

        # define attributes to be used for the database to save results
        self.out_format = out_format
        self.db_file = \
            self.model.out_f + '{}.SMART.{}.nc'.format(catchment, func) if self.out_format == 'netcdf' else \
            self.model.out_f + '{}.SMART.{}'.format(catchment, func)
        self.database = None

    def _init_db(self, compression=6):
        if self.out_format == 'netcdf':
            if Dataset:  # check that netCDF4 is installed and imported
                # create structure of NetCDF file
                with Dataset(self.db_file, 'w', parallel=self.p) as my_file:
                    # metadata
                    my_file.description = "Monte Carlo Simulation outputs with SMARTpy."
                    # dimensions
                    my_file.createDimension("NbSamples", len(self.p_map))
                    my_file.createDimension("NbParameters", len(self.model.parameters.names))
                    my_file.createDimension("NbObjFunctions", len(self.obj_fn_names))
                    # variables
                    params = my_file.createVariable("Parameters", np.float32, ("NbSamples", "NbParameters"),
                                                    zlib=True, complevel=compression)
                    params.units = ', '.join(self.model.parameters.names)
                    objfns = my_file.createVariable("ObjFunctions", np.float32, ("NbSamples", "NbObjFunctions"),
                                                    zlib=True, complevel=compression)
                    objfns.units = ', '.join(self.obj_fn_names)
                    if self.save_sim:
                        # dimension
                        my_file.createDimension("DateTime", None)  # Unlimited dimension
                        # coordinate variables
                        times = my_file.createVariable("DateTime", np.float64, ("DateTime",), zlib=True)
                        times.units = 'seconds since 1970-01-01 00:00:00.0'
                        # variable
                        simu = my_file.createVariable("Simulations", np.float32, ("NbSamples", "DateTime"),
                                                      zlib=True, complevel=compression)
                        simu.units = 'Discharge in m3/s'
                        # write simulation datetime series as Unix timestamp series
                        datetimes = [np.datetime64(dt) for dt in self.model.flow]
                        timestamps = (datetimes - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
                        my_file.variables['DateTime'][0:len(datetimes)] = timestamps
            else:
                raise Exception("The use of 'netcdf' as the output file format requires the package 'netCDF4', "
                                "please install it and retry, or choose another file format.")

        else:  # fall back on to default option (CSV file)
            self.database = open_csv_wb(self.db_file)
            simu_steps = [dt.strftime("%Y-%m-%d %H:%M:%S") for dt in self.model.flow] if self.save_sim else []
            # write header in database file
            self.database.write(','.join(self.obj_fn_names + self.param_names + simu_steps) + '\n')

    def parameters(self):
        return spotpy.parameter.generate(self.params)

    def run(self, compression=None):
        # if compression specified, NetCDF4 needs to know compression level (between 1 and 9) when creating the file
        if self.out_format == 'netcdf' and not isinstance(compression, bool) \
                and isinstance(compression, (int, long, float)):
            self._init_db(compression=compression)
        else:
            self._init_db()
        # run the Monte Carlo simulation
        sampler = spotpy.algorithms.mc(self, parallel=self.parallel)
        sampler.sample(len(self.p_map))
        # if compression specified, the CSV file created will be compressed
        if self.out_format == 'csv':
            self.database.close()
            if compression is True:
                with open(self.db_file, 'rb') as f_in:
                    with gzip.open(self.db_file + '.gz', 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                remove(self.db_file)

    def simulation(self, vector):
        simulations, constraint = self.model.simulate({
            'T': vector[0], 'C': vector[1], 'H': vector[2], 'D': vector[3], 'S': vector[4], 'Z': vector[5],
            'SK': vector[6], 'FK': vector[7], 'GK': vector[8], 'RK': vector[9]
        })
        # returns only simulations with observations
        return [
            [simulations[dt] for dt in self.model.flow],
            [constraint]
        ]

    def evaluation(self):
        return [
            np.asarray([observation for observation in self.model.flow.values()]),
            [self.constraints['gw']]
        ]

    def objectivefunction(self, simulation, evaluation):
        # select the series subset that have observations (i.e. not NaN)
        flow_eval = evaluation[0][~np.isnan(evaluation[0])]
        flow_simu = np.asarray(simulation[0])[~np.isnan(evaluation[0])]

        # calculate the objective functions
        obj1 = spotpy.objectivefunctions.nashsutcliffe(evaluation=flow_eval, simulation=flow_simu)
        obj2 = spotpy.objectivefunctions.lognashsutcliffe(evaluation=flow_eval, simulation=flow_simu)
        obj3 = sqrt_nash_sutcliffe(evaluation=flow_eval, simulation=flow_simu)
        obj4 = bounded_nash_sutcliffe(evaluation=flow_eval, simulation=flow_simu)
        obj5, obj5c, obj5a, obj5b = \
            spotpy.objectivefunctions.kge(evaluation=flow_eval, simulation=flow_simu, return_all=True)
        obj6 = spotpy.objectivefunctions.bias(evaluation=flow_eval, simulation=flow_simu)
        obj7 = spotpy.objectivefunctions.pbias(evaluation=flow_eval, simulation=flow_simu)
        obj8 = spotpy.objectivefunctions.rmse(evaluation=flow_eval, simulation=flow_simu)
        obj9 = spearman_rank_corr(evaluation=flow_eval, simulation=flow_simu)
        obj10 = mean_abs_rel_error(evaluation=flow_eval, simulation=flow_simu)
        obj11 = groundwater_constraint(evaluation=evaluation[1], simulation=simulation[1])

        if self.constraints['gw'] == -999.0:
            return [obj1, obj2, obj3, obj4, obj5, obj5c, obj5a, obj5b, obj6, obj7, obj8, obj9, obj10]
        else:
            return [obj1, obj2, obj3, obj4, obj5, obj5c, obj5a, obj5b, obj6, obj7, obj8, obj9, obj10, obj11]

    def save(self, obj_fns, parameters, simulations, *args, **kwargs):
        if self.out_format == 'netcdf':
            # it was already checked if netCDF4 was installed and imported, so just proceed
            with Dataset(self.db_file, 'a', parallel=self.p) as my_file:
                # convert the parameter values array to a list
                params = parameters.tolist()
                # get the index of the particular parameter set to write data in the right spot
                index = self.p_map[tuple(params)]
                # write data on parameter values and objective functions
                my_file.variables['Parameters'][index, 0:len(self.param_names)] = params
                my_file.variables['ObjFunctions'][index, 0:len(self.obj_fn_names)] = obj_fns
                if self.save_sim:
                    # write simulation discharge time series
                    my_s, my_e = 0, len(my_file.variables['DateTime'])
                    my_file.variables['Simulations'][index, my_s:my_e] = simulations[0]
        else:
            if self.save_sim:  # create a list of objectives functions + parameters + discharge series
                line = map(np.float32, obj_fns + parameters.tolist() + simulations[0])
            else:  # create a list of objectives functions + parameters
                line = map(np.float32, obj_fns + parameters.tolist())
                # write data in file
            self.database.write(','.join(map(lambda x: '%.6e' % x, line)) + '\n')

    def _get_sampled_sets_from_file(self, file_location, param_names, obj_fn_names, decompression_csv):
        obj_fns, params = list(), list()
        if self.out_format == 'netcdf':
            # collect parameter values and objective function values from netCDF
            if Dataset:
                with Dataset(file_location, 'r') as my_file:
                    for i in range(my_file.variables['Parameters'].shape[0]):
                        params.append(my_file.variables['Parameters'][i, :].tolist())
                        obj_fns.append(my_file.variables['ObjFunctions'][i, :].tolist())
            else:
                raise Exception("The use of 'netcdf' as the output file format requires the package 'netCDF4', "
                                "please install it and retry, or choose another file format.")
        else:
            if decompression_csv:
                with gzip.open(file_location + '.gz', 'rb') as my_file:
                    my_reader = DictReader(my_file)
                    obj_fns, params = list(), list()
                    for row in my_reader:
                        obj_fns.append([row[obj_fn] for obj_fn in obj_fn_names])
                        params.append([row[param] for param in param_names])
            # collect parameter values and objective function values from CSV
            else:
                with open_csv_rb(file_location) as my_file:
                    my_reader = DictReader(my_file)
                    obj_fns, params = list(), list()
                    for row in my_reader:
                        obj_fns.append([row[obj_fn] for obj_fn in obj_fn_names])
                        params.append([row[param] for param in param_names])
        # return parameter values and objective function values in two separate arrays
        return np.array(params, dtype=np.float32), np.array(obj_fns, dtype=np.float32)
