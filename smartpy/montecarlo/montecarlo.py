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

from builtins import map
from csv import DictReader
import numpy as np
from os import sep, remove, rename
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
from ..inout import get_dict_simulation_settings
from ..objfunctions import groundwater_constraint
from ..version import __version__


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
        self.model = SMART(catchment, c_area, start, end, delta_simu, delta_report, warm_up,
                           in_format, out_format, root_f,
                           g_area)

        # set the technical aspects of the simulation
        self.parallel = parallel  # using mpi to run on several cores
        self.p = True if parallel == 'mpi' else False
        self.save_sim = save_sim  # saving simulation timeseries alongside objective functions and parameter values

        # set the possible additional constraint for SMART on the portion of base flow in runoff
        self.constraints = {'gw': gw_constraint}

        # collect list of parameters and objective functions names that are common for all Monte Carlo classes
        self.param_names = self.model.parameters.names
        self.obj_fn_names = \
            ['NSE', 'KGE', 'KGEc', 'KGEa', 'KGEb', 'PBias', 'RMSE', 'GW'] \
            if self.constraints['gw'] else \
            ['NSE', 'KGE', 'KGEc', 'KGEa', 'KGEb', 'PBias', 'RMSE']

        # define attributes to be used for sample of model parameters
        self.p_map = None
        self.params = None

        # define attributes to be used for the database to save results
        self.out_format = out_format
        self.db_file = \
            self.model.out_f + '{}.SMART.{}.nc'.format(catchment, func) if self.out_format == 'netcdf' else \
            self.model.out_f + '{}.SMART.{}'.format(catchment, func)
        self.database = None

        # write out the observed discharge data used for the objective functions (that might have been shifted/rescaled)
        self.model.write_output_files(which='observed', parallel=self.p)

    def _init_db(self):
        if self.out_format == 'netcdf':
            if Dataset:  # check that netCDF4 is installed and imported
                self.database = Dataset(self.db_file, 'w', format='NETCDF4', parallel=self.p)
                # create structure of NetCDF file
                # metadata
                self.database.description = "Monte Carlo Simulation outputs with SMARTpy v{}.".format(__version__)
                # dimensions
                self.database.createDimension('NbSamples', len(self.p_map))
                self.database.createDimension('NbParameters', len(self.model.parameters.names))
                self.database.createDimension('NbObjFunctions', len(self.obj_fn_names))
                # variables
                params = self.database.createVariable('Parameters', np.float32, ('NbSamples', 'NbParameters'))
                params.units = ', '.join(self.model.parameters.names)
                objfns = self.database.createVariable('ObjFunctions', np.float32, ('NbSamples', 'NbObjFunctions'))
                objfns.units = ', '.join(self.obj_fn_names)
                if self.save_sim:
                    # dimension
                    self.database.createDimension('DateTime', len(self.model.flow))  # Unlimited dimension
                    # coordinate variables
                    times = self.database.createVariable('DateTime', np.float64, ('DateTime',))
                    times.units = "seconds since 1970-01-01 00:00:00.0"
                    # variable
                    simu = self.database.createVariable('Simulations', np.float32, ('NbSamples', 'DateTime'))
                    simu.units = "Discharge in m3/s"
                    # write simulation datetime series as Unix timestamp series
                    datetimes = [np.datetime64(dt) for dt in self.model.flow]
                    timestamps = (datetimes - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')
                    self.database.variables['DateTime'][0:len(datetimes)] = timestamps
            else:
                raise Exception("The use of 'netcdf' as the output file format requires the package 'netCDF4', "
                                "please install it and retry, or choose another file format.")

        else:  # fall back on to default option (CSV file)
            self.database = open(self.db_file, 'w', newline='', encoding='utf8')
            simu_steps = [dt.strftime('%Y-%m-%d %H:%M:%S') for dt in self.model.flow] if self.save_sim else []
            # write header in database file
            self.database.write(','.join(self.obj_fn_names + self.param_names + simu_steps) + '\n')

    def parameters(self):
        return spotpy.parameter.generate(self.params)

    def run(self, compression=None):
        """Run the simulations for the sample of parameter sets.

        :Parameters:

            compression: `bool` or `int`
                Whether the sample output file should be compressed.
                If the output file format is `'csv'`, a `bool` is
                expected (`True` if compression is requested). If the
                output file format is `'netcdf'`, a `bool` or an `int`
                is expected to determine the *complevel* of the Python
                package `netCDF4` which ranges between 1 and 9 (if
                `True`, set to compression level of 6, if `int` is
                provided, the integer value is the *complevel*). If
                not provided, set to default value `None` (i.e. no
                compression is performed).

        """
        # initialise the database (either CSV file or NETCDF file)
        self._init_db()
        # run the Monte Carlo simulation
        sampler = spotpy.algorithms.mc(self, dbformat='custom', parallel=self.parallel)
        sampler.sample(len(self.p_map))
        self.database.close()
        # if compression argument given, the file created will be compressed
        if self.out_format == 'netcdf':
            if compression is True:
                compression = 6
            if not isinstance(compression, bool) and isinstance(compression, (int, float)):
                with Dataset(self.db_file, 'r') as src, Dataset(self.db_file.replace('.nc', '_.nc'), 'w') as dst:
                    dst.description = src.description
                    for name, dimension in src.dimensions.items():
                        dst.createDimension(name, len(dimension))
                    for name, variable in src.variables.items():
                        v = dst.createVariable(name, variable.datatype, variable.dimensions,
                                               zlib=True, complevel=compression)
                        v.units = src.variables[name].units
                        dst.variables[name][:] = src.variables[name][:]
                remove(self.db_file)
                rename(self.db_file.replace('.nc', '_.nc'), self.db_file)
        elif self.out_format == 'csv':
            if compression is True:
                with open(self.db_file, 'rb') as f_in:
                    with gzip.open(self.db_file + '.gz', 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                remove(self.db_file)

    def simulation(self, vector):
        discharge, groundwater_component = self.model.simulate(
            {'T': vector[0], 'C': vector[1], 'H': vector[2], 'D': vector[3], 'S': vector[4], 'Z': vector[5],
             'SK': vector[6], 'FK': vector[7], 'GK': vector[8], 'RK': vector[9]}
        )
        return (
            discharge, [groundwater_component]
        )

    def evaluation(self):
        return (
            self.model.nd_flow, [self.constraints['gw']]
        )

    def objectivefunction(self, simulation, evaluation):
        # select the series subset that have observations (i.e. not NaN)
        flow_eval = evaluation[0][~np.isnan(evaluation[0])]
        flow_simu = np.asarray(simulation[0])[~np.isnan(evaluation[0])]

        # calculate the objective functions
        obj1 = spotpy.objectivefunctions.nashsutcliffe(evaluation=flow_eval, simulation=flow_simu)
        obj2, obj2c, obj2a, obj2b = \
            spotpy.objectivefunctions.kge(evaluation=flow_eval, simulation=flow_simu, return_all=True)
        obj3 = spotpy.objectivefunctions.pbias(evaluation=flow_eval, simulation=flow_simu)
        obj4 = spotpy.objectivefunctions.rmse(evaluation=flow_eval, simulation=flow_simu)

        if self.constraints['gw']:
            return [obj1, obj2, obj2c, obj2a, obj2b, obj3, obj4,
                    groundwater_constraint(evaluation=evaluation[1], simulation=simulation[1])]
        else:
            return [obj1, obj2, obj2c, obj2a, obj2b, obj3, obj4]

    def save(self, obj_fns, parameters, simulations, *args, **kwargs):
        if self.out_format == 'netcdf':
            # it was already checked if netCDF4 was installed and imported, so just proceed
            # convert the parameter values array to a list
            params = parameters.tolist()
            # get the index of the particular parameter set to write data in the right spot
            index = self.p_map[tuple(params)]
            # write data on parameter values and objective functions
            self.database.variables['Parameters'][index, 0:len(self.param_names)] = params
            self.database.variables['ObjFunctions'][index, 0:len(self.obj_fn_names)] = obj_fns
            if self.save_sim:
                # write simulation discharge time series
                my_s, my_e = 0, len(self.model.flow)
                self.database.variables['Simulations'][index, my_s:my_e] = simulations[0]
        else:
            if self.save_sim:  # create a list of objectives functions + parameters + discharge series
                line = map(np.float32, obj_fns + parameters.tolist() + simulations[0].tolist())
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
                with open(file_location, 'r', encoding='utf8') as my_file:
                    my_reader = DictReader(my_file)
                    obj_fns, params = list(), list()
                    for row in my_reader:
                        obj_fns.append([row[obj_fn] for obj_fn in obj_fn_names])
                        params.append([row[param] for param in param_names])
        # return parameter values and objective function values in two separate arrays
        return np.array(params, dtype=np.float32), np.array(obj_fns, dtype=np.float32)
