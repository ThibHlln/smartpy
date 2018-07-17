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

from csv import DictReader
import numpy as np
from builtins import zip

try:
    import spotpy
except ImportError:
    raise Exception('montecarlo.best requires the package spotpy to be installed.')

from .montecarlo import MonteCarlo


class Best(MonteCarlo):
    def __init__(self, catchment, target, nb_best, constraining, root_f, in_format, save_sim=False):
        MonteCarlo.__init__(self, catchment, root_f, in_format, save_sim=save_sim, func='{}best'.format(nb_best))

        # extract behavioural sets from sampling sets
        self.sampling_run_file = ''.join([self.model.out_f, catchment, '.SMART.lhs'])
        self.sampled_params, self.sampled_obj_fns = self.get_sampled_sets_from_file(self.sampling_run_file,
                                                                                    self.param_names,
                                                                                    self.obj_fn_names)

        try:
            self.target_fn_index = self.obj_fn_names.index(target)
        except ValueError:
            raise Exception("The objective function {} for conditioning in Best is not recognised."
                            "Please check for typos and case sensitive issues.".format(target))
        try:
            self.constraints_indices = [self.obj_fn_names.index(fn) for fn in constraining]
        except ValueError:
            raise Exception("One of the names of constraints in Best is not recognised."
                            "Please check for typos and case sensitive issues.")
        self.constraints_types = [constraining[fn][0] for fn in constraining]
        self.constraints_values = [constraining[fn][1] for fn in constraining]

        self.best_params = self.get_best_sets(self.sampled_params, self.sampled_obj_fns[:, [13]],
                                              self.constraints_values, self.constraints_types,
                                              self.sampled_obj_fns[:, [0]], nb_best)
        # give list of behavioural parameters
        self.params = [
            spotpy.parameter.List(self.param_names[0], self.best_params[:, 0]),
            spotpy.parameter.List(self.param_names[1], self.best_params[:, 1]),
            spotpy.parameter.List(self.param_names[2], self.best_params[:, 2]),
            spotpy.parameter.List(self.param_names[3], self.best_params[:, 3]),
            spotpy.parameter.List(self.param_names[4], self.best_params[:, 4]),
            spotpy.parameter.List(self.param_names[5], self.best_params[:, 5]),
            spotpy.parameter.List(self.param_names[6], self.best_params[:, 6]),
            spotpy.parameter.List(self.param_names[7], self.best_params[:, 7]),
            spotpy.parameter.List(self.param_names[8], self.best_params[:, 8]),
            spotpy.parameter.List(self.param_names[9], self.best_params[:, 9])
        ]

    def parameters(self):
        # returns an ndarray containing a tuple for each parameter (random value, step, optguess, minbound, maxbound)
        # parameter.generate() loops through each Object in self.params
        # and uses its astuple() method (inherited from the Base class)
        # that itself uses its __call__() (also inherited from the Base class)
        return spotpy.parameter.generate(self.params)

    def run(self, parallel):
        sampler = spotpy.algorithms.mc(self, parallel=parallel)
        sampler.sample(self.best_params.shape[0])

    @staticmethod
    def get_sampled_sets_from_file(file_location, param_names, obj_fn_names):
        with open(file_location) as my_file:
            my_reader = DictReader(my_file)
            obj_fns, params = list(), list()
            for row in my_reader:
                obj_fns.append([row[obj_fn] for obj_fn in obj_fn_names])
                params.append([row[param] for param in param_names])

        return np.array(params, dtype=np.float64), np.array(obj_fns, dtype=np.float64)

    @staticmethod
    def get_best_sets(params, constraints_fns, constraints_val, constraints_typ, sort_fn, nb_best):

        if constraints_fns.ndim != 2:
            raise Exception('The matrix containing the constraint functions is not 2D.')
        if params.ndim != 2:
            raise Exception('The matrix containing the parameters is not 2D.')
        if constraints_fns.shape[0] != params.shape[0]:
            raise Exception('The matrices containing constraint functions and parameters have different sample sizes.')
        if not ((constraints_fns.shape[1] == len(constraints_val)) and
                (constraints_fns.shape[1] == len(constraints_typ))):
            raise Exception('The constraint function matrix and the conditions matrices '
                            'do not have compatible dimensions.')

        if sort_fn.shape[0] != params.shape[0]:
            raise Exception('The matrices containing objective functions and parameters have different sample sizes.')
        if nb_best > params.shape[0]:
            raise Exception('The number of best models requested is higher than the sample size.')

        constrained = np.ones((constraints_fns.shape[0],), dtype=bool)
        for obj_fn, values, kind in zip(constraints_fns.T, constraints_val, constraints_typ):
            if kind == 'equal':
                if len(values) == 1:
                    selection = obj_fn == values[0]
                else:
                    raise Exception("The tuple for \"equal\" condition does not contain one and only one element.")
            elif kind == 'min':
                if len(values) == 1:
                    selection = obj_fn >= values[0]
                else:
                    raise Exception("The tuple for \"min\" condition does not contain one and only one element.")
            elif kind == 'max':
                if len(values) == 1:
                    selection = obj_fn <= values[0]
                else:
                    raise Exception("The tuple for \"max\" condition does not contain one and only one element.")
            elif kind == 'inside':
                if len(values) == 2:
                    if values[1] > values[0]:
                        selection = (obj_fn >= values[0]) & (obj_fn <= values[1])
                    else:
                        raise Exception("The two elements of the tuple for \"inside\" are inconsistent.")
                else:
                    raise Exception("The tuple for \"inside\" condition does not contain two and only two elements.")
            elif kind == 'outside':
                if len(values) == 2:
                    if values[1] > values[0]:
                        selection = (obj_fn <= values[0]) & (obj_fn >= values[1])
                    else:
                        raise Exception("The two elements of the tuple for \"outside\" are inconsistent.")
                else:
                    raise Exception("The tuple for \"outside\" condition does not contain two and only two elements.")
            else:
                raise Exception("The type of threshold \"{}\" is not in the database.".format(kind))

            constrained *= selection

        sort_fn_constrained = sort_fn[constrained, :]
        param_constrained = params[constrained, :]

        if nb_best > param_constrained.shape[0]:
            raise Exception('The number of best models requested is higher than the restrained sample size.')

        return param_constrained[sort_fn_constrained[:, 0].argsort()][-nb_best:]
