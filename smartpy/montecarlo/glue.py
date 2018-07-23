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

import numpy as np
from builtins import zip, range

try:
    import spotpy
except ImportError:
    raise Exception('montecarlo.glue requires the package spotpy to be installed.')

from .montecarlo import MonteCarlo


class GLUE(MonteCarlo):
    def __init__(self, catchment, root_f, in_format, out_format,
                 conditioning,
                 parallel='seq', save_sim=False, settings_filename=None, decompression_csv=False):
        MonteCarlo.__init__(self, catchment, root_f, in_format, out_format,
                            parallel=parallel, save_sim=save_sim, func='glue', settings_filename=settings_filename)

        # collect the sampling sets from the Monte Carlo simulation (LHS sampling)
        self.sampling_run_file = \
            ''.join([self.model.out_f, catchment, '.SMART.lhs.nc']) if self.out_format == 'netcdf' else \
            ''.join([self.model.out_f, catchment, '.SMART.lhs'])
        self.sampled_params, self.sampled_obj_fns = self._get_sampled_sets_from_file(self.sampling_run_file,
                                                                                     self.param_names,
                                                                                     self.obj_fn_names,
                                                                                     decompression_csv)
        # generate lists for possible condition(s) to distinguish behavioural and non-behavioural sets
        try:  # get the index for each condition
            self.objective_fn_indices = [self.obj_fn_names.index(fn) for fn in conditioning]
        except ValueError:
            raise Exception("One of the names of objective functions for conditioning in GLUE is not recognised."
                            "Please check for typos and case sensitive issues.")
        self.conditions_types = [conditioning[fn][0] for fn in conditioning]
        self.conditions_values = [conditioning[fn][1] for fn in conditioning]

        # extract behavioural sets from sampling sets
        self.behavioural_params = self._get_behavioural_sets(self.sampled_params,
                                                             self.sampled_obj_fns[:, self.objective_fn_indices],
                                                             self.conditions_values,
                                                             self.conditions_types)

        # create a map of parameter sets to give access to a unique index for each set
        self.p_map = {tuple(self.behavioural_params[r, :].tolist()): r for r in range(self.behavioural_params.shape[0])}

        # give list of behavioural parameters
        self.params = [
            spotpy.parameter.List(self.param_names[0], self.behavioural_params[:, 0]),
            spotpy.parameter.List(self.param_names[1], self.behavioural_params[:, 1]),
            spotpy.parameter.List(self.param_names[2], self.behavioural_params[:, 2]),
            spotpy.parameter.List(self.param_names[3], self.behavioural_params[:, 3]),
            spotpy.parameter.List(self.param_names[4], self.behavioural_params[:, 4]),
            spotpy.parameter.List(self.param_names[5], self.behavioural_params[:, 5]),
            spotpy.parameter.List(self.param_names[6], self.behavioural_params[:, 6]),
            spotpy.parameter.List(self.param_names[7], self.behavioural_params[:, 7]),
            spotpy.parameter.List(self.param_names[8], self.behavioural_params[:, 8]),
            spotpy.parameter.List(self.param_names[9], self.behavioural_params[:, 9])
        ]

    @staticmethod
    def _get_behavioural_sets(params, obj_fns, conditions_val, conditions_typ):
        """
        Based on the Generalized Likelihood Uncertainty Estimation (GLUE) methodology
        to identify behavioural models from Monte Carlo simulations

        See:
        Beven, K., Binley, A. (1992) The future of distributed models: Model calibration and uncertainty prediction.
        Hydrological Processes, 6, 279–298. https://doi.org/10.1002/HYP.3360060305

        Beven, K., Freer, J., 2001. Equifinality, data assimilation, and uncertainty estimation in mechanistic
        modelling of complex environmental systems using the GLUE methodology. Journal of Hydrology, 249, 11–29.
        https://doi.org/10.1016/S0022-1694(01)00421-8

        Beven, K., Binley, A., (2014) GLUE: 20 years on. Hydrological Processes, 28, 5897–5918.
        https://doi.org/10.1002/hyp.10082
        """
        # a few checks to make sure arguments given have compatible dimensions
        if obj_fns.ndim != 2:
            raise Exception('The matrix containing the objective functions is not 2D.')
        if params.ndim != 2:
            raise Exception('The matrix containing the parameters is not 2D.')
        if obj_fns.shape[0] != params.shape[0]:
            raise Exception('The matrices containing objective functions and parameters have different sample sizes.')
        if not ((obj_fns.shape[1] == len(conditions_val)) and (obj_fns.shape[1] == len(conditions_typ))):
            raise Exception('The objective function matrix and the conditions matrices '
                            'do not have compatible dimensions.')

        # generate a mask for the condition(s)
        behavioural = np.ones((obj_fns.shape[0],), dtype=bool)
        for obj_fn, values, kind in zip(obj_fns.T, conditions_val, conditions_typ):
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

            behavioural *= selection

        # apply the mask to the original parameters array and return the resulting array (possibly empty)
        return params[behavioural, :]
