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

try:
    import spotpy
except ImportError:
    raise Exception('montecarlo.best requires the package spotpy to be installed.')

from .montecarlo import MonteCarlo


class Total(MonteCarlo):
    def __init__(self, catchment, root_f, in_format, out_format,
                 parallel='seq', save_sim=False, settings_filename=None, decompression_csv=False):
        MonteCarlo.__init__(self, catchment, root_f, in_format, out_format,
                            parallel=parallel, save_sim=save_sim, func='total', settings_filename=settings_filename)

        # collect the sampling sets from the Monte Carlo simulation (LHS sampling)
        self.sampling_run_file = \
            ''.join([self.model.out_f, catchment, '.SMART.lhs.nc']) if self.out_format == 'netcdf' else \
            ''.join([self.model.out_f, catchment, '.SMART.lhs'])
        self.sampled_params, self.sampled_obj_fns = self._get_sampled_sets_from_file(self.sampling_run_file,
                                                                                     self.param_names,
                                                                                     self.obj_fn_names,
                                                                                     decompression_csv)
        # create a map of parameter sets to give access to a unique index for each set
        self.p_map = {tuple(self.sampled_params[r, :].tolist()): r for r in range(self.sampled_params.shape[0])}

        # give list of behavioural parameters
        self.params = [
            spotpy.parameter.List(self.param_names[0], self.sampled_params[:, 0]),
            spotpy.parameter.List(self.param_names[1], self.sampled_params[:, 1]),
            spotpy.parameter.List(self.param_names[2], self.sampled_params[:, 2]),
            spotpy.parameter.List(self.param_names[3], self.sampled_params[:, 3]),
            spotpy.parameter.List(self.param_names[4], self.sampled_params[:, 4]),
            spotpy.parameter.List(self.param_names[5], self.sampled_params[:, 5]),
            spotpy.parameter.List(self.param_names[6], self.sampled_params[:, 6]),
            spotpy.parameter.List(self.param_names[7], self.sampled_params[:, 7]),
            spotpy.parameter.List(self.param_names[8], self.sampled_params[:, 8]),
            spotpy.parameter.List(self.param_names[9], self.sampled_params[:, 9])
        ]
