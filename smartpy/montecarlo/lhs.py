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

import numpy as np
from scipy.stats import uniform
from builtins import range

try:
    import spotpy
except ImportError:
    raise Exception('montecarlo.lhs requires the package spotpy to be installed.')

from .montecarlo import MonteCarlo


class LHS(MonteCarlo):
    """LHS is the available to perform a sampling in the SMART parameter
    space using a Latin Hypercube sampling `McKay et al. (2000)
    <https:doi.org/10.1080/00401706.2000.10485979>`_.
    """
    def __init__(self, catchment, root_f, in_format, out_format,
                 sample_size,
                 parallel='seq', save_sim=False, settings_filename=None):
        """**Instantiation**

        :Parameters:

            catchment: `str`
                A name to identify the catchment of interest in the
                inputs and outputs directories.

            root_f: `str`
                The file path to the directory containing the model
                inputs and outputs. Note that a specific internal
                structure for this directory must be followed:

                .. code-block:: text

                   root_f
                   ├── in
                   │   └── catchment
                   │       ├── catchment.rain
                   │       ├── catchment.peva
                   │       ├── catchment.flow
                   │       └── catchment.sttngs
                   └── out

            in_format: `str`
                The input file format. It can either be `'csv'` or
                `'netcdf'`. Note that in either case, a specific file
                format must be followed.

            out_format: `str`
                The output file format. It can either be `'csv'` or
                `'netcdf'`.

            sample_size: `int`
                The number of parameter sets to sample in the Latin Hypercube.

            parallel: `str`, optional
                Whether the sampling is to performed in parallel (i.e.
                using MPI calls to run several simulations at the same
                time), or in serial (i.e. running simulations sequentially
                one after another). The options are:

                ===============  =======================================
                `'seq'`          Run the simulations one after another.
                `'mpi'`          Run several simulations at the same
                                 time. The number of simultaneous
                                 simulations is determined with the
                                 number of processes using `mpirun -np`.
                ===============  =======================================

                If not provided, set to default value `'seq'`.

            save_sim: `bool`, optional
                Whether to save the simulated discharge time series. If
                not provided, set to default value `False` (i.e. the
                simulated values are not recorded). Note that the sampled
                parameter values as well as a bundle of objective
                functions are always recorded in the sampling output
                file regardless of this argument.

            settings_filename: `str`, optional
                The name of the settings file to use to configure the
                SMART model. This argument is to be used when the
                settings file does not follow the specific file name
                expected by the model, i.e. *{root_f}/in/{catchment}.sttngs*.
                If not provided, set to the specific file name expected.
                Note that regardless of this argument, the settings file
                must be in the inputs folder for the given catchment
                (in other words, absolute paths are not supported here).

        """
        MonteCarlo.__init__(self, catchment, root_f, in_format, out_format,
                            parallel=parallel, save_sim=save_sim, func='lhs', settings_filename=settings_filename)

        # generate a sample of parameter sets from the using Latin Hypercube Sampling
        self.lhs_params = self._get_params_from_lh(sample_size)

        # create a map of parameter sets to give access to a unique index for each set
        self.p_map = {tuple(self.lhs_params[r, :].tolist()): r for r in range(self.lhs_params.shape[0])}
        
        # give list of parameters generated from Latin Hypercube Sampling
        self.params = [
            spotpy.parameter.List(self.param_names[0], self.lhs_params[:, 0]),
            spotpy.parameter.List(self.param_names[1], self.lhs_params[:, 1]),
            spotpy.parameter.List(self.param_names[2], self.lhs_params[:, 2]),
            spotpy.parameter.List(self.param_names[3], self.lhs_params[:, 3]),
            spotpy.parameter.List(self.param_names[4], self.lhs_params[:, 4]),
            spotpy.parameter.List(self.param_names[5], self.lhs_params[:, 5]),
            spotpy.parameter.List(self.param_names[6], self.lhs_params[:, 6]),
            spotpy.parameter.List(self.param_names[7], self.lhs_params[:, 7]),
            spotpy.parameter.List(self.param_names[8], self.lhs_params[:, 8]),
            spotpy.parameter.List(self.param_names[9], self.lhs_params[:, 9])
        ]
        
    def _get_params_from_lh(self, sample_size):
        """
        Based on the Latin Hypercube Sampling strategy first introduced in:
        McKay, M.D., Beckman, R.J., Conover, W.J. (1979) A Comparison of Three Methods for Selecting Values of Input
        Variables in the Analysis of Output from a Computer Code. Technometrics 21, 239. https://doi.org/10.2307/1268522
        """
        # get the limits for the uniform distributions of each parameter
        bounds = np.asarray(
            [[self.model.parameters.ranges[p][0], self.model.parameters.ranges[p][1]] for p in self.param_names],
            dtype=np.float64
        )

        # get number of parameters
        nb_params = len(self.param_names)

        # create a matrix of random values
        random_matrix = np.random.rand(sample_size, nb_params)

        # randomly permute the segment to be used for each parameter
        sampling_plan = np.zeros((sample_size, nb_params), dtype=np.float64)
        for p in range(nb_params):
            sampling_plan[:, p] = np.random.permutation(sample_size)

        # move away from (randomly selected) segment lower bound by a random value
        sampling_plan += random_matrix
        # standardise values to get values between 0 and 1
        sampling_plan /= sample_size

        # use sampling plan and inverse cumulative distribution function (CDF) of uniform dist. to get parameter values
        parameters = np.zeros((sample_size, nb_params), dtype=np.float64)
        for p in range(nb_params):
            parameters[:, p] = uniform.ppf(sampling_plan[:, p], bounds[p][0], bounds[p][1] - bounds[p][0])

        # return the matrix containing the sample of parameter sets
        return parameters
