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
from builtins import zip, range

try:
    import spotpy
except ImportError:
    raise Exception('montecarlo.best requires the package spotpy to be installed.')

from .montecarlo import MonteCarlo


class Best(MonteCarlo):
    """Best is the available to condition a sample of parameter sets
    by selecting the best performing parameter set(s) based on a given
    objective function used as target.

    .. important::

       A sampling must have already been performed prior to using this
       functionality, and the conditioned sample is then used to run
       the model on a different simulation period than the sampling
       simulation period.

    """
    def __init__(self, catchment, root_f, in_format, out_format,
                 target, nb_best, constraining=None,
                 parallel='seq', save_sim=False, settings_filename=None,
                 decompression_csv=False):
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

            target: `str`
                The name of the objective function to use as target for
                the determination of which parameter sets are the best
                at reproducing the observed river discharge. The options
                are:

                =========  ====================================================
                'NSE'      `Nash Sutcliffe Efficiency
                           <https://doi.org/10.1016/0022-1694(70)90255-6>`_.
                'KGE'      `Kling-Gupta Efficiency
                           <https://doi.org/10.1016/j.jhydrol.2009.08.003>`_.
                'KGEc'     *r* component of `Kling-Gupta Efficiency
                           <https://doi.org/10.1016/j.jhydrol.2009.08.003>`_.
                'KGEa'     *alpha* component of `Kling-Gupta Efficiency
                           <https://doi.org/10.1016/j.jhydrol.2009.08.003>`_.
                'KGEb'     *beta* component of `Kling-Gupta Efficiency
                           <https://doi.org/10.1016/j.jhydrol.2009.08.003>`_.
                'PBias'    Percent bias.
                'RMSE'     Root mean square error.
                =========  ====================================================

            nb_best: `int`
                The number of parameter sets to select as best. This can
                be any value between 1 and the *sample_size*.

            constraining: `dict`
                An option to filter the sample of parameter sets before
                selecting the best performing one(s). This is the same
                behaviour as *conditioning* in `montecarlo.GLUE`.

                *Parameter example:* ::

                    constraining={
                        'GW': ['equal', (1,)]  # GW = 1
                    }

                *Parameter example:* ::

                    constraining={
                        'PBias': ['inside', (-10, 10)],  # -10% <= PBias <= 10%
                        'KGE': ['min', (0.6,)]  # KGE >= 0.6
                    }

            parallel: `str`, optional
                Whether the sampling is to performed in parallel (i.e.
                using MPI calls to run several simulations at the same
                time), or in serial (i.e. running simulations sequentially
                one after another). The options are:

                ===============  =======================================
                Parallel         Description
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

            decompression_csv: `bool`, optional
                Whether the CSV files containing the sample of parameter
                sets is compressed or not. If it is, this must be set to
                `True` to decompress the CSV file as a pre-processing
                step. If not provided, set to default vaalue `False` (i.e.
                no decompression).

        """
        MonteCarlo.__init__(self, catchment, root_f, in_format, out_format,
                            parallel=parallel, save_sim=save_sim, func='{}best'.format(nb_best),
                            settings_filename=settings_filename)

        # collect the sampling sets from the Monte Carlo simulation (LHS sampling)
        self.sampling_run_file = \
            ''.join([self.model.out_f, catchment, '.SMART.lhs.nc']) if self.out_format == 'netcdf' else \
            ''.join([self.model.out_f, catchment, '.SMART.lhs'])
        self.sampled_params, self.sampled_obj_fns = self._get_sampled_sets_from_file(self.sampling_run_file,
                                                                                     self.param_names,
                                                                                     self.obj_fn_names,
                                                                                     decompression_csv)
        # get the index of the targeted objective functions
        try:
            self.target_fn_index = [self.obj_fn_names.index(target)]
        except ValueError:
            raise Exception("The objective function {} for conditioning in Best is not recognised."
                            "Please check for typos and case sensitive issues.".format(target))

        # generate lists for possible constraint(s) before selecting the best parameter sets
        if constraining:
            try:  # get the indices of the possible constraint(s)
                self.constraints_indices = [self.obj_fn_names.index(fn) for fn in constraining]
            except ValueError:
                raise Exception("One of the names of constraints in Best is not recognised."
                                "Please check for typos and case sensitive issues.")
            self.constraints_types = [constraining[fn][0] for fn in constraining]
            self.constraints_values = [constraining[fn][1] for fn in constraining]
        else:
            self.constraints_indices, self.constraints_types, self.constraints_values = [], [], []

        # extract the best parameter sets given the target and the possible constraint(s)
        self.best_params = self._get_best_sets(self.sampled_params, self.sampled_obj_fns[:, self.constraints_indices],
                                               self.constraints_values, self.constraints_types,
                                               self.sampled_obj_fns[:, self.target_fn_index], nb_best)

        # create a map of parameter sets to give access to a unique index for each set
        self.p_map = {tuple(self.best_params[r, :].tolist()): r for r in range(self.best_params.shape[0])}

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

    @staticmethod
    def _get_best_sets(params, constraints_fns, constraints_val, constraints_typ, sort_fn, nb_best):
        # a few checks to make sure arguments given have compatible dimensions
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

        # generate a mask for the possible constraint(s)
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

        # apply the constraint(s) mask to
        sort_fn_constrained = sort_fn[constrained, :]
        param_constrained = params[constrained, :]

        # check that there is enough parameter sets remaining after constraint to meet the sample size
        if nb_best > param_constrained.shape[0]:
            raise Exception('The number of best models requested is higher than the restrained sample size.')

        # return the best parameter sets on the given targeted objective functions
        return param_constrained[sort_fn_constrained[:, 0].argsort()][-nb_best:]
