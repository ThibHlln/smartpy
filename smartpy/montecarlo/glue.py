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
    raise Exception('montecarlo.glue requires the package spotpy to be installed.')

from .montecarlo import MonteCarlo


class GLUE(MonteCarlo):
    """GLUE is the available to condition a sample of parameter sets
    using a Generalized Likelihood Uncertainty Estimation (GLUE) approach
    (`Beven and Binley (1992) <https://doi.org/10.1002/HYP.3360060305>`_,
    `Beven and Freer (2001) <https://doi.org/10.1016/S0022-1694(01)00421-8>`_,
    `Beven and Binley (2014) <https://doi.org/10.1002/hyp.10082>`_). That
    is to say elicit which parameter sets in a sample are "behavioural"
    or not, based on some likelihood measure.

    .. important::

       A sampling must have already been performed prior to using this
       functionality, and the conditioned sample is then used to run
       the model on a different simulation period than the sampling
       simulation period.

    """
    def __init__(self, catchment, root_f, in_format, out_format,
                 conditioning,
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

            conditioning: `dict`
                The set of conditions to use to consider a parameter set
                as "behavioural" gathered in a dictionary. The objective
                functions to choose from as likelihood measures are:

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
                'GW'       Groundwater contribution to runoff. This is a
                           specific objective function for the SMART model.
                           This objective function acts as a filter to
                           eliminate parameter sets who produce ratios of
                           groundwater runoff, shallow and deep, over total
                           runoff that are not within :math:`\pm` 10% of
                           the observed/expected value. The objective
                           function evaluates as one when this is the
                           case, and as zero when this is not (i.e.
                           eliminated).
                =========  ====================================================

                *Parameter example:* ::

                    conditioning={
                        'GW': ['equal', (1,)],  # GW = 1
                        'NSE': ['min', (0.8,)]  # NSE >= 0.8
                        'PBias': ['max', (10,)],  # PBias <= 10%
                    }

                *Parameter example:* ::

                    conditioning={
                        'PBias': ['inside', (-10, 10)],  # -10% <= PBias <= 10%
                        'KGE': ['min', (0.6,)]  # KGE >= 0.6
                    }

                .. warning::

                   Depending on the conditions chosen, this may yield no
                   "behavioural" parameter set.

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
