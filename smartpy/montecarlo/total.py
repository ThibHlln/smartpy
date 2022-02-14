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

try:
    import spotpy
except ImportError:
    raise Exception('montecarlo.best requires the package spotpy to be installed.')

from .montecarlo import MonteCarlo


class Total(MonteCarlo):
    """Total is available to run an existing sample of parameter sets
    on a different period, without any conditioning.

    .. important::

       A sampling must have already been performed prior to using this
       functionality, and the unconditioned sample is then used to run
       the model on a different simulation period than the sampling
       simulation period.

    """
    def __init__(self, catchment, root_f, in_format, out_format,
                 parallel='seq', save_sim=False, settings_filename=None, decompression_csv=False):
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
