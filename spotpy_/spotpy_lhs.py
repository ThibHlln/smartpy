import spotpy
import numpy as np
import argparse
from os import path, getcwd, sep

from scripts.SMARTpy import SMART, valid_file_format
from scripts.SMARTfiles import get_dict_simulation_settings
from scripts.SMARTobjective import \
    groundwater_constraint, bounded_nash_sutcliffe, sqrt_nash_sutcliffe, spearman_rank_corr, mean_abs_rel_error


class SpotPySetUp(object):
    def __init__(self, catchment, root_f, in_fmt, save_sim=False):
        in_f = sep.join([root_f, 'in', catchment, sep])

        c_area, g_area, start, end, delta_simu, delta_report, warm_up, gw_constraint = \
            get_dict_simulation_settings(''.join([in_f, catchment, '.sttngs']))

        self.model = SMART(catchment, c_area, g_area, start, end, delta_simu, delta_report, warm_up, in_fmt, root_f)

        self.save_sim = save_sim

        self.constraints = {'gw': gw_constraint}

        self.param_names = ['T', 'C', 'H', 'D', 'S', 'Z', 'SK', 'FK', 'GK', 'RK']
        self.obj_fn_names = \
            ['NSE', 'lgNSE', 'rtNSE', 'C2M', 'KGE', 'KGEc', 'KGEa', 'KGEb', 'Bias', 'PBias', 'RMSE', 'Rho', 'MARE'] \
            if self.constraints['gw'] == -999.0 else \
            ['NSE', 'lgNSE', 'rtNSE', 'C2M', 'KGE', 'KGEc', 'KGEa', 'KGEb', 'Bias', 'PBias', 'RMSE', 'Rho', 'MARE', 'GW']

        # attribute params is a list of Objects (e.g. Uniform, Normal; all having Base as a parent class)
        # that are callable (where __call__ of Base picks what numpy.random function to use)
        self.params = [
            spotpy.parameter.Uniform('T', low=0.9, high=1.1),
            spotpy.parameter.Uniform('C', low=0.0, high=1.0),
            spotpy.parameter.Uniform('H', low=0.0, high=0.3),
            spotpy.parameter.Uniform('D', low=0.0, high=1.0),
            spotpy.parameter.Uniform('S', low=0.0, high=0.013),
            spotpy.parameter.Uniform('Z', low=15.0, high=150.0),
            spotpy.parameter.Uniform('SK', low=1.0, high=240.0),
            spotpy.parameter.Uniform('FK', low=48.0, high=1440.0),
            spotpy.parameter.Uniform('GK', low=1200.0, high=4800.0),
            spotpy.parameter.Uniform('RK', low=1.0, high=96.0)
        ]

        # set up a database to custom save results
        self.database = file(self.model.out_f + '{}.SMART.lhs'.format(catchment), 'wb')
        self.simu_steps = [dt.strftime("%Y-%m-%d %H:%M:%S") for dt in self.model.flow.iterkeys()] \
            if self.save_sim else []
        # write header in database file
        self.database.write(','.join(self.obj_fn_names + self.param_names + self.simu_steps) + '\n')

    def parameters(self):
        # returns an ndarray containing a tuple for each parameter (random value, step, optguess, minbound, maxbound)
        # parameter.generate() loops through each Object in self.params
        # and uses its astuple() method (inherited from the Base class)
        # that itself uses its __call__() (also inherited from the Base class)
        return spotpy.parameter.generate(self.params)

    def simulation(self, vector):
        simulations, constraint = self.model.simulate({
            'T': vector[0], 'C': vector[1], 'H': vector[2], 'D': vector[3], 'S': vector[4], 'Z': vector[5],
            'SK': vector[6], 'FK': vector[7], 'GK': vector[8], 'RK': vector[9]
        })
        # returns only simulations with observations
        return [
            [simulations[dt] for dt in self.model.flow.iterkeys()],
            [constraint]
        ]

    def evaluation(self):
        return [
            np.asarray([observation for observation in self.model.flow.itervalues()]),
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
        if self.save_sim:
            line = map(np.float32, obj_fns + parameters.tolist() + simulations[0])
        else:
            line = map(np.float32, obj_fns + parameters.tolist())
        self.database.write(','.join(map(str, line)) + '\n')


def spotpy_instructions(catchment, sample_size, parallel, in_format, root_f):

    spotpy_setup = SpotPySetUp(catchment, root_f, in_format, save_sim=False)

    sampler = spotpy.algorithms.lhs(spotpy_setup, parallel=parallel)

    sampler.sample(sample_size)


if __name__ == '__main__':
    # Define the root of the SMARTpy package
    if getcwd() == path.dirname(path.realpath(__file__)):  # execution from the directory where the script is
        smart_root = path.realpath('../..')  # move to parent of parent directory of this current python file
    else:  # execution not from the directory where the script is
        smart_root = getcwd()  # keep the current working directory

    # Collect the arguments to set up SPOTPY
    parser = argparse.ArgumentParser(description="simulate lumped catchment hydrology"
                                                 "for one catchment and one time period"
                                                 "using SPOTPY")
    parser.add_argument('catchment', type=str,
                        help="name of the catchment")
    parser.add_argument('sample_size', type=int,
                        help="size of the sample")
    parser.add_argument('-i', '--in_format', type=valid_file_format, default='csv',
                        help="format of input data files [csv or netcdf]")
    parser.add_argument('-s', '--sequence', dest='parallelisation', action='store_false',
                        help="compute each sample in sequence ")
    parser.add_argument('-p', '--parallel', dest='parallelisation', action='store_true',
                        help="parallel computing of the sample  ")
    parser.set_defaults(parallelisation=False)
    args = parser.parse_args()

    # send the relevant argument for parallelisation option
    if args.parallelisation:
        parallelisation = 'mpi'  # use MPI and parallel computing
    else:
        parallelisation = 'seq'  # use traditional sequential computing

    # Call main function containing SPOTPY instructions
    spotpy_instructions(args.catchment, args.sample_size, parallelisation, args.in_format, smart_root)
