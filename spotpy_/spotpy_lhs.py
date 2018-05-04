import spotpy
import argparse
from os import path, getcwd, sep

from scripts.SMARTpy import SMART
from scripts.SMARTfiles import get_dict_simulation_settings
from scripts.SMARTobjective import groundwater_constraint, bounded_nash_sutcliffe


class SpotPySetUp(object):
    def __init__(self, catchment, root_f):
        in_f = sep.join([root_f, 'in', catchment, sep])

        c_area, g_area, start, end, delta_simu, delta_report, warm_up, gw_constraint = \
            get_dict_simulation_settings(''.join([in_f, catchment, '.sttngs']))

        self.model = SMART(catchment, c_area, g_area, start, end, delta_simu, delta_report, warm_up, root_f)

        self.constraints = {'gw': gw_constraint}

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
        return [[simulations[dt] for dt in self.model.flow.iterkeys()], [constraint]]

    def evaluation(self):
        return [[observation for observation in self.model.flow.itervalues()], [self.constraints['gw']]]

    def objectivefunction(self, simulation, evaluation):
        obj1 = spotpy.objectivefunctions.nashsutcliffe(evaluation=evaluation[0], simulation=simulation[0])
        obj2 = spotpy.objectivefunctions.lognashsutcliffe(evaluation=evaluation[0], simulation=simulation[0])
        obj3 = spotpy.objectivefunctions.bias(evaluation=evaluation[0], simulation=simulation[0])
        obj4 = spotpy.objectivefunctions.pbias(evaluation=evaluation[0], simulation=simulation[0])
        obj5 = spotpy.objectivefunctions.kge(evaluation=evaluation[0], simulation=simulation[0])
        obj6 = spotpy.objectivefunctions.rmse(evaluation=evaluation[0], simulation=simulation[0])
        obj7 = bounded_nash_sutcliffe(evaluation=evaluation[0], simulation=simulation[0])
        obj8 = groundwater_constraint(evaluation=evaluation[1], simulation=simulation[1])
        if self.constraints['gw'] == -999.0:
            return [obj1, obj2, obj3, obj4, obj5, obj6, obj7]
        else:
            return [obj1, obj2, obj3, obj4, obj5, obj6, obj7, obj8]


def spotpy_instructions(catchment, sample_size, parallel, root_f):

    spotpy_setup = SpotPySetUp(catchment, root_f)

    sampler = spotpy.algorithms.lhs(spotpy_setup, dbname=spotpy_setup.model.out_f + '{}_LHS_SMART'.format(catchment),
                                    dbformat='csv', parallel=parallel, save_sim=False)
    sampler.sample(sample_size)

    results = sampler.getdata()


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
    spotpy_instructions(args.catchment, args.sample_size, parallelisation, smart_root)
