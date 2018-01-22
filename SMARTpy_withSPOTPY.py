import spotpy
from os import path
from datetime import datetime, timedelta

from SMARTpy import SMART


class SpotPySetUp(object):
    def __init__(self):
        catchment = 'Avonmore_l_IE_EA_10A050300'
        area = 230.41164787028104
        start = datetime(2011, 1, 1, 9, 0, 0)
        end = datetime(2016, 12, 31, 9, 0, 0)
        delta_simu = timedelta(minutes=60)
        delta_report = timedelta(minutes=1440)
        warm_up = 365
        root_f = path.realpath('..')

        self.model = SMART(catchment, area, start, end, delta_simu, delta_report, warm_up, root_f)

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
        simulations = self.model.simulate({
            'T': vector[0], 'C': vector[1], 'H': vector[2], 'D': vector[3], 'S': vector[4], 'Z': vector[5],
            'SK': vector[6], 'FK': vector[7], 'GK': vector[8], 'RK': vector[9]
        })
        # returns only simulations with observations
        return [simulations[dt] for dt in self.model.flow.iterkeys()]

    def evaluation(self):
        return [observation for observation in self.model.flow.itervalues()]

    def objectivefunction(self, simulation, evaluation):
        return - spotpy.objectivefunctions.rmse(evaluation, simulation)


if __name__ == '__main__':

    spotpy_setup = SpotPySetUp()

    sampler = spotpy.algorithms.lhs(spotpy_setup, dbname='C:/PycharmProjects/Python/SMARTpy/out/LHS_SMART',
                                    dbformat='csv', parallel='seq')

    sampler.sample(5)
    results = sampler.getdata()
