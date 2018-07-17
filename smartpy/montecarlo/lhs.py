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
    raise Exception('montecarlo.lhs requires the package spotpy to be installed.')

from .montecarlo import MonteCarlo


class LHS(MonteCarlo):
    def __init__(self, catchment, root_f, in_format, save_sim=False):
        MonteCarlo.__init__(self, catchment, root_f, in_format, save_sim=save_sim, func='lhs')

        # attribute params is a list of Objects (e.g. Uniform, Normal; all having Base as a parent class)
        # that are callable (where __call__ of Base picks what numpy.random function to use)
        self.params = [
            spotpy.parameter.Uniform(param,
                                     low=self.model.parameters.ranges[param][0],
                                     high=self.model.parameters.ranges[param][1])
            for param in self.param_names]

    def parameters(self):
        # returns an ndarray containing a tuple for each parameter (random value, step, optguess, minbound, maxbound)
        # parameter.generate() loops through each Object in self.params
        # and uses its astuple() method (inherited from the Base class)
        # that itself uses its __call__() (also inherited from the Base class)
        return spotpy.parameter.generate(self.params)

    def run(self, sample_size, parallel):
        sampler = spotpy.algorithms.lhs(self, parallel=parallel)
        sampler.sample(sample_size)
