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

from builtins import map
import numpy as np
from os import sep

try:
    import spotpy
except ImportError:
    raise Exception('montecarlo requires the package spotpy to be installed.')

from ..smart import SMART
from ..inout import get_dict_simulation_settings, open_csv_wb
from ..objfunctions import \
    groundwater_constraint, bounded_nash_sutcliffe, sqrt_nash_sutcliffe, spearman_rank_corr, mean_abs_rel_error


class MonteCarlo(object):
    def __init__(self, catchment, root_f, in_format, save_sim=False, func=None):
        in_f = sep.join([root_f, 'in', catchment, sep])

        c_area, g_area, start, end, delta_simu, delta_report, warm_up, gw_constraint = \
            get_dict_simulation_settings(''.join([in_f, catchment, '.sttngs']))

        self.model = SMART(catchment, c_area, g_area, start, end, delta_simu, delta_report, warm_up, in_format, root_f)

        self.save_sim = save_sim

        self.constraints = {'gw': gw_constraint}

        self.param_names = self.model.parameters.names
        self.obj_fn_names = \
            ['NSE', 'lgNSE', 'rtNSE', 'C2M', 'KGE', 'KGEc', 'KGEa', 'KGEb',
             'Bias', 'PBias', 'RMSE', 'Rho', 'MARE'] \
            if self.constraints['gw'] == -999.0 else \
            ['NSE', 'lgNSE', 'rtNSE', 'C2M', 'KGE', 'KGEc', 'KGEa', 'KGEb',
             'Bias', 'PBias', 'RMSE', 'Rho', 'MARE', 'GW']

        # set up a database to custom save results
        self.db_file = self.model.out_f + '{}.SMART.{}'.format(catchment, func)
        self.database = open_csv_wb(self.db_file)
        self.simu_steps = [dt.strftime("%Y-%m-%d %H:%M:%S") for dt in self.model.flow] \
            if self.save_sim else []

        # write header in database file
        self.database.write(','.join(self.obj_fn_names + self.param_names + self.simu_steps) + '\n')

    def simulation(self, vector):
        simulations, constraint = self.model.simulate({
            'T': vector[0], 'C': vector[1], 'H': vector[2], 'D': vector[3], 'S': vector[4], 'Z': vector[5],
            'SK': vector[6], 'FK': vector[7], 'GK': vector[8], 'RK': vector[9]
        })
        # returns only simulations with observations
        return [
            [simulations[dt] for dt in self.model.flow],
            [constraint]
        ]

    def evaluation(self):
        return [
            np.asarray([observation for observation in self.model.flow.values()]),
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
        self.database.write(','.join(map(lambda x: '%.8e' % x, line)) + '\n')
