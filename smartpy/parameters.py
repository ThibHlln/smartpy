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

from csv import DictReader

from .inout import open_csv_rb


class Parameters(object):
    def __init__(self):
        self.names = ['T', 'C', 'H', 'D', 'S', 'Z', 'SK', 'FK', 'GK', 'RK']
        self.ranges = {
            'T': (0.9, 1.1),
            'C': (0.0, 1.0),
            'H': (0.0, 0.3),
            'D': (0.0, 1.0),
            'S': (0.0, 0.013),
            'Z': (15.0, 150.0),
            'SK': (1.0, 240.0),
            'FK': (48.0, 1440.0),
            'GK': (1200.0, 4800.0),
            'RK': (1.0, 96.0),
        }
        self.values = dict()

    def get_parameters_from_file(self, file_location):
        my_dict_par = dict()
        try:
            with open_csv_rb(file_location) as my_file:
                my_reader = DictReader(my_file)
                for row in my_reader:
                    if row['PAR_NAME'] in self.names:
                        my_dict_par[row['PAR_NAME']] = float(row['PAR_VALUE'])
        except KeyError:
            raise Exception("There is 'PAR_NAME' or 'PAR_VALUE' column in {}.".format(file_location))
        except ValueError:
            raise Exception("There is at least one incorrect parameter value in {}.".format(file_location))
        except IOError:
            raise Exception("There is no parameters file at {}.".format(file_location))

        for param in self.names:
            try:
                self.values[param] = my_dict_par[param]
            except KeyError:
                raise Exception("The parameter {} is not available in the "
                                "parameters file at {}.".format(param, file_location))

