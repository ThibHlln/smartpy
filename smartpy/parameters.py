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

from csv import DictReader


class Parameters(object):
    def __init__(self):
        #: Return the SMART model parameter names as a `list`.
        self.names = ['T', 'C', 'H', 'D', 'S', 'Z', 'SK', 'FK', 'GK', 'RK']
        #: Return the typical SMART model parameter ranges as a `dict`.
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
            'RK': (1.0, 96.0)
        }
        #: Return the set of SMART model parameter values as a `dict`.
        self.values = dict()

    def set_parameters_with_file(self, file_location):
        """Assign the SMART model parameters values using a CSV file.

        :Parameters:

            file_location: `str`
                The absolute file path for the parameters file.

                *Parameter example:* ::

                    file_location='examples/in/ExampleDaily/ExampleDaily.parameters'

        """
        my_dict_par = dict()
        try:
            with open(file_location, 'r', encoding='utf8') as my_file:
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

    def set_parameters_with_dict(self, dictionary):
        """Assign the SMART model parameters values using a dictionary.

        :Parameters:

            dictionary: `dict`
                The dictionary containing the parameter values.

                *Parameter example:* ::

                    dictionary={
                        'T': 1.0,
                        'C': 1.0,
                        'H': 0.20845,
                        'D': 0.24606,
                        'S': 0.0001230,
                        'Z': 105.26,
                        'SK': 46.82,
                        'FK': 315.55,
                        'GK': 1066.73,
                        'RK': 10.64
                    }

        """
        for param in self.names:
            try:
                self.values[param] = dictionary[param]
            except KeyError:
                raise Exception("The parameter {} is not available in the dictionary provided.")
