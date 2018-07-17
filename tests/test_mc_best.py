import unittest
from csv import reader
from smartpy import inout, montecarlo


class TestMonteCarloBest(unittest.TestCase):
    maxDiff = None

    def setUp(self):
        self.setup = montecarlo.Best(
            catchment='ExampleDaily',
            target='NSE',
            nb_best=1,
            constraining={'GW': ['equal', (1.0,)]},
            root_f='examples/',
            in_format='csv',
            save_sim=False
        )

        self.expected_outcome = [
            '1.22914821e-01', '3.30838971e-02', '1.96008831e-01', '6.54817447e-02', '9.98601988e-02', '8.38074327e-01',
            '3.24638963e-01', '4.27357703e-01', '3.10289168e+00', '-5.72642288e+01', '5.17077208e+00', '8.63192499e-01',
            '5.74714243e-01', '1.00000000e+00', '9.39204156e-01', '7.15413690e-01', '5.50904917e-03', '3.95094156e-01',
            '2.18063840e-04', '1.01668373e+02', '6.14925117e+01', '1.26650793e+03', '4.35021484e+03', '8.41531677e+01'
        ]

    def test_compare_one_best(self):
        # run the SMART model
        self.setup.run(parallel='seq')

        # read in values from output folder
        my_res = list()
        with inout.open_csv_rb(self.setup.db_file) as my_file:
            my_reader = reader(my_file)
            # next(my_reader)  # skip first row that is the headers
            for row in my_reader:
                my_res = row

        # compare
        self.assertEqual(
            my_res,
            self.expected_outcome
        )


if __name__ == '__main__':
    unittest.main()
