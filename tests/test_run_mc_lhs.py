import unittest
from smartpy import montecarlo


class TestRunMonteCarloLHS(unittest.TestCase):
    maxDiff = None

    def setUp(self):
        self.catchment = 'Catchment'
        self.root_f = 'data/'
        self.sample_size = 5
        self.parallel = 'seq'
        self.save_sim = True

        self.extra = {
            'aar': 1200,
            'r-o_ratio': 0.45,
            'r-o_split': (0.10, 0.15, 0.15, 0.30, 0.30)
        }

    def test_run_lhs_terminates(self):
        for in_fmt in ('csv', 'netcdf'):
            for out_fmt in ('csv', 'netcdf'):
                with self.subTest(in_format=in_fmt, out_format=out_fmt):
                    # initialise the Monte Carlo LHS
                    setup = montecarlo.LHS(
                        catchment=self.catchment,
                        root_f=self.root_f,
                        in_format=in_fmt,
                        out_format=out_fmt,
                        sample_size=self.sample_size,
                        parallel=self.parallel,
                        save_sim=self.save_sim
                    )

                    # add some extra information to the model
                    setup.model.extra = self.extra

                    # run the simulation (with compression test)
                    setup.run(compression=True)


if __name__ == '__main__':
    test_loader = unittest.TestLoader()
    test_suite = unittest.TestSuite()

    test_suite.addTests(
        test_loader.loadTestsFromTestCase(TestRunMonteCarloLHS)
    )

    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(test_suite)

    if not result.wasSuccessful():
        exit(1)
