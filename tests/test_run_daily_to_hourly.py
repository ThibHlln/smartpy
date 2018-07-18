import unittest
from datetime import datetime, timedelta
import smartpy


class TestRunDaily2Hourly(unittest.TestCase):
    maxDiff = None

    def setUp(self):
        self.sm = smartpy.SMART(
            catchment='ExampleDaily',
            catchment_area_m2=175.46 * 1E6,
            start=datetime.strptime('01/01/2007 09:00:00', '%d/%m/%Y %H:%M:%S'),
            end=datetime.strptime('31/12/2016 09:00:00', '%d/%m/%Y %H:%M:%S'),
            time_delta_simu=timedelta(hours=1),
            time_delta_save=timedelta(days=1),
            warm_up_days=365,
            in_format='csv',
            root="examples/",
            gauged_area_m2=175.97 * 1E6
        )

        self.sm.extra = {'aar': 1200, 'r-o_ratio': 0.45, 'r-o_split': (0.10, 0.15, 0.15, 0.30, 0.30)}

        self.expected_outcome = {
            datetime.strptime('2007-01-01 09:00:00', '%Y-%m-%d %H:%M:%S'): 4.1350823716e+00,
            datetime.strptime('2007-01-02 09:00:00', '%Y-%m-%d %H:%M:%S'): 3.7619768927e+00,
            datetime.strptime('2007-01-03 09:00:00', '%Y-%m-%d %H:%M:%S'): 3.9061266625e+00,
            datetime.strptime('2007-01-04 09:00:00', '%Y-%m-%d %H:%M:%S'): 4.9793335365e+00,
            datetime.strptime('2007-01-05 09:00:00', '%Y-%m-%d %H:%M:%S'): 4.2704806967e+00,
            datetime.strptime('2007-01-06 09:00:00', '%Y-%m-%d %H:%M:%S'): 3.6685291375e+00,
            datetime.strptime('2007-01-07 09:00:00', '%Y-%m-%d %H:%M:%S'): 4.2292124270e+00,
            datetime.strptime('2007-01-08 09:00:00', '%Y-%m-%d %H:%M:%S'): 6.2207036611e+00,
            datetime.strptime('2007-01-09 09:00:00', '%Y-%m-%d %H:%M:%S'): 8.7830566298e+00,
            datetime.strptime('2007-01-10 09:00:00', '%Y-%m-%d %H:%M:%S'): 8.7895069192e+00,
            datetime.strptime('2007-01-11 09:00:00', '%Y-%m-%d %H:%M:%S'): 7.3850475756e+00,
            datetime.strptime('2007-01-12 09:00:00', '%Y-%m-%d %H:%M:%S'): 5.8383811983e+00,
            datetime.strptime('2007-01-13 09:00:00', '%Y-%m-%d %H:%M:%S'): 4.7872186105e+00,
            datetime.strptime('2007-01-14 09:00:00', '%Y-%m-%d %H:%M:%S'): 4.3139791338e+00,
            datetime.strptime('2007-01-15 09:00:00', '%Y-%m-%d %H:%M:%S'): 3.7081587731e+00,
            datetime.strptime('2007-01-16 09:00:00', '%Y-%m-%d %H:%M:%S'): 3.4653561053e+00,
            datetime.strptime('2007-01-17 09:00:00', '%Y-%m-%d %H:%M:%S'): 5.6181275183e+00,
            datetime.strptime('2007-01-18 09:00:00', '%Y-%m-%d %H:%M:%S'): 9.0997880749e+00,
            datetime.strptime('2007-01-19 09:00:00', '%Y-%m-%d %H:%M:%S'): 7.9722444608e+00,
            datetime.strptime('2007-01-20 09:00:00', '%Y-%m-%d %H:%M:%S'): 7.6816521225e+00,
            datetime.strptime('2007-01-21 09:00:00', '%Y-%m-%d %H:%M:%S'): 8.8364637433e+00,
            datetime.strptime('2007-01-22 09:00:00', '%Y-%m-%d %H:%M:%S'): 7.8055472483e+00,
            datetime.strptime('2007-01-23 09:00:00', '%Y-%m-%d %H:%M:%S'): 6.2620670255e+00,
            datetime.strptime('2007-01-24 09:00:00', '%Y-%m-%d %H:%M:%S'): 5.1492557676e+00,
            datetime.strptime('2007-01-25 09:00:00', '%Y-%m-%d %H:%M:%S'): 4.3241389832e+00,
            datetime.strptime('2007-01-26 09:00:00', '%Y-%m-%d %H:%M:%S'): 3.7301083641e+00,
            datetime.strptime('2007-01-27 09:00:00', '%Y-%m-%d %H:%M:%S'): 3.2939432928e+00,
            datetime.strptime('2007-01-28 09:00:00', '%Y-%m-%d %H:%M:%S'): 2.9589612419e+00,
            datetime.strptime('2007-01-29 09:00:00', '%Y-%m-%d %H:%M:%S'): 2.6897933617e+00,
            datetime.strptime('2007-01-30 09:00:00', '%Y-%m-%d %H:%M:%S'): 2.4648666979e+00,
            datetime.strptime('2007-01-31 09:00:00', '%Y-%m-%d %H:%M:%S'): 2.2709545106e+00,
            datetime.strptime('2013-12-29 09:00:00', '%Y-%m-%d %H:%M:%S'): 4.6507377154e+00,
            datetime.strptime('2013-12-30 09:00:00', '%Y-%m-%d %H:%M:%S'): 9.5016415784e+00,
            datetime.strptime('2013-12-31 09:00:00', '%Y-%m-%d %H:%M:%S'): 1.0176161763e+01,
            datetime.strptime('2014-01-01 09:00:00', '%Y-%m-%d %H:%M:%S'): 8.7378649742e+00,
            datetime.strptime('2014-01-02 09:00:00', '%Y-%m-%d %H:%M:%S'): 7.8691683656e+00,
            datetime.strptime('2014-01-03 09:00:00', '%Y-%m-%d %H:%M:%S'): 6.4532835261e+00,
            datetime.strptime('2014-01-04 09:00:00', '%Y-%m-%d %H:%M:%S'): 5.1871304107e+00,
            datetime.strptime('2014-01-05 09:00:00', '%Y-%m-%d %H:%M:%S'): 4.8977507445e+00,
            datetime.strptime('2014-01-06 09:00:00', '%Y-%m-%d %H:%M:%S'): 5.6701859532e+00,
            datetime.strptime('2014-01-07 09:00:00', '%Y-%m-%d %H:%M:%S'): 5.1945615951e+00,
            datetime.strptime('2014-01-08 09:00:00', '%Y-%m-%d %H:%M:%S'): 4.3213346888e+00,
            datetime.strptime('2014-01-09 09:00:00', '%Y-%m-%d %H:%M:%S'): 3.6363404567e+00,
            datetime.strptime('2014-01-10 09:00:00', '%Y-%m-%d %H:%M:%S'): 3.1948731408e+00,
            datetime.strptime('2014-01-11 09:00:00', '%Y-%m-%d %H:%M:%S'): 2.9398471219e+00,
            datetime.strptime('2014-01-12 09:00:00', '%Y-%m-%d %H:%M:%S'): 2.6553518068e+00,
            datetime.strptime('2014-01-13 09:00:00', '%Y-%m-%d %H:%M:%S'): 2.4427429699e+00,
            datetime.strptime('2014-01-14 09:00:00', '%Y-%m-%d %H:%M:%S'): 3.0575433650e+00,
            datetime.strptime('2014-01-15 09:00:00', '%Y-%m-%d %H:%M:%S'): 4.9930975160e+00,
            datetime.strptime('2014-01-16 09:00:00', '%Y-%m-%d %H:%M:%S'): 4.6855643808e+00,
            datetime.strptime('2014-01-17 09:00:00', '%Y-%m-%d %H:%M:%S'): 3.8042933811e+00,
            datetime.strptime('2014-01-18 09:00:00', '%Y-%m-%d %H:%M:%S'): 3.5317948321e+00,
            datetime.strptime('2016-11-23 09:00:00', '%Y-%m-%d %H:%M:%S'): 1.9328004620e-01,
            datetime.strptime('2016-11-24 09:00:00', '%Y-%m-%d %H:%M:%S'): 1.5444063791e-01,
            datetime.strptime('2016-11-25 09:00:00', '%Y-%m-%d %H:%M:%S'): 1.2978545372e-01,
            datetime.strptime('2016-11-26 09:00:00', '%Y-%m-%d %H:%M:%S'): 1.1368431010e-01,
            datetime.strptime('2016-11-27 09:00:00', '%Y-%m-%d %H:%M:%S'): 1.0275405485e-01,
            datetime.strptime('2016-11-28 09:00:00', '%Y-%m-%d %H:%M:%S'): 9.4972608597e-02,
            datetime.strptime('2016-11-29 09:00:00', '%Y-%m-%d %H:%M:%S'): 8.9130991154e-02,
            datetime.strptime('2016-11-30 09:00:00', '%Y-%m-%d %H:%M:%S'): 8.4505019744e-02,
            datetime.strptime('2016-12-01 09:00:00', '%Y-%m-%d %H:%M:%S'): 8.2564646362e-02,
            datetime.strptime('2016-12-02 09:00:00', '%Y-%m-%d %H:%M:%S'): 8.4631401458e-02,
            datetime.strptime('2016-12-03 09:00:00', '%Y-%m-%d %H:%M:%S'): 8.2604017773e-02,
            datetime.strptime('2016-12-04 09:00:00', '%Y-%m-%d %H:%M:%S'): 7.9642608244e-02,
            datetime.strptime('2016-12-05 09:00:00', '%Y-%m-%d %H:%M:%S'): 7.6788151599e-02,
            datetime.strptime('2016-12-06 09:00:00', '%Y-%m-%d %H:%M:%S'): 8.2671359651e-02,
            datetime.strptime('2016-12-07 09:00:00', '%Y-%m-%d %H:%M:%S'): 9.9691651793e-02,
            datetime.strptime('2016-12-08 09:00:00', '%Y-%m-%d %H:%M:%S'): 9.3218430758e-02,
            datetime.strptime('2016-12-09 09:00:00', '%Y-%m-%d %H:%M:%S'): 8.4345634836e-02,
            datetime.strptime('2016-12-10 09:00:00', '%Y-%m-%d %H:%M:%S'): 7.7756437465e-02,
            datetime.strptime('2016-12-11 09:00:00', '%Y-%m-%d %H:%M:%S'): 8.0946569189e-02,
            datetime.strptime('2016-12-12 09:00:00', '%Y-%m-%d %H:%M:%S'): 1.7667816907e-01,
            datetime.strptime('2016-12-13 09:00:00', '%Y-%m-%d %H:%M:%S'): 4.1854101173e-01,
            datetime.strptime('2016-12-14 09:00:00', '%Y-%m-%d %H:%M:%S'): 6.9752519595e-01,
            datetime.strptime('2016-12-15 09:00:00', '%Y-%m-%d %H:%M:%S'): 1.4274367506e+00,
            datetime.strptime('2016-12-16 09:00:00', '%Y-%m-%d %H:%M:%S'): 2.1898963789e+00,
            datetime.strptime('2016-12-17 09:00:00', '%Y-%m-%d %H:%M:%S'): 1.6407987667e+00,
            datetime.strptime('2016-12-18 09:00:00', '%Y-%m-%d %H:%M:%S'): 1.0488554343e+00,
            datetime.strptime('2016-12-19 09:00:00', '%Y-%m-%d %H:%M:%S'): 6.6712867060e-01,
            datetime.strptime('2016-12-20 09:00:00', '%Y-%m-%d %H:%M:%S'): 8.1267005754e-01,
            datetime.strptime('2016-12-21 09:00:00', '%Y-%m-%d %H:%M:%S'): 1.5448290230e+00,
            datetime.strptime('2016-12-22 09:00:00', '%Y-%m-%d %H:%M:%S'): 1.2790547619e+00,
            datetime.strptime('2016-12-23 09:00:00', '%Y-%m-%d %H:%M:%S'): 1.2916680138e+00,
            datetime.strptime('2016-12-24 09:00:00', '%Y-%m-%d %H:%M:%S'): 2.5498859282e+00,
            datetime.strptime('2016-12-25 09:00:00', '%Y-%m-%d %H:%M:%S'): 2.2885019253e+00,
            datetime.strptime('2016-12-26 09:00:00', '%Y-%m-%d %H:%M:%S'): 2.0765606493e+00,
            datetime.strptime('2016-12-27 09:00:00', '%Y-%m-%d %H:%M:%S'): 1.5479103702e+00,
            datetime.strptime('2016-12-28 09:00:00', '%Y-%m-%d %H:%M:%S'): 1.1243331842e+00,
            datetime.strptime('2016-12-29 09:00:00', '%Y-%m-%d %H:%M:%S'): 8.5075383672e-01,
            datetime.strptime('2016-12-30 09:00:00', '%Y-%m-%d %H:%M:%S'): 6.7547748371e-01,
            datetime.strptime('2016-12-31 09:00:00', '%Y-%m-%d %H:%M:%S'): 8.3091723923e-01
        }

    def test_compare_discharge_series(self):
        # get the parameter values
        self.sm.parameters.set_parameters_with_file(
            ''.join([self.sm.in_f, self.sm.catchment, '.parameters']))

        # run the SMART model
        self.sm.simulate(self.sm.parameters.values)

        # round all values to 12 decimals
        my_res = dict()
        my_ref = dict()
        for dt in self.expected_outcome:
            my_res[dt] = '%.6e' % self.sm.discharge[dt]
            my_ref[dt] = '%.6e' % self.expected_outcome[dt]

        # compare
        self.assertDictEqual(
            my_res,
            my_ref
        )


if __name__ == '__main__':
    unittest.main()
