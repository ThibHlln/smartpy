

class TimeFrame(object):
    """
    This class defines the temporal attributes of the simulation period. It contains the start and the end of the
    simulation as well as the lists of DateTime series for the simulation time steps and the reporting time steps (that
    can be identical or nested).

    N.B. 1: The simulation gap needs to be a multiple if the reporting gap and the simulation gap can ony
    be lower than or equal to the report gap
    N.B. 2: The start and the end of the simulation are defined by the user, the class always adds one data step
    prior to the start date in order to set the initial conditions, one or more simulation steps are added in
    consequence depending if the simulation step in a multiple of the data step or not (i.e. equal)
    """
    def __init__(self, datetime_start, datetime_end, simu_timedelta, report_timedelta):
        assert datetime_start <= datetime_end, "TimeFrame: Start > End"
        assert report_timedelta.total_seconds() % simu_timedelta.total_seconds() == 0, \
            "Reporting TimeDelta is not a multiple of Simulation TimeDelta."
        assert (datetime_end - datetime_start).total_seconds() % simu_timedelta.total_seconds() == 0, \
            "Simulation Period is not a multiple of Simulation TimeDelta."
        # DateTime of the start of the time period simulated
        self.start = datetime_start
        # DateTime of the end of the time period simulated
        self.end = datetime_end
        # TimeDelta of the simulation
        self.gap_simu = simu_timedelta
        # TimeDelta of the reporting
        self.gap_report = report_timedelta
        # List of DateTime for the reporting (i.e. list of time steps)
        self.series_report = TimeFrame._create_list_datetime(self, 'report')
        # List of DateTime for the simulation (i.e. list of time steps)
        self.series_simu = TimeFrame._create_list_datetime(self, 'simu')

    def _create_list_datetime(self, option):
        """
        This function returns a list of DateTime by using the start and the end of the simulation and the time gap
        (either the reporting time gap or the simulation time gap, using the option parameter to specify which one).

        N.B. For the initial conditions, the function always adds:
            - [if 'report' option] one data step prior to the reporting start date
            - [if 'simu' option] one (or more if reporting gap > simulation gap) simulation step(s)
            prior to the simulation start date

        :param option: choice to specify if function should work on reporting or on simulation series
        :type option: str()
        :return: a list of DateTime
        :rtype: list()
        """
        extent = self.end - self.start
        options = {'report': self.gap_report, 'simu': self.gap_simu}

        start_index = int((self.gap_report.total_seconds()) / (options[option].total_seconds()))
        end_index = int(extent.total_seconds() // (options[option].total_seconds())) + 1

        my_list_datetime = list()
        for factor in xrange(-start_index, end_index, 1):  # add one or more datetime before start
            my_datetime = self.start + factor * options[option]
            my_list_datetime.append(my_datetime)

        return my_list_datetime

    def get_gap_simu(self):
        return self.gap_simu

    def get_gap_report(self):
        return self.gap_report

    def get_series_simu(self):
        return self.series_simu

    def get_series_report(self):
        return self.series_report
