from os import path, makedirs, sep

from .timeframe import TimeFrame
from .parameters import Parameters
from .inout import \
    get_dict_rain_series_simu, get_dict_peva_series_simu, get_dict_discharge_series, write_flow_file_from_dict
from .structure import run


class SMART(object):
    def __init__(self, catchment, c_area_m2, g_area_m2, start, end,
                 time_delta_simu, time_delta_report, warm_up_days, in_fmt, root):
        # general information
        self.catchment = catchment
        self.area = c_area_m2
        # directory information
        self.root_f = root
        self.in_f = sep.join([self.root_f, 'in', self.catchment, sep])
        self.out_f = sep.join([self.root_f, 'out', self.catchment, sep])
        if not path.exists(self.out_f):
            makedirs(self.out_f)
        # temporal information
        self.start = start
        self.end = end
        self.delta_simu = time_delta_simu
        self.delta_report = time_delta_report
        self.timeframe = TimeFrame(self.start, self.end, self.delta_simu, self.delta_report)
        self.timeseries = self.timeframe.get_series_simu()
        self.timeseries_report = self.timeframe.get_series_report()
        self.warm_up = warm_up_days
        # physical information
        extra_ext = '.nc' if in_fmt == 'netcdf' else ''
        self.rain = get_dict_rain_series_simu(''.join([self.in_f, self.catchment, '.rain' + extra_ext]), in_fmt,
                                              self.timeseries[1], self.timeseries[-1], self.delta_simu)
        self.peva = get_dict_peva_series_simu(''.join([self.in_f, self.catchment, '.peva' + extra_ext]), in_fmt,
                                              self.timeseries[1], self.timeseries[-1], self.delta_simu)
        self.flow = get_dict_discharge_series(''.join([self.in_f, self.catchment, '.flow']),
                                              self.timeframe.get_series_report()[1],
                                              self.timeframe.get_series_report()[-1],
                                              c_area_m2, g_area_m2)
        # parameters
        self.parameters = Parameters()
        # model outputs
        self.outputs = None
        self.discharge = None

    def simulate(self, param):
        db = dict()
        self.outputs = run(self.area, self.delta_simu, self.rain, self.peva,
                           param, db, self.timeseries, self.timeseries_report,
                           warm_up=self.warm_up)
        self.discharge = self.outputs[0]
        return self.outputs

    def write_output_files(self):
        if self.outputs:
            write_flow_file_from_dict(self.timeframe, self.discharge,
                                      ''.join([self.out_f, self.catchment, '.mod.flow']),
                                      method='raw')

            write_flow_file_from_dict(self.timeframe, self.flow,
                                      ''.join([self.out_f, self.catchment, '.obs.flow']),
                                      method='raw')
        else:
            raise Exception("The output files cannot be written because the outputs attribute is unassigned."
                            "Please make sure to call the simulate method before writing the output files.")
