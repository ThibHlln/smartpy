import numpy as np
from datetime import datetime


def q5(nd_flow_series):
    return np.percentile(nd_flow_series, 95)


def q95(nd_flow_series):
    return np.percentile(nd_flow_series, 5)


def q5_norm(nd_flow_series):
    return np.percentile(nd_flow_series, 95) / np.mean(nd_flow_series)


def q95_norm(nd_flow_series):
    return np.percentile(nd_flow_series, 5) / np.mean(nd_flow_series)


def low_flow_var(nd_flow_series, nd_time_series):

    if not nd_flow_series.shape == nd_time_series.shape:
        raise Exception('Signature LowFlowVar: Flow and Time series are not the same length.')

    min_flow_per_year = list()
    mask_on_years = np.ones(nd_flow_series.shape, dtype=bool)

    start_year = nd_time_series[0].year
    end_year = nd_time_series[-1].year

    year = start_year
    while year <= end_year:
        subset_flow_series = nd_flow_series[
            (nd_time_series >= datetime.strptime('01/10/{} 00:00:00'.format(year), '%d/%m/%Y %H:%M:%S')) &
            (nd_time_series <= datetime.strptime('30/09/{} 00:00:00'.format(year + 1), '%d/%m/%Y %H:%M:%S'))
        ]

        if subset_flow_series.shape[0] >= 310:
            min_flow_per_year.append(subset_flow_series.min())
        else:
            mask_on_years *= ~(
                (nd_time_series >= datetime.strptime('01/10/{} 00:00:00'.format(year), '%d/%m/%Y %H:%M:%S')) &
                (nd_time_series <= datetime.strptime('30/09/{} 00:00:00'.format(year + 1), '%d/%m/%Y %H:%M:%S')))
        year += 1

    if min_flow_per_year:
        return np.mean(np.asarray(min_flow_per_year)) / np.median(nd_flow_series[mask_on_years])
    else:
        return 0.0


def high_flow_var(nd_flow_series, nd_time_series):

    if not nd_flow_series.shape == nd_time_series.shape:
        raise Exception('Signature HighFlowVar: Flow and Time series are not the same length.')

    min_flow_per_year = list()
    mask_on_years = np.ones(nd_flow_series.shape, dtype=bool)

    start_year = nd_time_series[0].year
    end_year = nd_time_series[-1].year

    year = start_year
    while year <= end_year:
        subset_flow_series = nd_flow_series[
            (nd_time_series >= datetime.strptime('01/10/{} 00:00:00'.format(year), '%d/%m/%Y %H:%M:%S')) &
            (nd_time_series <= datetime.strptime('30/09/{} 00:00:00'.format(year + 1), '%d/%m/%Y %H:%M:%S'))
            ]

        if subset_flow_series.shape[0] >= 310:
            min_flow_per_year.append(subset_flow_series.max())
        else:
            mask_on_years *= ~(
                    (nd_time_series >= datetime.strptime('01/10/{} 00:00:00'.format(year), '%d/%m/%Y %H:%M:%S')) &
                    (nd_time_series <= datetime.strptime('30/09/{} 00:00:00'.format(year + 1), '%d/%m/%Y %H:%M:%S')))
        year += 1

    if min_flow_per_year:
        return np.mean(np.asarray(min_flow_per_year)) / np.median(nd_flow_series[mask_on_years])
    else:
        return 0.0


def slope_fdc(nd_flow_series, q_start, q_end):
    return (np.percentile(nd_flow_series, 100 - q_end) - np.percentile(nd_flow_series, 100 - q_start)) / \
           (q_end - q_start) / np.mean(nd_flow_series)
