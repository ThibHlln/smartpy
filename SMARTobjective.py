import numpy as np
from scipy.stats import spearmanr


def nash_sutcliffe(evaluation, simulation):
    # convert list to nd arrays
    nda_flows_mod, nda_flows_obs = np.asarray(simulation), np.asarray(evaluation)

    # calculate Nash-Sutcliffe Efficiency
    nse = 1 - (np.sum((nda_flows_obs - nda_flows_mod) ** 2) / np.sum((nda_flows_obs - np.mean(nda_flows_obs)) ** 2))

    return nse


def log_nash_sutcliffe(evaluation, simulation):
    # convert list to nd arrays
    nda_flows_mod, nda_flows_obs = np.asarray(simulation), np.asarray(evaluation)

    # calculate Square Root Nash-Sutcliffe Efficiency
    sqrt_nse = 1 - (np.sum((np.log(nda_flows_obs) - np.log(nda_flows_mod)) ** 2) /
                    np.sum((np.log(nda_flows_obs) - np.mean(np.log(nda_flows_obs))) ** 2))

    return sqrt_nse


def sqrt_nash_sutcliffe(evaluation, simulation):
    # convert list to nd arrays
    nda_flows_mod, nda_flows_obs = np.asarray(simulation), np.asarray(evaluation)

    # calculate Square Root Nash-Sutcliffe Efficiency
    sqrt_nse = 1 - (np.sum((np.sqrt(nda_flows_obs) - np.sqrt(nda_flows_mod)) ** 2) /
                    np.sum((np.sqrt(nda_flows_obs) - np.mean(np.sqrt(nda_flows_obs))) ** 2))

    return sqrt_nse


def bounded_nash_sutcliffe(evaluation, simulation):
    # convert list to nd arrays
    nda_flows_mod, nda_flows_obs = np.asarray(simulation), np.asarray(evaluation)

    # calculate mean of observations
    mean_flow_obs = np.mean(nda_flows_obs)

    # calculate C2M, bounded formulation of NSE
    f = np.sum((nda_flows_obs - nda_flows_mod) ** 2)
    f0 = np.sum((nda_flows_obs - mean_flow_obs) ** 2)
    c2m = (1 - (f / f0)) / (1 + (f / f0))

    return c2m


def groundwater_constraint(evaluation, simulation):
    if (evaluation[0] - 0.1 <= simulation[0]) and (simulation[0] <= evaluation[0] + 0.1):
        return 1.0
    else:
        return 0.0


def kling_gupta(evaluation, simulation):
    # convert list to nd arrays
    nda_flows_mod, nda_flows_obs = np.asarray(simulation), np.asarray(evaluation)

    # calculate Kling-Gupta Efficiency
    cc = np.corrcoef(nda_flows_obs, nda_flows_mod)[0, 1]  # correlation coefficient (error in dynamics)
    alpha = np.std(nda_flows_mod) / np.std(nda_flows_obs)  # alpha (error in variability)
    beta = np.sum(nda_flows_mod) / np.sum(nda_flows_obs)  # central tendency beta (error in volume)
    kge = 1 - np.sqrt((cc - 1) ** 2 + (alpha - 1) ** 2 + (beta - 1) ** 2)

    return kge


def spearman_rank_corr(evaluation, simulation):
    # convert list to nd arrays
    nda_flows_mod, nda_flows_obs = np.asarray(simulation), np.asarray(evaluation)

    # return correlation only (rho)
    return spearmanr(nda_flows_mod, nda_flows_obs)[0]


def mean_abs_rel_error(evaluation, simulation):
    # convert list to nd arrays
    nda_flows_mod, nda_flows_obs = np.asarray(simulation), np.asarray(evaluation)

    # calculate mean absolute relative error (MARE)
    mare = np.sum(np.abs(nda_flows_obs - nda_flows_mod)) / np.sum(nda_flows_obs)

    return mare


def percent_bias(evaluation, simulation):
    # convert list to nd arrays
    nda_flows_mod, nda_flows_obs = np.asarray(simulation), np.asarray(evaluation)

    # calculate percent bias
    pbias = 100 * np.sum(nda_flows_obs - nda_flows_mod) / np.sum(nda_flows_obs)

    return pbias
