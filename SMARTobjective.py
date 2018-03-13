import numpy


def calculate_obj_fn(simulation, evaluation, method):

    if not len(simulation) == len(evaluation):
        raise Exception("Simulation and Evaluation series have different length.")
    if method == 'nse':
        return nash_sutcliffe(evaluation, simulation)
    elif method == '':
        return groundwater_constraint(evaluation, simulation)


def nash_sutcliffe(evaluation, simulation):
    # convert list to nd arrays
    nda_flows_mod, nda_flows_obs = numpy.asarray(simulation), numpy.asarray(evaluation)

    # calculate mean of observations
    mean_flow_obs = numpy.mean(nda_flows_obs)

    # calculate Nash-Sutcliffe Efficiency
    nse = 1 - (numpy.sum((nda_flows_obs - nda_flows_mod) ** 2) / numpy.sum((nda_flows_obs - mean_flow_obs) ** 2))

    return nse


def groundwater_constraint(evaluation, simulation):
    if (evaluation[0] - 0.1 <= simulation[0]) and (simulation[0] <= evaluation[0] + 0.1):
        return 1.0
    else:
        return 0.0


def bounded_nash_sutcliffe(evaluation, simulation):
    # convert list to nd arrays
    nda_flows_mod, nda_flows_obs = numpy.asarray(simulation), numpy.asarray(evaluation)

    # calculate mean of observations
    mean_flow_obs = numpy.mean(nda_flows_obs)

    # calculate C2M, bounded formulation of NSE
    f = numpy.sum((nda_flows_obs - nda_flows_mod) ** 2)
    f0 = numpy.sum((nda_flows_obs - mean_flow_obs) ** 2)
    c2m = (1 - (f / f0)) / (1 + (f / f0))

    return c2m
