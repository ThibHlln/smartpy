import numpy


def calculate_obj_fn(simulation, evaluation, method):

    if not len(simulation) == len(evaluation):
        raise Exception("Simulation and Evaluation series have different length.")
    if method == 'nse':
        return calculate_nse(simulation, evaluation)


def calculate_nse(flows_mod, flows_obs):

    # Remove the steps that are missing observations
    nda_flows_mod, nda_flows_obs = numpy.asarray(flows_mod), numpy.asarray(flows_obs)

    # calculate mean of observations
    mean_flow_obs = numpy.mean(nda_flows_obs)

    # calculate Nash-Sutcliffe Efficiency
    nse = 1 - (numpy.sum((nda_flows_obs - nda_flows_mod) ** 2) / numpy.sum((nda_flows_obs - mean_flow_obs) ** 2))

    return nse
