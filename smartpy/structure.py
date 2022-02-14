# This file is part of SMARTpy - An open-source rainfall-runoff model in Python
# Copyright (C) 2018-2022  Thibault Hallouin (1), Eva Mockler (1,2), Michael Bruen (1)
#
# (1) Dooge Centre for Water Resources Research, University College Dublin, Ireland
# (2) Environmental Protection Agency, Ireland
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

import numpy as np

try:
    import smartcpp
    smart_in_cpp = True
except ImportError:
    smartcpp = None
    smart_in_cpp = False


def run(area_m2, delta,
        nd_rain, nd_peva,
        nd_parameters, extra,
        timeseries, timeseries_report,
        report,
        **kwargs):
    """
    This function defines the structure of the database using model outputs and states names
    and defines the initial conditions of the different reservoirs in the hydrological model
    (either using a warm-up run if keyword argument 'warm_up', or with an educated guess if
    'extra' contains the necessary information)

    :param area_m2: catchment area in square meters
    :param delta: duration between two simulations time steps as a TimeDelta object
    :param nd_rain: numpy array containing the rainfall series
    :param nd_peva: numpy array containing the potential evapotranspiration series
    :param nd_parameters: numpy array containing the parameter values
    :param extra: dictionary containing optional information for educated guess of initial conditions
    :param timeseries: list of simulation time steps as DateTime objects
    :param timeseries_report: list of reporting time steps as DateTime objects
    :param report: reporting type (either 'summary' or 'raw')
    :param kwargs: additional keyword arguments
    :return: numpy array for discharge series and scalar for groundwater component to runoff
    """

    # determine whether SMARTcpp (C++ extension of SMART for Python) can be used or not
    if smart_in_cpp:
        try:  # all steps is available for SMARTcpp version >= 0.2.0
            all_steps = smartcpp.allsteps
        except AttributeError:  # i.e. SMARTcpp version < 0.2.0 (will try again with one step)
            all_steps = run_all_steps
    else:
        all_steps = run_all_steps

    # convert report to report type as an integer
    if report == 'summary':
        report_type = 1
    elif report == 'raw':
        report_type = 2
    else:
        raise Exception('Reporting type \'{}\' unknown.'.format(report))

    # determine temporal information
    simu_length = len(timeseries) - 1
    delta_sec = delta.total_seconds()
    report_gap = (len(timeseries) - 1) // (len(timeseries_report) - 1)

    # set up the data structures
    model_outputs = ['Q_aeva', 'Q_ove', 'Q_dra', 'Q_int', 'Q_sgw', 'Q_dgw', 'Q_out']
    model_states = ['V_ove', 'V_dra', 'V_int', 'V_sgw', 'V_dgw',
                    'V_ly1', 'V_ly2', 'V_ly3', 'V_ly4', 'V_ly5', 'V_ly6',
                    'V_river']
    model_variables = model_outputs + model_states

    nd_initial = np.zeros((len(model_variables),), dtype=np.float64)

    # set initial conditions
    if kwargs['warm_up'] != 0:  # either using warm-up run
        warm_up_length = int(kwargs['warm_up'] * 86400 / delta.total_seconds())
        
        if warm_up_length > simu_length:
            raise Exception("The warm-up duration (i.e. {} days) cannot exceed the length of the simulation period "
                            "because the beginning of the simulation period is used as made-up warm-up data for the "
                            "sake of model states initialisation. Please specify another warm-up duration to comply "
                            "with this requirement, or consider using actual warm-up data at the beginning of the "
                            "simulation period and set the warm-up period to 0.".format(kwargs['warm_up']))
            
        nd_initial_wu = np.zeros((len(model_variables),), dtype=np.float64)

        # start with non-empty linear reservoirs (1200mm/yr SAAR & 45% becomes runoff, split 60/30/10% GW/Soil/Surface)
        if extra:
            nd_initial_wu[model_variables.index('V_ove')] = \
                (extra['aar'] * extra['r-o_ratio']) * extra['r-o_split'][0] / 1000 * area_m2 / 8766 * nd_parameters[6]
            nd_initial_wu[model_variables.index('V_dra')] = \
                (extra['aar'] * extra['r-o_ratio']) * extra['r-o_split'][1] / 1000 * area_m2 / 8766 * nd_parameters[6]
            nd_initial_wu[model_variables.index('V_int')] = \
                (extra['aar'] * extra['r-o_ratio']) * extra['r-o_split'][2] / 1000 * area_m2 / 8766 * nd_parameters[7]
            nd_initial_wu[model_variables.index('V_sgw')] = \
                (extra['aar'] * extra['r-o_ratio']) * extra['r-o_split'][3] / 1000 * area_m2 / 8766 * nd_parameters[8]
            nd_initial_wu[model_variables.index('V_dgw')] = \
                (extra['aar'] * extra['r-o_ratio']) * extra['r-o_split'][4] / 1000 * area_m2 / 8766 * nd_parameters[8]
            nd_initial_wu[model_variables.index('V_river')] = \
                (extra['aar'] * extra['r-o_ratio']) / 1000 * area_m2 / 8766 * nd_parameters[9]

        # start with soil layers half full (six soil layers so half gives /12)
        nd_initial_wu[model_variables.index('V_ly1'):(model_variables.index('V_ly6') + 1)] = \
            (nd_parameters[5] / 12) / 1000 * area_m2

        nd_initial[:] = all_steps(area_m2, delta_sec, warm_up_length,
                                  nd_rain, nd_peva,
                                  nd_parameters, nd_initial_wu,
                                  report_type, report_gap)[2]

    else:  # or starting without warm-up run
        # start with non-empty linear reservoirs (1200mm/yr SAAR & 45% becomes runoff, split 10/30/60% Surface/Soil/GW)
        if extra:
            nd_initial[model_variables.index('V_ove')] = \
                (extra['aar'] * extra['r-o_ratio']) * extra['r-o_split'][0] / 1000 * area_m2 / 8766 * nd_parameters[6]
            nd_initial[model_variables.index('V_dra')] = \
                (extra['aar'] * extra['r-o_ratio']) * extra['r-o_split'][1] / 1000 * area_m2 / 8766 * nd_parameters[6]
            nd_initial[model_variables.index('V_int')] = \
                (extra['aar'] * extra['r-o_ratio']) * extra['r-o_split'][2] / 1000 * area_m2 / 8766 * nd_parameters[7]
            nd_initial[model_variables.index('V_sgw')] = \
                (extra['aar'] * extra['r-o_ratio']) * extra['r-o_split'][3] / 1000 * area_m2 / 8766 * nd_parameters[8]
            nd_initial[model_variables.index('V_dgw')] = \
                (extra['aar'] * extra['r-o_ratio']) * extra['r-o_split'][4] / 1000 * area_m2 / 8766 * nd_parameters[8]
            nd_initial[model_variables.index('V_river')] = \
                (extra['aar'] * extra['r-o_ratio']) / 1000 * area_m2 / 8766 * nd_parameters[9]
        # start with soil layers half full (six soil layers so half gives /12)
        nd_initial[model_variables.index('V_ly1'):(model_variables.index('V_ly6') + 1)] = \
            (nd_parameters[5] / 12) / 1000 * area_m2

    # run the actual simulation
    return all_steps(area_m2, delta_sec, simu_length,
                     nd_rain, nd_peva,
                     nd_parameters, nd_initial,
                     report_type, report_gap)[0:2]


def run_all_steps(area_m2, delta_sec, length_simu,
                  nd_rain, nd_peva,
                  nd_parameters, nd_initial,
                  report_type, report_gap):
    """
    This function calls the hydrological model (run_one_step) for all the time steps in the period in sequence.
    It returns the simulated discharge at the outlet for the simulation period, as well as the proportion of runoff
    contributed by the two groundwater pathways.

    :param area_m2: catchment area in square meters
    :param delta_sec: duration between two simulations time steps in seconds
    :param length_simu: number of simulation time steps (excluding initial step)
    :param nd_rain: numpy array containing the rainfall series
    :param nd_peva: numpy array containing the potential evapotranspiration series
    :param nd_parameters: numpy array containing the parameter values
    :param nd_initial: numpy array containing outputs and states (axis 1) for each simulation time step (axis 0)
    :param report_type: reporting type (either 1 for 'summary' or 2 for 'raw')
    :param report_gap: number of simulation time steps contained in one reporting time step
    :return: numpy array for discharge series and scalar for groundwater component to runoff
    """

    # determine whether SMARTcpp (C++ extension of SMART for Python) can be used or not
    if smart_in_cpp:  # i.e. SMARTcpp version < 0.2.0 (because SMARTcpp did not have all steps)
        one_step = smartcpp.onestep
    else:
        one_step = run_one_step

    # create a local data structure to store simulations
    nd_storage = np.zeros((length_simu + 1, nd_initial.shape[0],), dtype=np.float64)
    # assign initial conditions into storage array
    nd_storage[0, :] = nd_initial[:]
    # run simulations for each time step
    for i in range(1, length_simu + 1):
        nd_storage[i, :] = one_step(area_m2, delta_sec,
                                    nd_rain[i - 1], nd_peva[i - 1],
                                    nd_parameters[0], nd_parameters[1], nd_parameters[2], nd_parameters[3],
                                    nd_parameters[4], nd_parameters[5], nd_parameters[6], nd_parameters[7],
                                    nd_parameters[8], nd_parameters[9],
                                    *nd_storage[i - 1, 7:])
    # extract the reporting information from the simulations
    if report_type == 1:  # summary of the simulation steps corresponding to one reporting step
        discharge = np.mean(np.reshape(nd_storage[1:, 6], (-1, report_gap)), axis=-1)
        groundwater_component = np.sum(nd_storage[1:, [4, 5]]) / np.sum(nd_storage[1:, [1, 2, 3, 4, 5]])
    else:  # raw extraction of the value corresponding to the reporting datetime
        discharge = nd_storage[1:, 6][::-report_gap][::-1]
        groundwater_component = np.sum(nd_storage[1:, [4, 5]][::-report_gap, :][::-1, :]) / \
            np.sum(nd_storage[1:, [1, 2, 3, 4, 5]][::-report_gap, :][::-1, :])

    return discharge, groundwater_component, nd_storage[-1]


def run_one_step(area_m2, time_delta_sec,
                 c_in_rain, c_in_peva,
                 c_p_t, c_p_c, c_p_h, c_p_d, c_p_s, c_p_z, c_p_sk, c_p_fk, c_p_gk, r_p_rk,
                 c_s_v_h2o_ove, c_s_v_h2o_dra, c_s_v_h2o_int, c_s_v_h2o_sgw, c_s_v_h2o_dgw,
                 c_s_v_h2o_ly1, c_s_v_h2o_ly2, c_s_v_h2o_ly3, c_s_v_h2o_ly4, c_s_v_h2o_ly5, c_s_v_h2o_ly6,
                 r_s_v_riv
                 ):
    """
    This function links catchment model outputs and river routing model input.
    Altogether they form the hydrological model running one step.
    It returns all the outputs and states of the two models as a tuple.

    :param area_m2: catchment area in square meters
    :param time_delta_sec: simulation time step in seconds
    :param c_in_rain: cumulative rainfall in millimeters per time step for the simulation time step
    :param c_in_peva: cumulative potential evapotranspiration in millimeters per time step for the simulation time step
    :param c_p_t: SMART parameter T
    :param c_p_c: SMART parameter C
    :param c_p_h: SMART parameter H
    :param c_p_d: SMART parameter D
    :param c_p_s: SMART parameter S
    :param c_p_z: SMART parameter Z
    :param c_p_sk: SMART parameter SK
    :param c_p_fk: SMART parameter FK
    :param c_p_gk: SMART parameter GK
    :param r_p_rk: SMART parameter RK
    :param c_s_v_h2o_ove: SMART catchment state for overland flow reservoir
    :param c_s_v_h2o_dra: SMART catchment state for drain flow reservoir
    :param c_s_v_h2o_int: SMART catchment state for inter flow reservoir
    :param c_s_v_h2o_sgw: SMART catchment state for shallow groundwater flow reservoir
    :param c_s_v_h2o_dgw: SMART catchment state for deep groundwater flow reservoir
    :param c_s_v_h2o_ly1: SMART catchment state for first soil layer
    :param c_s_v_h2o_ly2: SMART catchment state for second soil layer
    :param c_s_v_h2o_ly3: SMART catchment state for third soil layer
    :param c_s_v_h2o_ly4: SMART catchment state for fourth soil layer
    :param c_s_v_h2o_ly5: SMART catchment state for fifth soil layer
    :param c_s_v_h2o_ly6: SMART catchment state for sixth soil layer
    :param r_s_v_riv: SMART river state for channel reservoir
    :return: SMART outputs and states
    """
    out_aeva, out_q_h2o_ove, out_q_h2o_dra, out_q_h2o_int, out_q_h2o_sgw, out_q_h2o_dgw, \
        s_v_h2o_ove, s_v_h2o_dra, s_v_h2o_int, s_v_h2o_sgw, s_v_h2o_dgw, \
        s_v_h2o_ly1, s_v_h2o_ly2, s_v_h2o_ly3, s_v_h2o_ly4, s_v_h2o_ly5, s_v_h2o_ly6 = \
        run_one_step_catchment(
            area_m2, time_delta_sec,
            c_in_rain, c_in_peva,
            c_p_t, c_p_c, c_p_h, c_p_d, c_p_s, c_p_z, c_p_sk, c_p_fk, c_p_gk,
            c_s_v_h2o_ove, c_s_v_h2o_dra, c_s_v_h2o_int, c_s_v_h2o_sgw, c_s_v_h2o_dgw,
            c_s_v_h2o_ly1, c_s_v_h2o_ly2, c_s_v_h2o_ly3, c_s_v_h2o_ly4, c_s_v_h2o_ly5, c_s_v_h2o_ly6
        )

    out_q_riv, s_v_riv = \
        run_one_step_river(
            time_delta_sec,
            out_q_h2o_ove + out_q_h2o_dra + out_q_h2o_int + out_q_h2o_sgw + out_q_h2o_dgw,
            r_p_rk,
            r_s_v_riv
        )

    return (
        out_aeva, out_q_h2o_ove, out_q_h2o_dra, out_q_h2o_int, out_q_h2o_sgw, out_q_h2o_dgw, out_q_riv,  # outputs
        s_v_h2o_ove, s_v_h2o_dra, s_v_h2o_int, s_v_h2o_sgw, s_v_h2o_dgw,   # states
        s_v_h2o_ly1, s_v_h2o_ly2, s_v_h2o_ly3, s_v_h2o_ly4, s_v_h2o_ly5, s_v_h2o_ly6,
        s_v_riv
    )


def run_one_step_catchment(area_m2, time_gap_sec,
                           c_in_rain, c_in_peva,
                           c_p_t, c_p_c, c_p_h, c_p_d, c_p_s, c_p_z, c_p_sk, c_p_fk, c_p_gk,
                           c_s_v_ove, c_s_v_dra, c_s_v_int, c_s_v_sgw, c_s_v_dgw,
                           c_s_v_ly1, c_s_v_ly2, c_s_v_ly3, c_s_v_ly4, c_s_v_ly5, c_s_v_ly6
                           ):
    """
    This function was written by Thibault Hallouin but is largely inspired by the work of Eva Mockler, namely for
    the work published in: Mockler, E., O’Loughlin, F., and Bruen, M.: Understanding hydrological flow paths in
    conceptual catchment models using uncertainty and sensitivity analysis, Computers & Geosciences, 90, 66–77,
    doi:10.1016/j.cageo.2015.08.015, 2016.

    Catchment model * c_ *
    _ Hydrology
    ___ Inputs * in_ *
    _____ c_in_rain         precipitation as rain [mm/time step]
    _____ c_in_peva         potential evapotranspiration [mm/time step]
    ___ Parameters * p_ *
    _____ c_p_t             T: rainfall aerial correction coefficient
    _____ c_p_c             C: evaporation decay parameter
    _____ c_p_h             H: quick runoff coefficient
    _____ c_p_d             D: drain flow parameter - fraction of saturation excess diverted to drain flow
    _____ c_p_s             S: soil outflow coefficient
    _____ c_p_z             Z: effective soil depth [mm]
    _____ c_p_sk            SK: surface routing parameter [hours]
    _____ c_p_fk            FK: inter flow routing parameter [hours]
    _____ c_p_gk            GK: groundwater routing parameter [hours]
    ___ States * s_ *
    _____ c_s_v_ove         volume of water in overland store [m3]
    _____ c_s_v_dra         volume of water in drain store [m3]
    _____ c_s_v_int         volume of water in inter store [m3]
    _____ c_s_v_sgw         volume of water in shallow groundwater store [m3]
    _____ c_s_v_dgw         volume of water in deep groundwater store [m3]
    _____ c_s_v_ly1         volume of water in first soil layer store [m3]
    _____ c_s_v_ly2         volume of water in second soil layer store [m3]
    _____ c_s_v_ly3         volume of water in third soil layer store [m3]
    _____ c_s_v_ly4         volume of water in fourth soil layer store [m3]
    _____ c_s_v_ly5         volume of water in fifth soil layer store [m3]
    _____ c_s_v_ly6         volume of water in sixth soil layer store [m3]
    ___ Outputs * out_ *
    _____ c_out_aeva        actual evapotranspiration [m3/s]
    _____ c_out_q_ove       overland flow [m3/s]
    _____ c_out_q_dra       drain flow [m3/s]
    _____ c_out_q_int       inter flow [m3/s]
    _____ c_out_q_sgw       shallow groundwater flow [m3/s]
    _____ c_out_q_dgw       deep groundwater flow [m3/s]
    """

    # # 1. Hydrology
    # # 1.0. Define internal constants
    nb_soil_layers = 6.0  # number of layers in soil column [-]

    # # 1.1. Convert non-SI units
    c_p_sk *= 3600.0  # convert hours in seconds
    c_p_fk *= 3600.0  # convert hours in seconds
    c_p_gk *= 3600.0  # convert hours in seconds

    # # 1.2. Hydrological calculations

    # /!\ all calculations in mm equivalent until further notice

    # calculate capacity Z and level LVL of each layer (assumed equal) from effective soil depth
    list_z_lyr = [
        0.0,  # artificial null value added to keep script clear later
        c_p_z / nb_soil_layers,  # Soil Layer 1
        c_p_z / nb_soil_layers,  # Soil Layer 2
        c_p_z / nb_soil_layers,  # Soil Layer 3
        c_p_z / nb_soil_layers,  # Soil Layer 4
        c_p_z / nb_soil_layers,  # Soil Layer 5
        c_p_z / nb_soil_layers  # Soil Layer 6
    ]

    list_lvl_lyr = [  # factor 1000 to convert m in mm
        0.0,  # artificial null value added to keep script clear later
        c_s_v_ly1 / area_m2 * 1e3,  # Soil Layer 1
        c_s_v_ly2 / area_m2 * 1e3,  # Soil Layer 2
        c_s_v_ly3 / area_m2 * 1e3,  # Soil Layer 3
        c_s_v_ly4 / area_m2 * 1e3,  # Soil Layer 4
        c_s_v_ly5 / area_m2 * 1e3,  # Soil Layer 5
        c_s_v_ly6 / area_m2 * 1e3  # Soil Layer 6
    ]

    # calculate cumulative level of water in all soil layers at beginning of time step (i.e. soil moisture)
    lvl_total_start = sum(list_lvl_lyr)

    # apply parameter T to rainfall data (aerial rainfall correction)
    rain = c_in_rain * c_p_t
    # calculate excess rainfall
    excess_rain = rain - c_in_peva
    # initialise actual evapotranspiration variable
    aeva = 0.0

    if excess_rain >= 0.0:  # excess rainfall available for runoff and infiltration
        # actual evapotranspiration = potential evapotranspiration
        aeva += c_in_peva
        # calculate surface runoff using quick runoff parameter H and relative soil moisture content
        h_prime = c_p_h * (lvl_total_start / c_p_z)
        overland_flow = h_prime * excess_rain  # excess rainfall contribution to quick surface runoff store
        excess_rain -= overland_flow  # remainder that infiltrates
        # calculate percolation through soil layers (from top layer [1] to bottom layer [6])
        for i in [1, 2, 3, 4, 5, 6]:
            space_in_lyr = list_z_lyr[i] - list_lvl_lyr[i]
            if excess_rain <= space_in_lyr:
                list_lvl_lyr[i] += excess_rain
                excess_rain = 0.0
            else:
                list_lvl_lyr[i] = list_z_lyr[i]
                excess_rain -= space_in_lyr
        # calculate saturation excess from remaining excess rainfall after filling layers (if not 0)
        drain_flow = c_p_d * excess_rain  # sat. excess contribution (if not 0) to quick interflow runoff store
        inter_flow = (1.0 - c_p_d) * excess_rain  # sat. excess contribution (if not 0) to slow interflow runoff store
        # calculate leak from soil layers (i.e. piston flow becoming active during rainfall events)
        s_prime = c_p_s * (lvl_total_start / c_p_z)
        # leak to interflow
        for i in [1, 2, 3, 4, 5, 6]:  # soil moisture outflow reducing exponentially downwards
            leak_interflow = list_lvl_lyr[i] * (s_prime ** i)
            if leak_interflow < list_lvl_lyr[i]:
                inter_flow += leak_interflow  # soil moisture outflow contribution to slow interflow runoff store
                list_lvl_lyr[i] -= leak_interflow
        # leak to shallow groundwater flow
        shallow_flow = 0.0
        for i in [1, 2, 3, 4, 5, 6]:  # soil moisture outflow reducing linearly downwards
            leak_shallow_flow = list_lvl_lyr[i] * (s_prime / i)
            if leak_shallow_flow < list_lvl_lyr[i]:
                shallow_flow += leak_shallow_flow  # soil moisture outflow contribution to slow shallow GW runoff store
                list_lvl_lyr[i] -= leak_shallow_flow
        # leak to deep groundwater flow
        deep_flow = 0.0
        for i in [6, 5, 4, 3, 2, 1]:  # soil moisture outflow reducing exponentially upwards
            leak_deep_flow = list_lvl_lyr[i] * (s_prime ** (7 - i))
            if leak_deep_flow < list_lvl_lyr[i]:
                deep_flow += leak_deep_flow  # soil moisture outflow contribution to slow deep GW runoff store
                list_lvl_lyr[i] -= leak_deep_flow
    else:  # no excess rainfall (i.e. potential evapotranspiration not satisfied by available rainfall)
        overland_flow = 0.0  # no soil moisture contribution to quick overland flow runoff store
        drain_flow = 0.0  # no soil moisture contribution to quick drain flow runoff store
        inter_flow = 0.0  # no soil moisture contribution to quick + leak interflow runoff store
        shallow_flow = 0.0  # no soil moisture contribution to shallow groundwater flow runoff store
        deep_flow = 0.0  # no soil moisture contribution to deep groundwater flow runoff store

        deficit_rain = excess_rain * (-1.0)  # excess is negative => excess is actually a deficit
        aeva += rain
        for i in [1, 2, 3, 4, 5, 6]:  # attempt to satisfy PE from soil layers (from top layer [1] to bottom layer [6]
            if list_lvl_lyr[i] >= deficit_rain:  # i.e. all moisture required available in this soil layer
                list_lvl_lyr[i] -= deficit_rain  # soil layer is reduced by the moisture required
                aeva += deficit_rain  # this moisture contributes to the actual evapotranspiration
                deficit_rain = 0.0  # the full moisture still required has been met
            else:  # i.e. not all moisture required available in this soil layer
                aeva += list_lvl_lyr[i]  # takes what is available in this layer for evapotranspiration
                # effectively reduce the evapotranspiration demand for the next layer using parameter C
                # i.e. the more you move down through the soil layers, the less AET can meet PET (exponentially)
                deficit_rain = c_p_c * (deficit_rain - list_lvl_lyr[i])
                list_lvl_lyr[i] = 0.0  # soil layer is now empty

    # /!\ all calculations in S.I. units now (i.e. mm converted into cubic metres)

    # calculate actual evapotranspiration as a flux
    c_out_aeva = aeva / 1e3 * area_m2 / time_gap_sec  # [m3/s]

    # route overland flow (quick surface runoff)
    c_out_q_h2o_ove = c_s_v_ove / c_p_sk  # [m3/s]
    c_s_v_ove += (overland_flow / 1e3 * area_m2) - (c_out_q_h2o_ove * time_gap_sec)  # [m3] - [m3]
    if c_s_v_ove < 0.0:
        c_s_v_ove = 0.0
    # route drain flow (quick interflow runoff)
    c_out_q_h2o_dra = c_s_v_dra / c_p_sk  # [m3/s]
    c_s_v_dra += (drain_flow / 1e3 * area_m2) - (c_out_q_h2o_dra * time_gap_sec)  # [m3] - [m3]
    if c_s_v_dra < 0.0:
        c_s_v_dra = 0.0
    # route interflow (slow interflow runoff)
    c_out_q_h2o_int = c_s_v_int / c_p_fk  # [m3/s]
    c_s_v_int += (inter_flow / 1e3 * area_m2) - (c_out_q_h2o_int * time_gap_sec)  # [m3] - [m3]
    if c_s_v_int < 0.0:
        c_s_v_int = 0.0
    # route shallow groundwater flow (slow shallow GW runoff)
    c_out_q_h2o_sgw = c_s_v_sgw / c_p_gk  # [m3/s]
    c_s_v_sgw += (shallow_flow / 1e3 * area_m2) - (c_out_q_h2o_sgw * time_gap_sec)  # [m3] - [m3]
    if c_s_v_sgw < 0.0:
        c_s_v_sgw = 0.0
    # route deep groundwater flow (slow deep GW runoff)
    c_out_q_h2o_dgw = c_s_v_dgw / c_p_gk  # [m3/s]
    c_s_v_dgw += (deep_flow / 1e3 * area_m2) - (c_out_q_h2o_dgw * time_gap_sec)  # [m3] - [m3]
    if c_s_v_dgw < 0.0:
        c_s_v_dgw = 0.0

    # # 1.3. Return outputs and states
    return (
        c_out_aeva, c_out_q_h2o_ove, c_out_q_h2o_dra, c_out_q_h2o_int, c_out_q_h2o_sgw, c_out_q_h2o_dgw,
        c_s_v_ove, c_s_v_dra, c_s_v_int, c_s_v_sgw, c_s_v_dgw,
        list_lvl_lyr[1] / 1e3 * area_m2, list_lvl_lyr[2] / 1e3 * area_m2, list_lvl_lyr[3] / 1e3 * area_m2,
        list_lvl_lyr[4] / 1e3 * area_m2, list_lvl_lyr[5] / 1e3 * area_m2, list_lvl_lyr[6] / 1e3 * area_m2
    )


def run_one_step_river(time_gap_sec,
                       r_in_q_riv, r_p_rk, r_s_v_riv):
    """
    This function was written by Thibault Hallouin but is largely inspired by the work of Eva Mockler, namely for
    the work published in: Mockler, E., O’Loughlin, F., and Bruen, M.: Understanding hydrological flow paths in
    conceptual catchment models using uncertainty and sensitivity analysis, Computers & Geosciences, 90, 66–77,
    doi:10.1016/j.cageo.2015.08.015, 2016.

    River model * r_ *
    _ Hydrology
    ___ Inputs * in_ *
    _____ r_in_q_riv      flow at inlet [m3/s]
    ___ Parameters * p_ *
    _____ r_p_rk          linear factor k for water where Storage = k.Flow [hours]
    ___ States * s_ *
    _____ r_s_v_riv       volume of water in store [m3]
    ___ Outputs * out_ *
    _____ r_out_q_riv     flow at outlet [m3/s]
    """
    # # 1. Hydrology
    # # 1.0. Define internal constants
    r_p_rk *= 3600.0  # convert hours in seconds

    # # 1.1. Hydrological calculations

    # calculate outflow, at current time step
    r_out_q_riv = r_s_v_riv / r_p_rk
    # calculate storage in temporary variable, for next time step
    r_s_v_h2o_old = r_s_v_riv
    r_s_v_h2o_temp = r_s_v_h2o_old + (r_in_q_riv - r_out_q_riv) * time_gap_sec
    # check if storage has gone negative
    if r_s_v_h2o_temp < 0.0:  # temporary cannot be used
        # constrain outflow: allow maximum outflow at 95% of what was in store
        r_out_q_riv = 0.95 * (r_in_q_riv + r_s_v_h2o_old / time_gap_sec)
        # calculate final storage with constrained outflow
        r_s_v_riv += (r_in_q_riv - r_out_q_riv) * time_gap_sec
    else:
        r_s_v_riv = r_s_v_h2o_temp  # temporary storage becomes final storage

    # # 1.2. Return output and state
    return (
        r_out_q_riv, r_s_v_riv
    )
