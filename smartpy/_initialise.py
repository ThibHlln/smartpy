

def initialise(
        # parameters
        theta_t, theta_c, theta_h, theta_d, theta_s,
        theta_z, theta_sk, theta_fk, theta_gk, theta_rk,
        # states
        soil_layers_amounts, overland_reservoir_amount,
        drain_reservoir_amount, inter_reservoir_amount,
        shallow_gw_reservoir_amount, deep_gw_reservoir_amount,
        river_reservoir_amount,
        # constants
        timedelta, drainage_area, rho_water,
):
    # initialise soil layers to be half-full
    soil_layers_amounts[0] = theta_z / 6. / 2.

    # initialise linear reservoirs with following assumptions:
    # - 1200mm/yr of annual average rainfall
    # - 45% of which becomes runoff
    # - 10/30/60% split between surface/soil/groundwater runoff
    aar =  (
        1200  # [mm yr-1] is equivalent to [kg m-2 yr-1]
        / (365.25 * 60 * 60)  # convert [kg m-2 yr-1] to [kg m-2 s-1]
    )

    overland_reservoir_amount[0] = aar * 0.45 * 0.10 * theta_sk
    drain_reservoir_amount[0] = aar * 0.45 * 0.15 * theta_sk
    inter_reservoir_amount[0] = aar * 0.45 * 0.15 * theta_fk
    shallow_gw_reservoir_amount[0] = aar * 0.45 * 0.30 * theta_gk
    deep_gw_reservoir_amount[0] = aar * 0.45 * 0.30 * theta_gk

    river_reservoir_amount[0] = aar * theta_sk
