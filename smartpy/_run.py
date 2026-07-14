import numpy as np


def run(
        # inputs
        rainfall_flux, potential_evapotranspiration_flux,
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
        # outputs
        actual_evapotranspiration_flux, river_discharge_flux,
        # internals
        shallow_gw_flux, deep_gw_flux,
        overland_runoff_flux, drain_runoff_flux,
        inter_runoff_flux, shallow_gw_runoff_flux,
        deep_gw_runoff_flux,
        # number of timesteps
        nt
):
    for i in range(nt):
        update(
            # inputs
            rainfall_flux[i:i+1], potential_evapotranspiration_flux[i:i+1],
            # parameters
            theta_t, theta_c, theta_h, theta_d, theta_s,
            theta_z, theta_sk, theta_fk, theta_gk, theta_rk,
            # states
            soil_layers_amounts[i:i+2], overland_reservoir_amount[i:i+2],
            drain_reservoir_amount[i:i+2], inter_reservoir_amount[i:i+2],
            shallow_gw_reservoir_amount[i:i+2], deep_gw_reservoir_amount[i:i+2],
            river_reservoir_amount[i:i+2],
            # constants
            timedelta, drainage_area, rho_water,
            # outputs
            actual_evapotranspiration_flux[i:i+1], river_discharge_flux[i:i+1],
            # internals
            shallow_gw_flux, deep_gw_flux,
            overland_runoff_flux, drain_runoff_flux,
            inter_runoff_flux, shallow_gw_runoff_flux,
            deep_gw_runoff_flux
        )


def update(
        # inputs
        rainfall_flux, potential_evapotranspiration_flux,
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
        # outputs
        actual_evapotranspiration_flux, river_discharge_flux,
        # internals
        shallow_gw_flux, deep_gw_flux,
        overland_runoff_flux, drain_runoff_flux,
        inter_runoff_flux, shallow_gw_runoff_flux,
        deep_gw_runoff_flux
):
    _update_production(
        rainfall_flux, potential_evapotranspiration_flux,
        theta_t, theta_c, theta_h, theta_d, theta_s,
        theta_z, theta_sk, theta_fk, theta_gk,
        soil_layers_amounts, overland_reservoir_amount,
        drain_reservoir_amount, inter_reservoir_amount,
        shallow_gw_reservoir_amount, deep_gw_reservoir_amount,
        timedelta,
        actual_evapotranspiration_flux, overland_runoff_flux,
        drain_runoff_flux, inter_runoff_flux, shallow_gw_runoff_flux,
        deep_gw_runoff_flux,
        shallow_gw_flux, deep_gw_flux
    )

    _update_routing(
        overland_runoff_flux, drain_runoff_flux, inter_runoff_flux,
        shallow_gw_runoff_flux, deep_gw_runoff_flux,
        theta_rk,
        river_reservoir_amount,
        timedelta, drainage_area, rho_water,
        river_discharge_flux
    )

def _update_production(
        # inputs
        rainfall_flux,
        potential_evapotranspiration_flux,
        # parameters
        theta_t,
        theta_c,
        theta_h,
        theta_d,
        theta_s,
        theta_z,
        theta_sk,
        theta_fk,
        theta_gk,
        # states
        soil_layers_amounts,
        overland_reservoir_amount,
        drain_reservoir_amount,
        inter_reservoir_amount,
        shallow_gw_reservoir_amount,
        deep_gw_reservoir_amount,
        # constants
        timedelta,
        # outputs
        actual_evapotranspiration_flux,
        overland_runoff_flux,
        drain_runoff_flux,
        inter_runoff_flux,
        shallow_gw_runoff_flux,
        deep_gw_runoff_flux,
        # internals
        shallow_gw_flux,
        deep_gw_flux
):
    # apply parameter T to rainfall data (aerial rainfall correction)
    corrected_rainfall_flux = rainfall_flux * theta_t

    # determine limiting conditions
    rainfall_minus_evapotranspiration_flux = (
            corrected_rainfall_flux - potential_evapotranspiration_flux
    )
    is_energy_limited = rainfall_minus_evapotranspiration_flux > 0.0
    is_water_limited = ~is_energy_limited

    # calculate total antecedent soil moisture
    soil_amount = np.sum(soil_layers_amounts[0, ...], axis=-1)

    # ------------------------------------------------------------------
    # under energy-limited conditions
    # >>> --------------------------------------------------------------

    effective_rainfall_flux = np.where(
        is_energy_limited, rainfall_minus_evapotranspiration_flux, 0.0
    )

    # -------------------------------------------------------------- <<<

    # ------------------------------------------------------------------
    # under water-limited conditions
    # >>> --------------------------------------------------------------

    # ignore cells where there is rain excess
    unmet_evapotranspiration_flux = np.where(
        is_water_limited, -rainfall_minus_evapotranspiration_flux, 0.0
    )

    # provisionally set soil evaporation as total available moisture
    max_soil_evaporation_flux = np.where(
        is_water_limited, soil_amount / timedelta, 0.0
    )

    # limit contribution to unmet ET where there is moisture excess
    soil_evaporation_flux = np.where(
        max_soil_evaporation_flux >= unmet_evapotranspiration_flux,
        unmet_evapotranspiration_flux,
        max_soil_evaporation_flux
    )

    # -------------------------------------------------------------- <<<

    # calculate actual evapotranspiration
    actual_evapotranspiration_flux[...] = np.where(
        is_energy_limited,
        potential_evapotranspiration_flux,
        corrected_rainfall_flux + soil_evaporation_flux
    )

    # determine excess rain amount from effective rainfall flux
    excess_rainfall_amount = effective_rainfall_flux * timedelta

    # initialise current soil layers to their level at previous step
    soil_layers_amounts[1, ...] = soil_layers_amounts[0, ...]

    # ------------------------------------------------------------------
    # under energy-limited conditions
    # >>> --------------------------------------------------------------

    # calculate surface runoff using quick runoff parameter H and
    # relative soil moisture content
    theta_h_prime = theta_h * (soil_amount / theta_z)
    # excess rainfall contribution to quick surface runoff store
    overland_flow = theta_h_prime * excess_rainfall_amount
    # remainder that infiltrates
    excess_rainfall_amount -= overland_flow

    # calculate percolation through soil layers
    # (from top layer [1st] to bottom layer [6th])
    layer_capacity = theta_z / 6.
    for i in range(6):
        layer_level = soil_layers_amounts[1, ..., i]

        # determine space in layer before reaching full capacity
        layer_space = layer_capacity - layer_level

        has_enough_space = excess_rainfall_amount <= layer_space

        # enough space in layer to hold entire excess rain
        layer_level[...] = np.where(
            is_energy_limited & has_enough_space,
            layer_level + excess_rainfall_amount,
            layer_level
        )
        excess_rainfall_amount[...] = np.where(
            is_energy_limited & has_enough_space,
            0.,
            excess_rainfall_amount
        )

        # not enough space in layer to hold entire excess rain
        layer_level[...] = np.where(
            is_energy_limited & ~has_enough_space,
            layer_capacity,
            layer_level
        )
        excess_rainfall_amount[...] = np.where(
            is_energy_limited & ~has_enough_space,
            excess_rainfall_amount - layer_space,
            excess_rainfall_amount
        )

    # calculate saturation excess from remaining excess rainfall
    # sat. excess contrib. (if not 0) to quicker soil runoff store
    drain_flow = theta_d * excess_rainfall_amount
    # sat. excess contrib. (if not 0) to slower soil runoff store
    inter_flux = (1.0 - theta_d) * excess_rainfall_amount

    # -------------------------------------------------------------- <<<

    # calculate leak from soil layers
    # (i.e. piston flow becoming active during rainfall events)
    theta_s_prime = theta_s * (soil_amount / theta_z)

    # calculate soil moisture contributions to runoff stores
    for i in range(6):
        layer_level = soil_layers_amounts[1, ..., i]

        # leak to interflow
        leak_inter_flux = np.where(
            is_energy_limited,
            # soil moisture outflow reducing exponentially downwards
            layer_level * (theta_s_prime ** (i + 1)),
            # no soil moisture contribution to runoff store
            0.
        )
        inter_flux += leak_inter_flux
        layer_level[...] = layer_level - leak_inter_flux

        # leak to shallow groundwater flow
        leak_shallow_gw_flux = np.where(
            is_energy_limited,
            # soil moisture outflow reducing linearly downwards
            layer_level * (theta_s_prime / (i + 1)),
            # no soil moisture contribution to runoff store
            0
        )
        shallow_gw_flux += leak_shallow_gw_flux
        layer_level[...] = layer_level - leak_shallow_gw_flux

        # leak to deep groundwater flow
        leak_deep_gw_flux = np.where(
            is_energy_limited,
            # soil moisture outflow reducing exponentially upwards
            layer_level * (theta_s_prime ** (6 - i)),
            # no soil moisture contribution to runoff store
            0
        )
        deep_gw_flux += leak_deep_gw_flux
        layer_level[...] = layer_level - leak_deep_gw_flux

    # ------------------------------------------------------------------
    # under water-limited conditions
    # >>> --------------------------------------------------------------

    # attempt to satisfy PE from soil layers
    # (from top layer [1st] to bottom layer [6th])
    for i in range(6):
        layer_level = soil_layers_amounts[1, ..., i]

        enough_moisture = unmet_evapotranspiration_flux <= layer_level

        # enough soil moisture in layer
        layer_level[...] = np.where(
            is_water_limited & enough_moisture,
            layer_level - unmet_evapotranspiration_flux,
            layer_level
        )
        unmet_evapotranspiration_flux[...] = np.where(
            is_water_limited & enough_moisture,
            0.,
            unmet_evapotranspiration_flux
        )

        # not enough soil moisture in layer
        layer_level[...] = np.where(
            is_water_limited & ~enough_moisture,
            0.,
            layer_level
        )
        unmet_evapotranspiration_flux[...] = np.where(
            is_water_limited & ~enough_moisture,
            theta_c * (unmet_evapotranspiration_flux - layer_level),
            unmet_evapotranspiration_flux
        )

    # -------------------------------------------------------------- <<<

    # route overland runoff
    overland_runoff_flux[...] = overland_reservoir_amount[0] / theta_sk
    overland_reservoir_amount[1, ...] = (
        overland_reservoir_amount[0] + overland_flow
        - overland_runoff_flux * timedelta
    )
    overland_reservoir_amount[1, ...] *= overland_reservoir_amount[1] > 0

    # route drain runoff
    drain_runoff_flux[...] = drain_reservoir_amount[0] / theta_sk
    drain_reservoir_amount[1, ...] = (
        drain_reservoir_amount[0] + drain_flow
        - drain_runoff_flux * timedelta
    )
    drain_reservoir_amount[1, ...] *= drain_reservoir_amount[1] > 0

    # route inter runoff
    inter_runoff_flux[...] = inter_reservoir_amount[0] / theta_fk
    inter_reservoir_amount[1, ...] = (
        inter_reservoir_amount[0] + inter_flux
        - inter_runoff_flux * timedelta
    )
    inter_reservoir_amount[1, ...] *= inter_reservoir_amount[1] > 0

    # route shallow groundwater runoff
    shallow_gw_runoff_flux[...] = shallow_gw_reservoir_amount[0] / theta_gk
    shallow_gw_reservoir_amount[1, ...] = (
        shallow_gw_reservoir_amount[0] + shallow_gw_flux
        - shallow_gw_runoff_flux * timedelta
    )
    shallow_gw_reservoir_amount[1, ...] *= shallow_gw_reservoir_amount[1] > 0

    # route deep groundwater runoff
    deep_gw_runoff_flux[...] = deep_gw_reservoir_amount[0] / theta_gk
    deep_gw_reservoir_amount[1, ...] = (
        deep_gw_reservoir_amount[0] + deep_gw_flux
        - deep_gw_runoff_flux * timedelta
    )
    deep_gw_reservoir_amount[1, ...] *= deep_gw_reservoir_amount[1] > 0

def _update_routing(
        # inputs
        overland_runoff_flux,
        drain_runoff_flux,
        inter_runoff_flux,
        shallow_gw_runoff_flux,
        deep_gw_runoff_flux,
        # parameters
        theta_rk,
        # states
        river_reservoir_amount,
        # constants
        timedelta,
        drainage_area,
        rho_water,
        # outputs
        river_discharge_flux
):
    total_runoff_flux = (
        overland_runoff_flux + drain_runoff_flux + inter_runoff_flux
        + shallow_gw_runoff_flux + deep_gw_runoff_flux
    )

    # provisionally calculate river flow
    discharge_flux = river_reservoir_amount[0] / theta_rk

    # provisionally calculate new river store state
    reservoir_amount = (
        river_reservoir_amount[0]
        + (total_runoff_flux - discharge_flux) * timedelta
    )

    # check whether store has gone negative
    discharge_flux = np.where(
        reservoir_amount < 0,
        # allow max outflow at 95% of what was in store
        0.95 * (total_runoff_flux + river_reservoir_amount[-1] / timedelta),
        discharge_flux
    )
    river_reservoir_amount[1, ...] = (
        river_reservoir_amount[0]
        + (total_runoff_flux - discharge_flux) * timedelta
    )

    # convert [kg m-2 s-1] to [m3 s-1]
    river_discharge_flux[...] = (
        discharge_flux / rho_water * drainage_area
    )
