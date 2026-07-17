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
        soil_evaporation_amount, shallow_gw_amount, deep_gw_amount,
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
            soil_evaporation_amount, shallow_gw_amount, deep_gw_amount,
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
        soil_evaporation_amount, shallow_gw_amount, deep_gw_amount,
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
        soil_evaporation_amount, shallow_gw_amount, deep_gw_amount
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
        soil_evaporation_amount,
        shallow_gw_amount,
        deep_gw_amount
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
    soil_total_amount = np.sum(soil_layers_amounts[0, ...], axis=-1)

    # initialise current soil layers to their level at previous step
    soil_layers_amounts[1, ...] = soil_layers_amounts[0, ...]

    # ------------------------------------------------------------------
    # under energy-limited conditions
    # >>> --------------------------------------------------------------
    effective_rainfall_flux = np.where(
        is_energy_limited, rainfall_minus_evapotranspiration_flux, 0.0
    )

    # determine excess rain amount from effective rainfall flux
    excess_rainfall_amount = effective_rainfall_flux * timedelta

    # calculate surface runoff using quick runoff parameter H and
    # relative soil moisture content
    theta_h_prime = theta_h * (soil_total_amount / theta_z)
    # excess rainfall contribution to quick surface runoff store
    overland_amount = theta_h_prime * excess_rainfall_amount
    # remainder that infiltrates
    excess_rainfall_amount -= overland_amount

    # calculate percolation through soil layers
    # (from top layer [1st] to bottom layer [6th])
    layer_max_amount = theta_z / 6.
    for i in range(6):
        layer_amount = soil_layers_amounts[1, ..., i]

        # determine space in layer before reaching full capacity
        layer_free_amount = layer_max_amount - layer_amount

        has_enough_space = excess_rainfall_amount <= layer_free_amount

        # enough space in layer to hold entire excess rain
        layer_amount[...] = np.where(
            is_energy_limited & has_enough_space,
            layer_amount + excess_rainfall_amount,
            layer_amount
        )
        excess_rainfall_amount[...] = np.where(
            is_energy_limited & has_enough_space,
            0.,
            excess_rainfall_amount
        )

        # not enough space in layer to hold entire excess rain
        layer_amount[...] = np.where(
            is_energy_limited & ~has_enough_space,
            layer_max_amount,
            layer_amount
        )
        excess_rainfall_amount[...] = np.where(
            is_energy_limited & ~has_enough_space,
            excess_rainfall_amount - layer_free_amount,
            excess_rainfall_amount
        )

    # calculate saturation excess from remaining excess rainfall
    # sat. excess contrib. (if not 0) to quicker soil runoff store
    drain_amount = theta_d * excess_rainfall_amount
    # sat. excess contrib. (if not 0) to slower soil runoff store
    inter_amount = (1.0 - theta_d) * excess_rainfall_amount

    # calculate leak from soil layers
    # (i.e. piston flow becoming active during rainfall events)
    theta_s_prime = theta_s * (soil_total_amount / theta_z)

    # calculate soil moisture contributions to runoff stores
    for i in range(6):
        layer_amount = soil_layers_amounts[1, ..., i]

        # leak to interflow
        leak_inter_amount = np.where(
            is_energy_limited,
            # soil moisture outflow reducing exponentially downwards
            layer_amount * (theta_s_prime ** (i + 1)),
            # no soil moisture contribution to runoff store
            0.
        )
        inter_amount += leak_inter_amount
        layer_amount[...] = layer_amount - leak_inter_amount

    shallow_gw_amount[:] = 0
    for i in range(6):
        layer_amount = soil_layers_amounts[1, ..., i]

        # leak to shallow groundwater flow
        leak_shallow_gw_amount = np.where(
            is_energy_limited,
            # soil moisture outflow reducing linearly downwards
            layer_amount * (theta_s_prime / (i + 1)),
            # no soil moisture contribution to runoff store
            0
        )
        shallow_gw_amount += leak_shallow_gw_amount
        layer_amount[...] = layer_amount - leak_shallow_gw_amount

    deep_gw_amount[:] = 0
    for i in range(5, -1, -1):
        layer_amount = soil_layers_amounts[1, ..., i]

        # leak to deep groundwater flow
        leak_deep_gw_amount = np.where(
            is_energy_limited,
            # soil moisture outflow reducing exponentially upwards
            layer_amount * (theta_s_prime ** (6 - i)),
            # no soil moisture contribution to runoff store
            0
        )
        deep_gw_amount += leak_deep_gw_amount
        layer_amount[...] = layer_amount - leak_deep_gw_amount

    # -------------------------------------------------------------- <<<

    # ------------------------------------------------------------------
    # under water-limited conditions
    # >>> --------------------------------------------------------------
    unmet_evapotranspiration_flux = np.where(
        is_water_limited, -rainfall_minus_evapotranspiration_flux, 0.0
    )

    # attempt to satisfy PE from soil layers
    # (from top layer [1st] to bottom layer [6th])
    soil_evaporation_amount[:] = 0
    for i in range(6):
        layer_amount = soil_layers_amounts[1, ..., i]

        enough_soil_moisture = (
            unmet_evapotranspiration_flux * timedelta <= layer_amount
        )

        # enough soil moisture in layer
        layer_amount[...] = np.where(
            is_water_limited & enough_soil_moisture,
            layer_amount - unmet_evapotranspiration_flux * timedelta,
            layer_amount
        )
        soil_evaporation_amount[...] = np.where(
            is_water_limited & enough_soil_moisture,
            soil_evaporation_amount
            + unmet_evapotranspiration_flux * timedelta,
            soil_evaporation_amount
        )
        unmet_evapotranspiration_flux[...] = np.where(
            is_water_limited & enough_soil_moisture,
            0.,
            unmet_evapotranspiration_flux
        )

        # not enough soil moisture in layer
        soil_evaporation_amount[...] = np.where(
            is_water_limited & ~enough_soil_moisture,
            soil_evaporation_amount + layer_amount / timedelta,
            soil_evaporation_amount
        )
        unmet_evapotranspiration_flux[...] = np.where(
            is_water_limited & ~enough_soil_moisture,
            theta_c * (
                unmet_evapotranspiration_flux - layer_amount / timedelta
            ),
            unmet_evapotranspiration_flux
        )
        layer_amount[...] = np.where(
            is_water_limited & ~enough_soil_moisture,
            0.,
            layer_amount
        )

    # -------------------------------------------------------------- <<<

    # calculate actual evapotranspiration
    actual_evapotranspiration_flux[...] = np.where(
        is_energy_limited,
        potential_evapotranspiration_flux,
        corrected_rainfall_flux + soil_evaporation_amount * timedelta
    )

    # route overland runoff
    overland_runoff_flux[...] = overland_reservoir_amount[0] / theta_sk
    overland_reservoir_amount[1, ...] = (
        overland_reservoir_amount[0] + overland_amount
        - overland_runoff_flux * timedelta
    )
    overland_reservoir_amount[1, ...] *= overland_reservoir_amount[1] > 0

    # route drain runoff
    drain_runoff_flux[...] = drain_reservoir_amount[0] / theta_sk
    drain_reservoir_amount[1, ...] = (
        drain_reservoir_amount[0] + drain_amount
        - drain_runoff_flux * timedelta
    )
    drain_reservoir_amount[1, ...] *= drain_reservoir_amount[1] > 0

    # route inter runoff
    inter_runoff_flux[...] = inter_reservoir_amount[0] / theta_fk
    inter_reservoir_amount[1, ...] = (
        inter_reservoir_amount[0] + inter_amount
        - inter_runoff_flux * timedelta
    )
    inter_reservoir_amount[1, ...] *= inter_reservoir_amount[1] > 0

    # route shallow groundwater runoff
    shallow_gw_runoff_flux[...] = shallow_gw_reservoir_amount[0] / theta_gk
    shallow_gw_reservoir_amount[1, ...] = (
            shallow_gw_reservoir_amount[0] + shallow_gw_amount
            - shallow_gw_runoff_flux * timedelta
    )
    shallow_gw_reservoir_amount[1, ...] *= shallow_gw_reservoir_amount[1] > 0

    # route deep groundwater runoff
    deep_gw_runoff_flux[...] = deep_gw_reservoir_amount[0] / theta_gk
    deep_gw_reservoir_amount[1, ...] = (
            deep_gw_reservoir_amount[0] + deep_gw_amount
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
