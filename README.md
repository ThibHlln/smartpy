# SMARTpy - The official open-source of the Soil Moisture Accounting and Routing for Transport (SMART) rainfall-runoff model

SMARTpy is an open-source hydrological catchment model in Python. It is licensed under GNU GPL-3.0 (see license file provided). SMART is a top-down rainfall-runoff model composed of a soil moisture accounting component and a routing component. It requires rainfall and potential evapotranspiration time series, it features a set of ten parameters, and it yields a discharge time series.

## Model Specifications

### Model Inputs

* aerial rainfall time series [mm/time step]
* aerial potential evapotranspiration time series [mm/time step]

### Model Parameters

* T: rainfall aerial correction coefficient [-]
* C: evaporation decay parameter [-]
* H: quick runoff coefficient [-]
* D: drain flow parameter - fraction of saturation excess diverted to drain flow [-]
* S: soil outflow coefficient [-]
* Z: effective soil depth [mm]
* SK: surface routing parameter [hours]
* FK: inter flow routing parameter [hours]
* GK: groundwater routing parameter [hours]
* RK: river channel routing parameter [hours]

### Model Outputs

* discharge time series at catchment outlet [m<sup>3</sup>/s]

## Version History

* 1.0.0 [01 Jul 2018]: First version of SMARTpy