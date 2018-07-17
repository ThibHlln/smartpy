[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![PyPI Version](https://badge.fury.io/py/smartpy.svg)](https://pypi.python.org/pypi/smartpy)
[![Travis CI Build Status](https://www.travis-ci.org/ThibHlln/smartpy.svg?branch=master)](https://www.travis-ci.org/ThibHlln/smartpy)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/ThibHlln/smartpy?branch=master&svg=true)](https://ci.appveyor.com/project/ThibHlln/smartpy)

# SMARTpy - An open-source version of the rainfall-runoff model SMART in Python

SMARTpy is an open-source hydrological catchment model in Python. It is licensed under GNU GPL-3.0 (see [licence file](LICENCE.md) provided). SMART (Soil Moisture Accounting and Routing for Transport) is a top-down rainfall-runoff model composed of a soil moisture accounting component and linear routing components. It requires rainfall and potential evapotranspiration time series as inputs, it features a set of ten parameters, and it yields a discharge time series.

Mockler, E., O’Loughlin, F., and Bruen, M.: Understanding hydrological flow paths in conceptual catchment models using uncertainty and sensitivity analysis, *Computers & Geosciences*, 90, 66–77,[doi:10.1016/j.cageo.2015.08.015](https://dx.doi.org/10.1016/j.cageo.2015.08.015), 2016

## How to Install

SMARTpy is available on PyPI, so you can simply use pip:

    python -m pip install smartpy

Alternatively, you can download the source code (*i.e.* this repository) and use the command:

    python setup.py install

## Dependencies

SMARTpy requires the popular Python package `numpy` and `scipy` to be installed on the Python implementation where `smartpy` is installed. For Python 2 and 3 compatibilities, the package `future` is also required.
Additional optional dependencies include `netCDF4` if one wishes to use NetCDF files as input, and `smartcpp` if one wishes to use an accelerator module for the `SMART` model (it gives access to a [C++ extension for the SMART model](https://github.com/ThibHlln/smartcpp)).

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

* 0.1.1 [17 Jul 2018]: Version with proper PyPI display
	* Fixes display issue of README.md on PyPI
* 0.1.0 [17 Jul 2018]: First version of SMARTpy

## Acknowledgment

This tool was developed with the financial support of Ireland's Environmental Protection Agency (Grant Number 2014-W-LS-5).