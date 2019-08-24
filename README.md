[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![PyPI Version](https://badge.fury.io/py/smartpy.svg)](https://pypi.python.org/pypi/smartpy)
[![Travis CI Build Status](https://www.travis-ci.org/ThibHlln/smartpy.svg?branch=master)](https://www.travis-ci.org/ThibHlln/smartpy)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/ThibHlln/smartpy?branch=master&svg=true)](https://ci.appveyor.com/project/ThibHlln/smartpy)
[![DOI](https://zenodo.org/badge/118467753.svg)](https://zenodo.org/badge/latestdoi/118467753)

# SMARTpy - An open-source version of the rainfall-runoff model SMART in Python

SMARTpy is an open-source hydrological catchment model in Python. It is licensed under GNU GPL-3.0 (see [licence file](LICENCE.md) provided). SMART (Soil Moisture Accounting and Routing for Transport) is a top-down rainfall-runoff model composed of a soil moisture accounting component and linear routing components. It requires rainfall and potential evapotranspiration time series as inputs, it features a set of ten parameters, and it yields a discharge time series.

## How to Install

SMARTpy is available on PyPI, so you can simply use pip:

    python -m pip install smartpy

You can also use a link to the GitHub repository directly:

	python -m pip install git+https://github.com/ThibHlln/smartpy.git

Alternatively, you can download the source code (*i.e.* the GitHub repository) and, from the downloaded directory itself, run the command:

    python setup.py install

## How to Use

A tutorial in the form of a [Jupyter notebook](examples/api_usage_example.ipynb) is available to get started with the usage of SMARTpy's API. The input files required for the tutorial are all provided in the `examples/` folder.

## How to Cite

If you are using SMARTpy, please consider citing the software as follows (click on the link to get the DOI of a specific version):
* Hallouin, T., Mockler, E., Bruen, M. (XXXX). SMARTpy: Conceptual Rainfall-Runoff Model (Version X.X.X). Zenodo. https://doi.org/10.5281/zenodo.2564041

## Dependencies

SMARTpy requires the popular Python packages `numpy` and `scipy` to be installed on the Python implementation where `smartpy` is installed. The package `future` is also required (used for Python 2 to 3 compatibilities). Additional optional dependencies include `netCDF4` if one wishes to use NetCDF files as input or output, and `smartcpp` if one wishes to use an accelerator module for the SMART model ([C++ extension for the SMART model](https://github.com/ThibHlln/smartcpp)).

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

### References

Mockler, E., O’Loughlin, F., and Bruen, M.: Understanding hydrological flow paths in conceptual catchment models using uncertainty and sensitivity analysis, *Computers & Geosciences*, 90, 66–77,[doi:10.1016/j.cageo.2015.08.015](https://dx.doi.org/10.1016/j.cageo.2015.08.015), 2016

## Input/Output File Formats

SMARTpy is designed to read CSV (Comma-Separated Values) files and NetCDF (Network Common Data Form) input files (for rain, peva, and flow), as well as to write CSV and NetCDF output files (for discharge series, and monte carlo simulation results). However, the use of NetCDF files requires the Python package `netCDF4` to be installed on the Python implementation where SMARTpy is installed (specific pre-requisites prior the installation of `netCDF4` exist and can be found at [unidata.github.io/netcdf4-python](http://unidata.github.io/netcdf4-python/)).

## Monte Carlo Simulations

The `montecarlo` suite of classes that comes with `smartpy` gives access to various options for Monte Carlo simulations. The parameter space of the SMART model can be explored using Latin Hypercube Sampling (`LHS`) for any sample size required. Once the sampling is complete for on a given simulation period, the performance of the whole set of parameter sets can be evaluated on another simulation period with `Total`, or the set of parameter sets can be conditioned according to their own performances against observed discharge data on any of the objective function(s) calculated by SMARTpy (using `GLUE` to distinguish from behavioural and non-behavioural parameter sets, or using `Best` to retain a pre-defined number of best performing samples) and the resulting subset of parameter sets can be evaluated on another simulation period.

## Parallel Computing

If Monte Carlo simulations are required, it is important to make use of the available computer power to reduce the runtime. Personal Computers now commonly feature several processor cores that can be used to run as many runs of the SMART model in parallel (*i.e.* at the same time), not to mention High Performance Clusters, where the benefits of parallel computing will be even more significant. The `montecarlo` classes of `smartpy` are using the `spotpy` package to give access to an easy way to run simulations in parallel. `spotpy` itself requires `mpi4py` to operate, which applies to `smartpy` by extension. So before using `montecarlo` with `parallel='mpi'`, a Message Passing Interface (MPI) library (*e.g.* Open MPI) and `mpi4py` need to be installed on your machine, and `spotpy` needs to be installed too. Any of the `montecarlo` classes can take an optional argument parallel, its default value is set to 'seq' (for sequential computing), but can be set to 'mpi' if your setup allows it (for parallel computing).

## Version History

* 0.2.1 [24 Aug 2019]: [General enhancements](https://github.com/ThibHlln/smartpy/releases/tag/v0.2.1)
* 0.2.0 [16 Nov 2018]: [Speed improvement by making use of new version of SMARTcpp](https://github.com/ThibHlln/smartpy/releases/tag/v0.2.0)
* 0.1.4 [12 Nov 2018]: [General enhancements](https://github.com/ThibHlln/smartpy/releases/tag/v0.1.4)
* 0.1.3 [24 Jul 2018]: [Version improved for Monte Carlo simulations with parallel computing](https://github.com/ThibHlln/smartpy/releases/tag/v0.1.3)
* 0.1.2 [18 Jul 2018]: [Version functioning without evaluation data](https://github.com/ThibHlln/smartpy/releases/tag/v0.1.2)
* 0.1.1 [17 Jul 2018]: [Version with proper PyPI display](https://github.com/ThibHlln/smartpy/releases/tag/v0.1.1)
* 0.1.0 [17 Jul 2018]: [First version of SMARTpy](https://github.com/ThibHlln/smartpy/releases/tag/v0.1.0)

## Acknowledgment

This tool was developed with the financial support of Ireland's Environmental Protection Agency (Grant Number 2014-W-LS-5).