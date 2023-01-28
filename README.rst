An implementation of the rainfall-runoff model SMART in Python
==============================================================

.. image:: https://img.shields.io/pypi/v/smartpy?style=flat-square
   :target: https://pypi.python.org/pypi/smartpy
   :alt: PyPI version
.. image:: https://img.shields.io/badge/dynamic/json?url=https://zenodo.org/api/records/2564041&label=doi&query=doi&style=flat-square
   :target: https://zenodo.org/badge/latestdoi/118467753
   :alt: DOI
.. image:: https://img.shields.io/badge/License-GPL%20v3-green.svg?style=flat-square
   :target: https://www.gnu.org/licenses/gpl-3.0
   :alt: Licence GPL-3.0
.. image:: https://img.shields.io/github/actions/workflow/status/thibhlln/smartpy/tests.yml?style=flat-square&label=tests
   :target: https://github.com/ThibHlln/smartpy/actions/workflows/tests.yml
   :alt: Tests Status

`smartpy` is an open-source Python version of the catchment model SMART
(Soil Moisture Accounting and Routing for Transport). SMART is a bucket-style
rainfall-runoff model composed of a soil moisture accounting component and
linear routing components. It requires rainfall and potential evapotranspiration
time series as input, it features a set of ten parameters, and it yields
a discharge time series as output.
