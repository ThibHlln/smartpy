An implementation of the rainfall-runoff model SMART in Python
==============================================================

.. image:: https://img.shields.io/pypi/v/smartpy?style=flat-square
   :target: https://pypi.python.org/pypi/smartpy
   :alt: PyPI version
.. image:: https://img.shields.io/badge/dynamic/json?url=https://zenodo.org/api/records/2564041&label=doi&query=doi&style=flat-square
   :target: https://zenodo.org/badge/latestdoi/118467753
   :alt: DOI
.. image:: https://img.shields.io/badge/License-GPL%20v3-blue.svg?style=flat-square
   :target: https://www.gnu.org/licenses/gpl-3.0
   :alt: Licence GPL-3.0
.. image:: https://www.travis-ci.org/ThibHlln/smartpy.svg?branch=master&style=flat-square
   :target: https://www.travis-ci.org/ThibHlln/smartpy
   :alt: Test status linux+macos
.. image:: https://ci.appveyor.com/api/projects/status/github/ThibHlln/smartpy?branch=master&svg=true
   :target: https://ci.appveyor.com/project/ThibHlln/smartpy
   :alt: Test status windows

`smartpy` is an open-source hydrological catchment model in Python. It is
licensed under GNU GPL-3.0. SMART (Soil Moisture Accounting and Routing
for Transport) is a bucket-style rainfall-runoff model composed of a
soil moisture accounting component and linear routing components. It
requires rainfall and potential evapotranspiration time series as inputs,
it features a set of ten parameters, and it yields a discharge time series.