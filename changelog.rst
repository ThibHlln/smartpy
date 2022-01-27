.. default-role:: obj

latest
------

Yet to be versioned and released. Only available from *dev* branch until then.

.. rubric:: Big fixes

* fix bug in output writing for `montecarlo` due to change in behaviour in `spotpy`
  (`e2f695b <https://github.com/ThibHlln/smartpy/commit/e2f695baa1634a5e371cfe1ccc7705709660a97f>`_)
* fix mistake in routing parameter used for "educated guess" as initial conditions
  (`5345b1d <https://github.com/ThibHlln/smartpy/commit/5345b1df012f23e883fe48130fa29f8e991353be>`_)

.. rubric:: Documentation

* add documentation build with `sphinx`
  (`#3 <https://github.com/thibhlln/smartpy/pull/3>`_)


v0.2.1
------

Released on 2018-08-24.

.. rubric:: Algorithm

* check for the duration of the warm-up period against the simulation period
  (to avoid the former being greater than the latter)

.. rubric:: Bug fixes

* change test on existence of observed streamflow data in `get_evaluation_array`
  method (which was raising a ValueError since v0.2.0)

.. rubric:: Documentation fixes

* update tutorial notebook that was not up-to-date since v0.2.0 and that
  contained errors
* update readme file to remove mention of objective functions that do not
  exist anymore since v0.2.0


v0.2.0
------

Released on 2018-11-16.

.. rubric:: Algorithm

* redesign of `run` and `run_all_steps` to allow for the use of the
  version 0.2.0 of `smartcpp` where a method `allsteps` makes it possible
  to go through the whole simulation time series in C++

v0.1.4
------

Released on 2018-11-12.

.. rubric:: Functionality

* make gauged area optional (default to catchment area)
* allow choice on output files to write out
* indicate `smartpy` version in output files

.. rubric:: Algorithm

* use more versatile way to convert Epoch timestamps to DateTime objects

.. rubric:: Bug fixes

* use parallel IO to write NetCDF output files in parallel Monte Carlo simulations
* move netCDF compression to the end of simulation (for netCDF4 v1.4.2 compatibility)

v0.1.3
------

Released on 2018-07-24.

.. rubric:: Algorithm

* add built-in Latin hypercube sampling functionality for Monte Carlo experiments

.. rubric:: Functionality

* add a class `Total` to run existing sampling on an evaluation period
* allow use of netCDF files for outputs (for both single and Monte Carlo runs)
* allow use of netCDF files for sampling (for Monte Carlo runs)
* allow use of netCDF files for evaluation data
* add optional compression when writing CSV files

v0.1.2
------

Released on 2018-07-18.

.. rubric:: Functionality

* allow use of `smartpy` without the need for evaluation data
* add possibility to set parameter values using a dictionary

v0.1.1
------

Released on 2018-07-17.

.. rubric:: Bug fixes

* fix improper display of README on PyPI due to use of outdated `setuptools`,
  `twine`, and `wheel` when building the distribution

v0.1.0
------

Released on 2018-07-17.

* first release
