.. currentmodule:: smartpy
.. default-role:: obj

Installation
============

If you wish to install the most recent stable version of `smartpy`,
it is available on the Python Package Index (PyPI), simply run:

.. code-block:: bash

   pip install smartpy

If you need the latest, potentially unstable, features listed in the
:doc:`change log <changelog>`, please use the *dev* branch on the
GitHub repository:

.. code-block:: bash

   pip install git+https://github.com/ThibHlln/smartpy.git@dev


.. rubric:: Requirements

The following packages are required to use `smartpy`:

.. literalinclude:: ../../requirements.txt
   :language: none

.. rubric:: Optional dependency for input/output

`smartpy` is designed to read/write CSV (Comma-Separated Values) files
and NetCDF (Network Common Data Form) input/output files. However, only
CSV files can natively be read/written in Python.

To read/write NetCDF files, the `netCDF4` python package is required (specific
pre-requisites prior the installation of `netCDF4` exist and can be found
at `<http://unidata.github.io/netcdf4-python/>`_).

.. rubric:: Optional dependencies for Monte Carlo experiments

The module `smartpy.montecarlo` requires the `spotpy` Python package
to perform the sampling required for Monte Carlo experiments
(`<https://spotpy.readthedocs.io/>`_).

In addition, if this sampling is to be done in parallel, the `mpi4py`
Python package is also required (`<https://mpi4py.readthedocs.io/>`_)

.. rubric:: Optional dependency for performance improvement

`smartpy` features a separate accelerator extension written in C++,
`smartcpp` (`<https://github.com/ThibHlln/smartcpp>`_) that will be
automatically used if it is installed alongside `smartpy`. Using
`smartcpp` can significantly improve the execution time.
