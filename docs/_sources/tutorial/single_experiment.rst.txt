.. currentmodule:: smartpy
.. default-role:: obj

Single experiment
=================

Setting up
----------

After having installed `smartpy` (following the
:doc:`installation <installation>` instructions), the first thing that
needs to be done is to import the `smartpy` package to start using its
modules, classes, and functions.

.. code-block:: python

   import smartpy


Creating an instance of the SMART class
'''''''''''''''''''''''''''''''''''''''

Then, you need to create an instance of the class SMART for the catchment
you desire to model. By doing so, all the input data will be collected
from the input files. All arguments are mandatory except *gauged_area_m2*.
If *gauged_area_m2* is provided, it means that discharge measurements are
available, and `smartpy` is expecting a *.flow* file in the input folder.
These measurements can be used later on in post-processing to evaluate
the performance of the simulation.

.. code-block:: python

   from datetime import datetime, timedelta

   sm = smartpy.SMART(
       catchment='ExampleDaily',
       catchment_area_m2=175.46e6,
       start=datetime.strptime('01/01/2007 09:00:00', '%d/%m/%Y %H:%M:%S'),
       end=datetime.strptime('31/12/2016 09:00:00', '%d/%m/%Y %H:%M:%S'),
       time_delta_simu=timedelta(hours=1),
       time_delta_save=timedelta(days=1),
       warm_up_days=365,
       in_format='csv',
       out_format='csv',
       root='examples/',
       gauged_area_m2=175.97e6
   )

.. seealso:: `smartpy.SMART`


Giving parameter values to the model instance
'''''''''''''''''''''''''''''''''''''''''''''

Now it is time to give the model its parameter values. The recommended
way to go is to use a file containing the parameter values, following
the template in the example file used for this tutorial. The given values
will be accessible in the attribute `{SMART_instance}.parameters.values`
(where *SMART_instance* is `sm` in this tutorial, that it is to say the
name binding to the instance of SMART used).


.. code-block:: python

   sm.parameters.set_parameters_with_file('examples/in/ExampleDaily/ExampleDaily.parameters')


Alternatively, you can avoid the creation of a *.parameters* file and
directly provide them as a dictionary, as follows:

.. code-block:: python

   sm.parameters.set_parameters_with_dict(
       {
           'T': 1.0,
           'C': 1.0,
           'H': 0.20845,
           'D': 0.24606,
           'S': 0.0001230,
           'Z': 105.26,
           'SK': 46.82,
           'FK': 315.55,
           'GK': 1066.73,
           'RK': 10.64
       }
   )

.. seealso:: `smartpy.parameters.Parameters.set_parameters_with_file`,
             `smartpy.parameters.Parameters.set_parameters_with_dict`


Defining the initial conditions
'''''''''''''''''''''''''''''''

If not initialised, the model states (i.e. levels of soil moisture layers,
and levels of linear reservoirs) are set to zero (i.e. empty), which will
give poor streamflow predictions for the beginning of your simulation
period. It is good practice to use a warm-up period to initialise the model
states. This is readily available with `smartpy`, and we have already defined
a duration for it in this tutorial (i.e. 365 days) when creating our SMART
instance with its *warm_up_days* parameter. If this duration is not set
to zero, `smartpy` will use the first X days (here the first 365 days) of
the model input time series (i.e. precipitation and potential
evapotranspiration) to run the model for the given duration prior the
*start* date of the simulation in order to start the actual simulation
with model layers and reservoirs that are not empty (the last levels of
the warm-up period are used as the assumed initial levels for the
simulation period).

If one does not want to proceed this way, it is possible to set the warm-up
period to 0 days. In this case, the model states are set to zero, and it is
suggested that the user includes the warm-up period as part of the simulation
period, and discards this period in the model outputs afterwards.

In addition to the warm-up period, `smartpy` has an additional feature to
avoid starting with empty layers and reservoirs. This one is usable if the
user has some general information about the hydrology in their catchment:

* *aar*: an estimate of the annual average rainfall in millimetres,
* *r-o_ratio*: an estimate of the runoff ratio,
* *r-o_split*: an estimate of the runoff split between overland flow/drain
  flow/interflow/shallow groundwater flow/deep ground water flow.

If this is the case, `smartpy` will estimate the model states from this
information to make an "educated guess" of the levels of the soil layers
and the routing reservoirs. The information can be provided to the model
as follows:


.. code-block:: python

   sm.extra = {
       'aar': 1200,
       'r-o_ratio': 0.45,
       'r-o_split': (0.10, 0.15, 0.15, 0.30, 0.30)
   }

.. warning::

   This option is more advanced. If the user has access to accurate estimates
   for these values, they may prefer to use this option as initialisation
   method in place of a warm-up period, or in addition to the warm-up period.

Simulating
----------

The SMART instance has now all it needs to start the simulation for the
catchment over the time period set up. One last instruction and the
simulation will be under way. The method simulate needs to be given
the parameter values, they are directly accessible from the SMART
instance attributes.

.. code-block:: python

   sm.simulate(sm.parameters.values)

.. seealso:: `smartpy.SMART.simulate`


Retrieving the outputs
----------------------

The outputs provided by `smartpy` are the simulated streamflow series and
the observed streamflow series. The latter is provided if and only if the
*gauged_area_m2* was provided in the SMART instance. If this is the case,
the *.flow* file provided in the input folder has already been read and
the data is directly available from the SMART instance (see below).

.. warning::

   In cases where *catchment_area_m2* and *gauged_area_m2* are different,
   the observed streamflow series provided in the the *.flow* file is
   rescaled accordingly assuming proportionality.


Writing the outputs to file
---------------------------

If one wants to save the model outputs in files for later use. This can
be done using the following instruction:

.. code-block:: python

   sm.write_output_files(which='both')

.. seealso:: `smartpy.SMART.write_output_files`

The *which* argument can take the value 'modelled', 'observed', or 'both',
to write the simulated streamflow series, the observed streamflow series,
or both series, respectively. The format of the files used is the one
provided in the *out_format* attribute of the SMART instance (either 'csv'
or 'netcdf').

Retrieving the outputs as arrays
--------------------------------

If one wants to directly get the model outputs in the current Python session,
two getter methods of the SMART instance are available to do so:

.. code-block:: python

   evaluation = sm.get_evaluation_array()
   simulation = sm.get_simulation_array()


Post-processing (example)
-------------------------

The evaluation and simulation arrays returned have the same length.
However, it is frequent to have days with missing discharge measurements,
that would have been absent from the *.flow* file provided in the input
folder. `smartpy` assigned a `numpy.nan` (not a number) value for any day
with missing data, this is why the two arrays have the same length in the
model outputs. However, when comparing them using an objective function,
only days with available observations can be used. To deal with this, a
suggested solution is as follows:

.. code-block:: python

   import numpy as np

   mask = ~np.isnan(evaluation)
   evaluation = evaluation[mask]
   simulation = simulation[mask]


Now that the arrays have no missing value, but still the same length,
they can be compared using any objective function. A suggested approach
is to use the `hydroeval` Python package (`<https://thibhlln.github.io/hydroeval>`_)
to get access to some of the most commonly used objective functions in hydrology:

.. code-block:: python

   import hydroeval

   # calculate the Nash-Sutcliffe efficiency
   nse = hydroeval.evaluator(hydroeval.nse, simulation, evaluation)

   # calculate the original Kling-Gupta efficiency
   kge, r, alpha, beta = hydroeval.evaluator(hydroeval.kge, simulation, evaluation)


.. note::

   This example uses made-up data, this is why the efficiency of the model
   simulation is low, e.g. with an NSE value is of 0.39.
