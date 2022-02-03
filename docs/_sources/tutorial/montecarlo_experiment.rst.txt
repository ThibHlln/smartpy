.. currentmodule:: smartpy
.. default-role:: obj

Monte Carlo experiment
======================

The Monte Carlo experiment presented here consists of two steps: first the
:ref:`sampling <sampling>` in the SMART model parameter space, then the
:ref:`conditioning <conditioning>` of the sample to select those parameter
sets deemed as satisfactorily reproducing the observed river discharge.

In `smartpy`, the sampling generates at random a number of parameter sets
that are used to run the SMART model on a given simulation period (the
"calibration" period). The performance of the SMART model structure with
each parameter set is assessed using one or more objective functions.
These objective function values (and optionally the corresponding
streamflow simulations) correspond to the sampling output.

Then, in `smartpy`, the conditioning selects a subset of parameter sets
from the sampling based on their objective functions values. This subset
of parameter sets is then used to run the SMART model on another given
simulation period (the "evaluation" period, necessarily independent from
the sampling period). The performance of the SMART model with this subset
of parameter sets is then assessed again using one or more objective
functions (not necessarily the same ones as for the conditioning). These
objective function values (and optionally the corresponding streamflow
simulations) correspond to the conditioning output.

.. _sampling:

Sampling
--------

`smartpy` only offers a Latin Hypercube approach (`McKay et al., 2000
<https:doi.org/10.1080/00401706.2000.10485979>`_) for sampling the SMART
model parameter space.

Setting up
''''''''''

The configuration of the sampling is done as follows:

.. code-block:: python

   from smartpy import montecarlo

   lhs = montecarlo.LHS(
       catchment='ExampleDaily',
       root_f='examples/',
       in_format='csv,
       out_format='csv,
       sample_size=10,
       parallel=False,
       save_sim=True,
       settings_filename='ExampleDaily.sttngs'
   )

.. seealso:: `smartpy.montecarlo.LHS`

.. note::

   The configuration of the SMART model itself, unlike for the single
   experiment, is done through a settings file. An example of such a
   settings file is available as a template in the *examples/* directory.

Adjusting the parameter space
`````````````````````````````

The default model parameter ranges used to define the boundaries of the
Latin hypercube are accessible as follows:

.. code-block:: python

   >>> print(lhs.model.parameters.ranges)
   {'T': (0.9, 1.1), 'C': (0.0, 1.0), 'H': (0.0, 0.3), 'D': (0.0, 1.0), 'S': (0.0, 0.013), 'Z': (15.0, 150.0), 'SK': (1.0, 240.0), 'FK': (48.0, 1440.0), 'GK': (1200.0, 4800.0), 'RK': (1.0, 96.0)}

And they can be adjusted by assigning a new dictionary to this attribute:

.. code-block:: python

   lhs.parameter.ranges = {
       'T': (0.9, 1.1),
       'C': (0.0, 1.0),
       'H': (0.0, 0.3),
       'D': (0.0, 1.0),
       'S': (0.0, 0.013),
       'Z': (15.0, 150.0),
       'SK': (1.0, 240.0),
       'FK': (48.0, 1440.0),
       'GK': (1200.0, 4800.0),
       'RK': (1.0, 96.0)
   }


Defining the initial conditions
```````````````````````````````

Like in the single experiment example, the "educated guess" approach can be
used as follows:

.. code-block:: python

   lhs.model.extra = {
       'aar': 1200,
       'r-o_ratio': 0.45,
       'r-o_split': (0.10, 0.15, 0.15, 0.30, 0.30)
   }


Also, a warm up period can be chosen in the settings file where the model
is configured.

Simulating
''''''''''

Simulating in serial
````````````````````

Running the sampling in the parameter space sequentially (i.e. simulating
with one parameter set at a time using one process) can be done as follows:

.. code-block:: python

   lhs.run(compression=False)


.. seealso:: `smartpy.montecarlo.LHS.run`

Simulating in parallel
``````````````````````

To start the sampling in parallel (i.e. simulating with multiple parameter
sets at a time using multiple processes, the parameter *parallel* must be
set as `True` in the set up step above, and all the steps presented above
(including the `lhs.run(...)`) must be saved in a Python script, say
*lhs_sampling.py*. Then, this script must be executed via a `mpirun`
command as follows (e.g. using 4 processors to initialise 4 processes):

.. code-block:: bash

   mpirun -np 4 python lhs_sampling.py

Outputs
```````

The sampling performed above will have produced an output file
*{root_f}/out/{catchment}.SMART.lhs* (if `'csv'` was chosen as the output
format or *{root_f}/out/{catchment}.SMART.lhs.nc* (if `'netcdf'` was
chosen as the output format).

.. _conditioning:

Conditioning
------------

.. warning::

   A sampling output file is required to perform the conditioning step.

`smartpy` offers three conditioning approaches which differ in the way
the parameter sets are elicited as satisfactory:

- `GLUE` is based on the selection of one or more objective functions
  and the a priori definition of conditions (i.e. thresholds) to achieve
  with the model simulation on each objective function to deem the
  parameter set as "behavioural";
- `Best` is based on the selection of one target objective function
  and the a priori choice of the number of "best" performing parameter
  sets to retain;
- `Total` is not a conditioning approach per se, as it retains all
  parameter sets for the simulation on the "evaluation" period. This is
  provided as a convenient way to offer to the user the opportunity
  to perform their own conditioning approach a posteriori using the
  streamflow simulations on the "calibration" and the "evaluation" periods
  obtained from the sampling and conditioning steps, respectively.

.. seealso:: `smartpy.montecarlo.GLUE`, `smartpy.montecarlo.Best`,
             `smartpy.montecarlo.Total`

The remainder of this tutorial is based on the `Best` approach.

Setting up
''''''''''

The configuration of the conditioning is done as follows:

.. code-block:: python

   from smartpy import montecarlo

   montecarlo.Best(
       catchment='ExampleDaily',
       root_f='examples/',
       in_format='csv',
       out_format='csv',
       target='KGE,
       nb_best='2',
       parallel=False,
       save_sim=True
   )

For example, here the two best performing parameter sets amongst the
ten sets sampled are selected based on their Kling-Gupta Efficiency (KGE)
values.

.. important::

   The *catchment* and *root_f* must be the same between the sampling
   and the conditioning steps. In addition, the *out_format* of the
   sampling step must be the sme as the *in_format* of the conditioning
   step.

Defining the initial conditions
```````````````````````````````

The initial conditions can be defined the same way they are in the
sampling step above.

Simulating
''''''''''

The simulations can be performed the same way they are in the sampling
step above.

Outputs
```````

The conditioning performed above will have produced an output file
*{root_f}/out/{catchment}.SMART.best* (if `'csv'` was chosen as the output
format or *{root_f}/out/{catchment}.SMART.best.nc* (if `'netcdf'` was
chosen as the output format).
