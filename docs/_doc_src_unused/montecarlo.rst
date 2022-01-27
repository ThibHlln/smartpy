.. currentmodule:: smartpy
.. default-role:: obj

Monte Carlo Experiments
=======================

The `montecarlo` suite of classes that comes with `smartpy` gives access
to various options for Monte Carlo simulations. The parameter space of
the SMART model can be explored using Latin Hypercube Sampling (`LHS`)
for any sample size required. Once the sampling is complete for on a
given simulation period, the performance of the whole set of parameter
sets can be evaluated on another simulation period with `Total`, or the
set of parameter sets can be conditioned according to their own
performances against observed discharge data on any of the objective
function(s) calculated by `smartpy` (using `GLUE` to distinguish from
behavioural and non-behavioural parameter sets, or using `Best` to
retain a pre-defined number of best performing samples) and the
resulting subset of parameter sets can be evaluated on another
simulation period.
