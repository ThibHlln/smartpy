.. currentmodule:: smartpy
.. default-role:: obj

Parallel Computing
==================

If Monte Carlo simulations are required, it is important to make use of
the available computer power to reduce the runtime. Personal Computers
now commonly feature several processor cores that can be used to run as
many runs of the SMART model in parallel (*i.e.* at the same time), not
to mention High Performance Clusters, where the benefits of parallel
computing will be even more significant. The `montecarlo` classes of
`smartpy` are using the `spotpy` package to give access to an easy way
to run simulations in parallel. `spotpy` itself requires `mpi4py` to
operate, which applies to `smartpy` by extension. So before using
`montecarlo` with `parallel='mpi'`, a Message Passing Interface (MPI)
library (*e.g.* Open MPI) and `mpi4py` need to be installed on your
machine, and `spotpy` needs to be installed too. Any of the `montecarlo`
classes can take an optional argument parallel, its default value is set
to 'seq' (for sequential computing), but can be set to 'mpi' if your
setup allows it (for parallel computing).
