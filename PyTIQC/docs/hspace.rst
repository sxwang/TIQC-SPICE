======
Hspace
======

A hilbert space object contains the variables and operators that are generically used to specify an ion-phonon-laser system. An ``hspace`` object is initialized with the following arguments::

    hspace = sim.hspace(nuions=1, levels=2, maxphonons=0, densmat=0)

where:

* ``nuion`` = number of ions
* ``levels`` = atomic levels per ion. Fully developed for 2-level system, but some operators are also defined for 3-level systems.
* ``maxphonons`` = max number of phonons. The choice of this variable depends on the type of simulation desired:

 - for simulations involving thermal states, ``maxphonons`` should be chosen large enough to fit the thermal state (this is not explicitly checked.)
 - likewise for simulations of heating rate. A warning will be printed if the heating rate is set high enough that the phonon space would be exceeded.
 - for the Molmer-Sorensen gate on 2 ions, ``maxphonons=7`` is needed to get 3 digits of precision on the populations of the resulting GHZ state. 
 - otherwise, for most simulations focusing on the qubit states, the ``maxphonons`` can be set to 0.

* ``densmat`` = use density matrix formalism. NOT IMPLEMENTED.

The hilbert space dimension is calculated from these variables as::

    hspace.dimensions = (maxphonons+1) * levels ** nuions

Initial state
-------------

The ``hspace`` object has a function ``initial_state`` which can be used to specify an initial state for the simulation. The default state is the ground state, ``|SS..S,0>``. 

A thermal state can be set as the initial state::

    hspace.initial_state('thermal', nBar=0.5)

where ``nBar`` is the mean thermal state occupation in a normalized poisson distribution. 

A quantum state can be specified as follows. In this example, there are 5 ions and one of them is set to the state D::

    hspace.initial_state('quantum', qstate="DSSSS")

Quantum states are represented by an array of complex values, where the length of the array is ``hspace.dimensions``. The initial state is stored in the variable ``hspace.y0``.
