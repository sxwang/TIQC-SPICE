=========================
To-do list for TIQC-SPICE
=========================

For developers: please contribute!

Speed
-----

* core computation speedups

  * sparse matrices
  * save H or U between steps
  * make Hamiltonian (dynamically) local (subsystems only)
  * dynamically resize phonon space
  * provide different solvers for timedependent hamiltonians -> figure out the best/fastest one

* parallelization

  * switch to a better managed/maintained tool: Condor
  * store/manage data in an SQL DB and enable interrupted simulations

development
-----------

* feature/design issues

  * add Rred (red sideband pulse)
  * densitymatrixformalism (master equation approach to solving SE) 
  * load params and pulseseq from file
  * parsers for interfacing with pulse compiler

* testsuite:

  * tests for bichromatic

    * MS(0.5) should give ghz state for even number of qubits (requires additional 0.5 rotation to map into ghz when it's an odd number)

  * tests for decoherence - add heating rate, etc

* state/process tomographies module InvestigatePulseSeq

  * provide both addressed as well as global approach for tomographies
  * include alternative methods of doing partial tomographies (compressed, verification ala hofmann, ...)
  * also provide tools for evaluation

other error sources
-------------------

* quadrupole shift
* spontaneous decay during measurement <-- include in run script
* wrongly set initial params
* long-term drift
* separate dephase_magfield and dephase_laser <-- only needed if levels>2

optional features
-----------------

* for larger lamb-dicke parameter: do the calculations in position/momentum space using a grid of points in that phase space. Instead of |elec state, phonon>, use |elec state, (x,k)>
