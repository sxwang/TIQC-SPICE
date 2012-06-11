==========
Tutorial
==========

Getting started
---------------

TIQC-SPICE was developed and tested using the python IDE IPython_. It uses the library matplotlib_ for plotting. Both IPython and matplotlib can be installed via apt on ubuntu/debian-like systems. 

.. _IPython: http://ipython.org/
.. _matplotlib: http://matplotlib.sourceforge.net/

TIQC-SPICE modules and submodules can be imported if a python path is set to locate the PyTIQC directory. At the minimum, the following two modules need to be imported::

    import PyTIQC.simtools as sim
    import PyTIQC.qctools as qc


Conventions
-----------

All times are in units of microseconds and frequencies are in units of MHz (including the factor of 2pi. 

For example, a secular frequency of 1 MHz is stored in the parameter ``omz`` as ``omz = 2*pi``.

Ion ordering: qubits are labeled with 0-indexing. In a standard quantum circuit drawing, the top qubit is labeled 0. In a qubit string in TIQC-SPICE, ions are counted from right to left. E.g. in a string "DSSSS", ion 4 is in D. 

For labelling states, 1 = S, 0 = D. In a state vector (represented by a numpy complex array), the first element thus corresponds to the state ``|DD..D>`` and the last element is ``|SS..S>``. For example, the ground state of a 2-qubit system is ``np.array([0, 0, 0, 1])``.


Simulations
-----------

Running a simulation with TIQC-SPICE involves the following major steps:

* define the hilbert space object ``hspace``, the parameters object ``parameters``, and the decoherence object ``decoherence``
* define a pulse sequence ``PulseSequence``
* run the simulation
* perform desired further processing using the returned ``database`` object

Each of these steps will be described in detail in this tutorial.

.. toctree::
   :maxdepth: 0

   hspace
   parameters
   decoherence
   pulseseq
   simulateevolution
   database

