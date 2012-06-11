==========
core.gates
==========

This is a self-contained module which defines the ideal unitaries, for both the pulse types used in the experiment, and the standard gates used in quantum computing (CNOTs, SWAP etc). It can be use independently to e.g. generate the overall unitary for a short algorithm or pulse sequence. It is used by ``simtools`` to generate ideal unitaries for the pulses used in the full simulation.

basic matrices
--------------

These are the Pauli matrices and similar small building-blocks::

  X = np.array([[0,1],[1,0]])
  Y = np.array([[0,-1j],[1j,0]])
  Z = np.array([[1,0],[0,-1]])
  H = 1/np.sqrt(2)*np.array([[1,1],[1,-1]])
  I = np.array([[1,0],[0,1]])
  S = np.array([[1,0],[0,1j]])
  T = np.array([[1,0],[0,np.exp(1j*np.pi/4)]])

composite gates
---------------

Use these to generate circuits consisting of mostly single- and two-qubit gates.

.. autofunction:: PyTIQC.core.gates.dotN

.. autofunction:: PyTIQC.core.gates.kronN

.. autofunction:: PyTIQC.core.gates.R

.. autofunction:: PyTIQC.core.gates.G

.. autofunction:: PyTIQC.core.gates.control

.. autofunction:: PyTIQC.core.gates.controlMN

.. autofunction:: PyTIQC.core.gates.swap

.. autofunction:: PyTIQC.core.gates.CNOT

.. autofunction:: PyTIQC.core.gates.Toffoli

Two helper functions:

.. autofunction:: PyTIQC.core.gates.dec2bin

.. autofunction:: PyTIQC.core.gates.bin2dec

Gates with arbitrary angles
---------------------------

These are global and individual gates used by the experiment.

.. autofunction:: PyTIQC.core.gates.Rg

.. autofunction:: PyTIQC.core.gates.Rx

.. autofunction:: PyTIQC.core.gates.Ry

.. autofunction:: PyTIQC.core.gates.Ug

.. autofunction:: PyTIQC.core.gates.Ri

.. autofunction:: PyTIQC.core.gates.Ux

.. autofunction:: PyTIQC.core.gates.Uy
.. autofunction:: PyTIQC.core.gates.Uz
.. autofunction:: PyTIQC.core.gates.MS
.. autofunction:: PyTIQC.core.gates.MSx
.. autofunction:: PyTIQC.core.gates.MSy

calculateevolution
------------------

.. autofunction:: PyTIQC.core.gates.calculateevolution

functions for evaluation
------------------------

.. autofunction:: PyTIQC.core.gates.fidelity

.. autofunction:: PyTIQC.core.gates.jozsafid

.. autofunction:: PyTIQC.core.gates.sso

.. autofunction:: PyTIQC.core.gates.tracedist

These can be called directly or by using the module object fidelities_dict, which has the strings ``jozsafid, tracedist-rho, tracedist-pop, sso`` as keys.

.. autodata:: PyTIQC.core.gates.fidelities_dict

plotting and display
--------------------

.. autofunction:: PyTIQC.core.gates.dispmtx
.. autofunction:: PyTIQC.core.gates.displaystates
.. autofunction:: PyTIQC.core.gates.displaytracedstates
.. autofunction:: PyTIQC.core.gates.plotBar3D
