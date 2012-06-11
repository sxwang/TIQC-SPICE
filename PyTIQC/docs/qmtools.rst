============
core.qmtools
============

This module is used internally by :doc:`qctools` and contains the core objects, ``Hamiltonian`` and ``Noise``.

Hamiltonian
-----------

The time-independent Hamiltonian is constructed by calculating the coupling terms and energies separately. Coupling terms are defined element-wise (for historical reasons. Question: is this faster/slower than constructing it with ``hspace.operator_dict``?)

The time-dependent Hamiltonian is constructed with operators in ``hspace.operator_dict``.

.. autoclass:: PyTIQC.core.qmtools.Hamilton
   :members:

Noise
-----

This object contains a "noise dictionary". The keys to the dictionary are the noise types we want to model. The values are a list of three elements: a multiplicative part (a scalar representing e.g. intensity fluctuation), an additive part (an additional matrix with the same dimensions as H), and a projective part (a matrix with the same dimensions as H, but it's a non-unitary operator e.g. spontaneous decay). Any of these values/matrices may actually be a call to a class method which returns the appropriate matrix when given a time as input.

.. autoclass:: PyTIQC.core.qmtools.Noise
   :members:

other helper functions
----------------------

.. autofunction:: PyTIQC.core.qmtools.indexToState

.. autofunction:: PyTIQC.core.qmtools.stateToIndex

.. autofunction:: PyTIQC.core.qmtools.indexToExcitations

.. autofunction:: PyTIQC.core.qmtools.SEsolver
