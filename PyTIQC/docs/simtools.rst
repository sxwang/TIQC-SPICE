=============
core.simtools
=============

hspace
------

.. autoclass:: PyTIQC.core.simtools.hspace
   :members:

parameters
----------

This class contains a large number of parameters related to the experiment, and a few functions to manipulate them and calculate the inter-dependent ones. All members are publicly accessible and have default values.

.. autoclass:: PyTIQC.core.simtools.parameters
   :members:

decoherence 
-----------

This class contains the time vectors that parametrizes each error source, and the functions that generate them. It has a dictionary, accessible via ``decoherence.dict``, which contains the following keys. An individual error source is turned on or off by setting the boolean value of the value associated with that key::

  'none': True,
  'all': False, # turn everything on
  'addressing': False, # addressing error - constant during all evolution
  'dephase': False, # dephasing
  'heating': False, # heating
  'intensity': False, # intensity fluctuations
  'initerr': False, # qubit initialization error
  'specmode': False, # spectator mode coupling
  'spontdecay': False, # spont. emission from the qubit

.. autoclass:: PyTIQC.core.simtools.decoherence
   :members:

pulse and pulse sequences
-------------------------

The class ``pulse`` is a base class from which all other pulse types are inherited.

Most of the class member variables are parameters which gets passed to the function that constructs the Hamiltonian. Pulses are specified by their type, their rotation angles *theta* and *phi*, and addressed ion (default is global). The optional parameter *use_ideal* switches to using the ideal unitary matrix calculated by the functions in PyTIQC.ideal.gates for the state evolution.

.. autoclass:: PyTIQC.core.simtools.pulse
   :members:

.. autoclass:: PyTIQC.core.simtools.Rcar
   :members:

.. autoclass:: PyTIQC.core.simtools.Rblue
   :members:

.. autoclass:: PyTIQC.core.simtools.Rac
   :members:

.. autoclass:: PyTIQC.core.simtools.RMS
   :members:

.. autoclass:: PyTIQC.core.simtools.Delay
   :members:

.. autoclass:: PyTIQC.core.simtools.Hide
   :members:

.. autoclass:: PyTIQC.core.simtools.MeasInit
   :members:

A pulse sequence is a class which contains a Python list of pulses, plus some functions to operate on them.

.. autoclass:: PyTIQC.core.simtools.PulseSequence
   :members:

database
--------

The *database* class is used to store data from the simulation. Basic inputs are *T*, the time vector, and *Y*, the list of states at each timestep. It contains functions to do basic processing of the data. For example, if the *pulseseq* parameter is include, it calculates the states, populations, and density matrices at the end of each pulse.

.. autoclass:: PyTIQC.core.simtools.database
   :members:

other functions
---------------

.. autofunction:: PyTIQC.core.simtools.DecayedPopulations_CCD

.. autofunction:: PyTIQC.core.simtools.generate_decay_list

.. autofunction:: PyTIQC.core.simtools.DecayedExcitations_PMT
