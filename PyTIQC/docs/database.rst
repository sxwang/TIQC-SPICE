========
Database
========

The result of a simulation is stored in a ``database`` object. The object is instantiated as follows::

    data = sim.database(T, Y, hspace, pulseseq=None, register=[], statesvalid=True)

where

* ``T`` is a numpy array of times, from 0 to the end of the pulse sequence, where states are stored
* ``Y`` is the numpy array of states at the times in ``T``
* ``hspace`` is the hilbert space object for the simulation
* ``pulseseq`` is the pulse sequence that was used to run the simulation. This optional argument, if provided, allows selection of the values in Y that corresponds to states at the end of each pulse.
* ``register`` is the classical register used to store intermediate measurements in a pulse sequence (at each occurrence of the ``MeasInit`` pulse).
* ``statesvalid`` refers to whether the ``Y`` values are actually states or a vector of something else (e.g. populations instead). Setting ``statesvalid=False`` means that some functions in ``database`` that processes the data are not executed by default. 

The states ``Y`` are by default post-processed to provide populations and density matrices which can be accessed directly by the user. Some of the most useful ones are listed in the following table.

+-------------+-----------------------------------------------+
| Member      | Description                                   |
+=============+===============================================+
| T           | vector of times in the state evolution        |
+-------------+-----------------------------------------------+
| Y   	      | state vector vs time                          |
+-------------+-----------------------------------------------+
| YP   	      | population vector vs time                     |
+-------------+-----------------------------------------------+
| YtrN 	      | all state populations traced over phonons     |
+-------------+-----------------------------------------------+
| TP, YP, YRP | time, states, and populations at the end of   |
|      	      |	each pulse in a pulse sequence                |
+-------------+-----------------------------------------------+
| YRPN 	      | ``YRP`` traced over motional states           |
+-------------+-----------------------------------------------+
| RhoPN       | density matrices traced over motional states  |
|             | at the end of each pulse                      |
+-------------+-----------------------------------------------+
| RhoPNAll    | all density matrices for all MC instances. The|
|             | shape of this matrix is (number of MC         |
|             | instances, number of pulses,                  |
|             | ``hspace.dimensions``, ``hspace.dimensions``  |
+-------------+-----------------------------------------------+

``database`` objects have the functions ``add`` and ``mean`` defined, which are used internally by ``simulateevolution`` to combine multiple MC instance results. For example::

    data = data1 + data2

Here, the populations and density matrices are summed. Then::

    data = mean(data)

Here, the (summed) populations and density matrices are divided by the total number of runs to get the averaged result of a group of MC instances.

Visualization
~~~~~~~~~~~~~

The variables in a ``database`` can be visualized in several ways, many of which maps to actual plots one might see in the experiment. Here, the numerical argument is the figure number as displayed by matplotlib. If the argument is set to 0, no figure is plotted. 

* ``statesatpulse(1)`` plots the population of each state at the end of each pulse in a pulse sequence, as a function of pulse number. For an N-ion system, there would be 2^N lines plotted, each line containing a number of points equal to the number of pulses + 1 (the initial state is also plotted). Phonons are traced out.
* ``endpopulation`` prints the population of each state (including phonons) at the end of the sequence.
* ``tracedpopulation(1)`` plots the D state population of each ion as a function of time. Phonons are traced out. This is similar to the data we get from the experiment from the CCD camera (assuming the CCD has single-ion resolution).
* ``displaypopulation(1)`` is similar to ``tracedpopulation`` except that phonons are not traced out.
* ``displayPMTpopulations(1)`` plots the data we get from the experiment as measured by the PMT. The PMT can only measure the total number of bright ions, but not which ions are bright. For an N-ion system, there would be N+1 lines, each line corresponding to there being 0..N bright ions.
* ``displayPhaseSpace(1)`` plots the state evolution in phase space (x, p). Useful for investigating the MS gate.

Saving data
~~~~~~~~~~~

The ``save`` function saves a ``database`` object to a pickle file ::

    data.save(filename)

If filename is omitted, by default it is a timestamp appended with ``-data.pk``.

The entire simulation including the pulseseq, params, and dec objects can be saved using the following function in the ``qc`` module::

    qc.saveRun(pulseseq, params, dec, data, filename, clear=True)

The saved file is a Python shelve object with 4 keys, which are the pulseseq, params, dec, and data objects. Beware that these files can easily get very large. The argument ``clear`` removes the largest database components ``T``, ``Y``, and ``YR`` to save some space. 
The shelve objects can be loaded with::

    pulseseq, params, dec, data = qc.loadfile(filename).
