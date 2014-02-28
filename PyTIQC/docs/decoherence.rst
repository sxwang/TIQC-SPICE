===========
Decoherence
===========

Various error sources in the experiment are simulated using a quantum Monte-Carlo (MC) approach. A random trajectory of a particular parameter is generated, and a single simulation realizes an individual evolution with a single trajectory. The ensemble average is obtained by repeating the simulation with different random trajectories and taking the average of the populations and density matrices. Errors are represented in the simulation as modifications to the Hamiltonian or to global parameters.

The decoherence class is instantiated with a parameter object as an argument::

    dec = sim.decoherence(params)

A decoherence object must be created even if no errors are to be turned on for any particular simulation. By default, the simulation of decoherence effects are turned off. To turn it on and to control the number of MC instances, use the following::

    dec.doRandNtimes = 0

Setting ``doRandNtimes`` to a non-zero value turns on the noise objects, but unless the noise sources are turned on, all the Hamiltonian modifications are set to 0. 

By default, all noise sources are off. To turn on specific noise sources or all noise sources, set the relevant variable in the noise dictionary::

    dec.dict['all'] = True # turn on all noise sources

The available keys, each mapping to a noise source, are::

    addressing, dephase, heating, hiding, intensity, initerr, specmode, spontdecay, switchtime

To run multiple MC instances, much can be gained in execution time by enabling parallel processing via ``Parallel Python``. By default, PP is disabled. To enable it::

    dec.doPP = True
    dec.doPPprintstats = True

``doPPprintstats`` controls whether to print server statistics after a block of MC instances have been executed. It is enabled by default. 

To show a progressbar for the MC instances (enabled by default)::

    dec.progbar = True

The progressbar increments only when an instance completes and is received by the host. The estimated time remaining may not be representative (unless there is a very large number of instances so the average would be approximately correct), and should be ignored. 
