=================
Simulateevolution
=================

After defining the ``hspace``, ``parameters``, ``decoherence``, and ``pulsesequence`` objects, we can run the simulation::

    data = qc.simulateevolution(pulseseq, params, dec)

The printed output of this function during its execution depends on the state of the variables ``params.progbar``, ``params.printpulse``, and ``dec.progbar``. Usually it will be a progressbar displaying the progress and estimated time remaining.

This function returns a ``database`` object which contains the output of the simulation, processed multiple ways to map to various metrics we might use in the experiment and analysis, and includes functions to further process and visualize the data. This object is described in the next page of this tutorial.
