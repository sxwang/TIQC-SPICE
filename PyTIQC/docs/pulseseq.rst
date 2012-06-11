===============
Pulse Sequences
===============

Each simulation requires a pulse sequence object, which consists of pulse objects.

Pulse
~~~~~

A pulse is an object defined with the following properties::

    theta # rotation angle
    phi # phase of the rotation
    ion # addressed ion (can be global)

Available pulse objects are: Rcar, Rblue, Rac, RMS. These corresponds to actual operations: carrier rotation, blue sideband rotation, AC stark shift pulse, MS gate. These pulses are instantiated by::

    sim.Rcar(params, theta, phi, ion)
    sim.Rblue(params, theta, phi, ion)
    sim.Rac(params, theta, phi, ion)
    sim.RMS(params, theta, phi)

For Rcar, Rblue, and RMS, ``ion`` is set to be ``-1`` by default, indicating a global rotation. 

The angle parameters ``theta`` and ``phi`` need to have ``pi`` included in their arguments. For example, a global pi carrier rotation about y (Ry pulse) is::

    sim.Rcar(params, pi, pi/2)

There are also pulse objects which represent special events which are not actual laser pulses:

Delay is the absense of a laser pulse, for example a delay time in the Ramsey experiment.This is specified as::

    sim.Delay(params, duration)

where ``duration`` is time in microseconds as usual.

Hide is an instruction to hide or unhide an ion. In the simulation, hiding works by setting the corresponding element in the addressing matrix to 0. For example, to hide ion 0::

    sim.Hide(params, ion=0, hide=True)

MeasInit is an instruction to perform a measurement and reinitialize on one ion, and store the result in a classical register::

    sim.MeasInit(params, ion=0)

Options
-------

All pulses have a ``use_ideal`` switch which selects whether the ideal unitary operation should be used (to compute the state evolution directly using the ideal rotation operator, and not evaluate the Hamiltonian), or do the full Hamiltonian evolution. By default this is set to False to use the full Hamiltonian evolution. To turn it on::

    pulse.use_ideal = True

All pulses have a function ``calculateIdealUnitary`` which uses the functions defined in the module ``PyTIQC.core.gates`` to compute the ideal unitary operator.

To examine the parameters of a pulse object, simply use ``print``::

    print pulse

The angle and phase will be displayed as numerical values (pi = 3.142). 

Pulse Sequence
~~~~~~~~~~~~~~

A pulse sequence is a class that maintains a list of pulses, and the functions that operaes on it. To define a pulse sequence, provide a list of pulses as an argument. For example, a pulse sequence consisting of two pulses is as follows::

    pulseseq = sim.PulseSequence( [
        sim.Rcar(params, pi, 0),
	sim.Rac(params, pi, 0, 0) ])

Pulse sequence definitions can be nested. A ``pulseseq`` object can be one element in the list that is used to instantiate another ``pulseseq``.

At the end of a pulseseq's initialization routine, the function ``makepulsesequence`` is called. This function goes through the list of pulses and increments a global time counter, and sets each pulse's ``starttime``, ``endtime``, and ``duration`` to be consistent with their order and rotation angles. This function and these variables avoid having to declare time as a global variable in the state evolution. However, note that ``makepulsesequence`` only sets these variables but does not check that the resulting times are consistent. In particular, if a pulse object appears twice in a pulse sequence (i.e. two pointers to the same pulse object are found at two different locations in the pulse sequence list), then the later occurrence overwrites the times set for the earlier pulse and the resulting state evolution will likely fail with hard-to-predict errors. The consistency of the times in a pulse sequence is, however, explicitly checked as each pulse is encountered in the actual simulation (in ``qc.simulateevolution``).

``pulseseq`` can be treated directly as a list, via the following mapping to built-in functions::

    len(pulseseq)       # returns the length of the list
    pulseseq[i]         # indexing operator is defined
    pulse in pulseseq   # the 'in' operator is defined
    print pulseseq      # display the list of pulses
    pulseseq + pulseseq1  # concatenate sequences

Custom functions are defined to modify the sequence itself::

    pulseseq.append(pulses)
    pulseseq.prepend(pulses)
    pulseseq.extend(pulseseq)

The argument to ``append`` and ``prepend`` are lists of pulses. The argument to ``extend`` is a PulseSequence object.

A pulse sequence consisting of global and addressed pulses can be modified to perform specific operations on different ions, simply by changing the ion indexing in a systematic way. For example, a sequence that performs a CNOT gate on ions 0 and 1 out of 3 ions can be changed to perform the gate on ions 0 and 2 simply by replacing all occurrences of ion=1 with ion=2. The function ``changeions`` is provided for this purpose. In the example case here::

    pulseseq.changeions(params, ionindices=[0,2,1])

The argument ``ionindices`` is a list that represents the permutation of the default sequence indices 0,1,2... etc.

Several functions are provided to visualize and summarize a pulse sequence.

* ``getdigitalrepresentation`` return a time vector with 1's where a laser is on and 0's where no laser is on. This is similar to what we see in the lab with an oscilloscope looking at the pulse signals sent to the AOMs.

* ``counttypes`` prints a summary of the numbers of pulse of each type in a pulse sequence.

* ``totalUnitary`` returns the ideal unitary that is implemented by the pulse sequence.

