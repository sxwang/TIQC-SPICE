""" The simplest possible script, showing a simulation of Rabi flops. """

import numpy as np
import PyTIQC.core.simtools as sim
import PyTIQC.core.qctools as qc

pi = np.pi

NumberOfIons = 1
NumberofPhonons = 10

hspace = sim.hspace(NumberOfIons,2,NumberofPhonons,0)
hspace.initial_state('thermal', nBar=5)  # use thermal state
params = sim.parameters(hspace)
params.eta = np.array([0.3,0.3])
dec = sim.decoherence(params)

params.use_servers([ 'anna' ])

dec.doRandNtimes = 8
dec.dict['all'] = True
dec.doPP = True

dec.progbar = False
params.progbar= True
params.printpulse = True
params.LDapproximation = False
params.addressingerr = 0.03

pulseseq = sim.PulseSequence( [ \
    sim.Rcar(params, 10*pi, 0, 0), \
    ] )
#pulseseq[0].dotimedepPulse = True

data = qc.simulateevolution(pulseseq, params, dec)

data.tracedpopulation(2)
