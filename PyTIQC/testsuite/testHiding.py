""" Test hiding/unhiding and associated errors """

import numpy as np
import PyTIQC.core.simtools as sim
import PyTIQC.core.qctools as qc

pi = np.pi

NumberOfIons = 4
NumberofPhonons = 0

hspace = sim.hspace(NumberOfIons,2,NumberofPhonons,0)
params = sim.parameters(hspace)
dec = sim.decoherence(params)

params.hidingerr = 0.75

dec.doRandNtimes = 4
dec.dict['hiding'] = True
#dec.doPP = True

dec.progbar = False
#params.progbar= True
#params.printpulse = True

pulseseq = sim.PulseSequence( [ \
    sim.Hide(params, 0, True),
    sim.Hide(params, 0, False),
    sim.Rcar(params, 2*pi, 0), \
    sim.Hide(params, 2, True),
    sim.Hide(params, 2, False),
    sim.Rcar(params, 2*pi, 0), \
    sim.Hide(params, 2, True),
    sim.Hide(params, 2, False),
    sim.Rcar(params, 2*pi, 0), \
    ] )

data = qc.simulateevolution(pulseseq, params, dec)

data.tracedpopulation(1)
