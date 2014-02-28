''' Cirac-Zoller CNOT '''

import numpy as np
import PyTIQC.core.simtools as sim
import PyTIQC.core.qctools as qc

pi = np.pi

hspace = sim.hspace(2,2,2,0)
params = sim.parameters(hspace)
dec = sim.decoherence(params)

params.y0[qc.stateToIndex('SS,0', hspace)] = 1
params.stepsize = 10
params.ignorelightshift = 1
#params.addressing = np.array([[1,0.1],[0.1,1]])

ang = np.arccos(-np.real(np.exp(pi/2*1j*np.sqrt(2))));
pulseseq = sim.PulseSequence( [ \
    sim.Rblue(params, pi, 0, 0), \
    sim.Rcar(params, pi/2, 0, 1), \
    sim.Rblue(params, pi/2, 0, 1), \
    sim.Rblue(params, np.sqrt(2)*pi, pi/2, 1), \
    sim.Rblue(params, pi/2, pi, 1), \
    sim.Rcar(params, pi/2, pi + ang, 1), \
    sim.Rblue(params, pi, 0, 0) \
    ] )

data = qc.simulateevolution(pulseseq, params, dec)

data.displaypopulation(1)
data.tracedpopulation(2)


