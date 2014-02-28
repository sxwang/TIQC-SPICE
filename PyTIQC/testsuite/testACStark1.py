''' test script for AC Stark pulse.
do a single AC stark pulse and show the state evolution. '''

import numpy as np
import PyTIQC.core.simtools as sim
import PyTIQC.core.qctools as qc
import matplotlib.pyplot as pl
import scipy.optimize as spop
pi = np.pi

hspace = sim.hspace(1,2,0,0)
params = sim.parameters(hspace)
dec = sim.decoherence(params)

params.y0[qc.stateToIndex('S'+',0', hspace)] = 1
params.stepsize = 1
params.printpulse = False
params.progbar = False

# change this to check amplitude of off-resonant excitations
#params.ACintenFac = 60
#params.recalculate()

#params.solver = "Cheby" # seems to not work for this one

pulseseq = sim.PulseSequence( [ \
  sim.Rcar(params, pi/2, 0, 0), \
  sim.Rac(params,pi, 0, 0), \
  sim.Rcar(params, pi/2, 0, 0) \
  ] )

# this line turns on the ODE solver for the AC pulse
# then the total population exceeds 1, growing linearly w/ time
#pulseseq[1].dotimedepPulse = True
#pulseseq[1].duration = 40

data = qc.simulateevolution(pulseseq, params, dec)

data.tracedpopulation(1)

#plot the total population for all states vs time
#pl.figure(2)
#pl.plot(data.T, np.sum(data.YR,axis=1))
#pl.xlabel('time [us]')
#pl.ylabel('total population')
#pl.show()


def data_optimize(inputparam, params, dec):
    ''' use this function to search for params.ACcorr given params.omac.  '''
    params.ACcorr = inputparam
    data = qc.simulateevolution(pulseseq, params, dec)
    data.tracedpopulation(1)

    return data.YtrI[-1]
    
print spop.fmin(data_optimize, 0.09, args=(params, dec), full_output=True)
print "for params.omac=", params.omac
