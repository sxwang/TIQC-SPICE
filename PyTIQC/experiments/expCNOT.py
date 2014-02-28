''' a 2-ion CNOT in a 3-ion chain based on pulse sequence from Nebendahl (2009) '''

import numpy as np
import PyTIQC.core.simtools as sim
import PyTIQC.core.qctools as qc
import matplotlib.pyplot as pl

import time

pi = np.pi

NumberOfIons = 3
NumberOfPhonons = 7
hspace = sim.hspace(NumberOfIons,2,NumberOfPhonons,0)
params = sim.parameters(hspace)

params.stepsize = 100
params.printpulse = True

dec = sim.decoherence(params)

pulseseq = sim.PulseSequence( [ \
    sim.Rcar(params, pi/2, 0), 
    sim.Rac(params, pi/2, 0, 0), 
    sim.RMS(params, pi/4, 0), 
    sim.Rcar(params, pi/4, 0), 
    sim.Rac(params, pi, 0, 2), 
    sim.Rcar(params, pi/4, 0), 
    sim.RMS(params, pi/4, 0), 
    sim.Rac(params, pi/2, 0, 0), 
    sim.Rcar(params, pi/2, 0), 
    sim.Rac(params, pi, 0, 2)
    ] )

#params.y0[qc.stateToIndex('DDD' + ',0', hspace)] =  0.25-0.25j
#params.y0[qc.stateToIndex('DDS' + ',0', hspace)] = -0.25-0.25j
#params.y0[qc.stateToIndex('DSD' + ',0', hspace)] = -0.25-0.25j
#params.y0[qc.stateToIndex('DSS' + ',0', hspace)] = -0.25+0.25j
#params.y0[qc.stateToIndex('SDD' + ',0', hspace)] =  0.25-0.25j
#params.y0[qc.stateToIndex('SDS' + ',0', hspace)] = -0.25-0.25j
#params.y0[qc.stateToIndex('SSD' + ',0', hspace)] = -0.25-0.25j
#params.y0[qc.stateToIndex('SSS' + ',0', hspace)] = -0.25+0.25j

tic = time.clock()
data = qc.simulateevolution(pulseseq, params, dec)
toc = time.clock()
print "runtime: ", toc-tic, " sec"

data.statesatpulse(2)
data.endpopulation()
