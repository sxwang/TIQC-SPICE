''' test of MS gates '''

import numpy as np
import PyTIQC.core.simtools as sim
import PyTIQC.core.qctools as qc
import time
import matplotlib.pyplot as pl
import PyTIQC.tools.progressbar as progbar
import scipy.optimize as spop

pi = np.pi

NumberOfIons = 2
PhononOverhead = 7
hspace = sim.hspace(NumberOfIons,2,NumberOfIons+PhononOverhead,0)
params = sim.parameters(hspace)
dec = sim.decoherence(params)

#dec.doRandNtimes = 1
#dec.dict['all'] = True
#dec.doPP = False
#params.use_servers( ['local'])
params.progbar = True

#params.shortestMS = 4
params.calcPulseParams()

pulseseq = sim.PulseSequence( [ \
#    sim.Hide(params, 3, True),
#    sim.Hide(params, 4, True), 
    sim.RMS(params, pi/2, 0),
    #sim.Rcar(params, pi/2, pi/2)
    ] )

data = qc.simulateevolution(pulseseq, params, dec)

data.displaypopulation(0)
data.displayPMTpopulations(1)
#data.displayPhaseSpace(0)

####################################################
# extra post-run diagnostic functions


def data_optimize(inputparam, params, dec):
    ''' use this function to search for params.MScorr for every detuning. 
    How-to:
    - in the root-level code above, enter params.MSdelta and params.shortestMS
    - here put 'inputparam' in place of the number you want to optimize
    - run the file and enter the result in above
    - repeat for every number that appears in params.MScorr '''
    #params.MScorr = {
    #    'duration': 1.1642,
    #    'omrabi': 1.0268,
    #    'shortPulseFactorPos': {2: 1.000, 4: 1.0106, 8: inputparam, 16: 1.01835},
    #    'shortPulseFactorNeg': {2: 1.010, 4: 1.0214, 8: 1.0271, 16: 1.02986}
    #    }
    pulseseq = sim.PulseSequence( [
        sim.Hide(params, 0, True),
        sim.Hide(params, 1, True),
        sim.RMS(params, pi/2,0),
        sim.Rcar(params, pi/2, pi/2)] )
    pulseseq[2].omrabi_r = params.omc_ms * inputparam
    pulseseq[2].omrabi_b = params.omc_ms * inputparam

    data = qc.simulateevolution(pulseseq, params, dec)
    data.displayPMTpopulations(0)

    inbalance = abs(data.P_PMT_end[2]-data.P_PMT_end[5])
    unwanted = abs(data.P_PMT_end[0]+data.P_PMT_end[1]+data.P_PMT_end[3]+data.P_PMT_end[4])
    
    print inputparam, inbalance+unwanted

    return inbalance+unwanted
    
#print spop.fmin(data_optimize, 1., args=(params, dec), full_output=True)
#print "omrabi factor =", 
#print "for params.MSdelta=", params.MSdelta, "shortestMS=", params.shortestMS


def parity(data, params):
    ''' check the parity of the pi/2 MS gate.
    must run expMS first. call: parity(data.Yend) '''

    phase = np.linspace(0,2*pi,20)
    params.y0 = data.Yend
    P_PMT = data.P_PMT
    parity = np.zeros(len(phase))
    params.printpulse = False
    dec = sim.decoherence(params)

    widgets = [progbar.Percentage(), ' ', progbar.Bar(),' ', progbar.ETA()]
    pbar = progbar.ProgressBar(widgets=widgets).start()

    for i in range(len(phase)):
        pulseseq = sim.PulseSequence( [ sim.Rcar(params, pi/2, phase[i], -1) ])
        data = qc.simulateevolution(pulseseq, params, dec)
        data.tracedpopulation(0)
        parity[i] = data.YtrN[-1,0] + data.YtrN[-1,3] - data.YtrN[-1,1] - data.YtrN[-1,2]
        pbar.update(int(1.*(i+1)*100/(len(phase))))

    def sinefunc(x, a, b):
        return -a*np.sin(2*x)+b

    [p,pcov] = spop.curve_fit(sinefunc, phase, parity)

    population = P_PMT[-1,0] + P_PMT[-1,-1]
    coherence = abs(p[0])
    fidelity = (population+coherence)/2

    xphase = np.linspace(0, 2*pi, 200)
    parityfit = sinefunc(xphase, p[0], p[1])
    pl.plot(phase, parity,'.', xphase, parityfit, '-')
    pl.xlabel('phase [rad]')
    pl.ylabel('parity')
    pl.show()

    print "population = ", population
    print "coherence contrast = ", coherence
    print "fidelity = ", fidelity
    return phase, parity, pulseseq

def MSpulseshaped():
    ''' demonstration of pulse shaping removing fast oscillations '''
    tic = time.time()

    params.shape = 5
    pulselen = np.linspace(10./86*pi/2, 2*pi, 4*160)
    realpulselen = np.zeros(len(pulselen))
    result = np.zeros([len(pulselen), 2])

    for i in range(len(pulselen)):
        print "pulse length ", pulselen[i]
        pulseseq = sim.PulseSequence( [ sim.RMS(params, pulselen[i], 0, -1) ] )
        realpulselen[i] = pulseseq[0].duration
        data = qc.simulateevolution(pulseseq, params, dec)
        data.tracedpopulation(0)
        result[i,:] = 1-data.YtrI[-1]

    toc = time.time()
 
    pl.plot(realpulselen, result)
    pl.show()

    print "MSpulseshaped() runtime ", toc-tic, "sec"

    return realpulselen, result


def timestepcheck(params, dec):
    ''' check that timestep chosen for ODE is converged '''
    timesteps = np.logspace(0, -3, 10)
    result = np.zeros([len(timesteps), 2])

    params.shape = 5
    params.printpulse = False
    pulseseq = sim.PulseSequence( [ sim.RMS(params, pi/2, 0, -1) ] )

    for i in range(len(timesteps)):
        params.ODEtimestep = timesteps[i]
        print "timestep = ", params.ODEtimestep
        data = qc.simulateevolution(pulseseq, params, dec)
        result[i,0] = abs(data.Yend[0])**2
        result[i,1] = abs(data.Yend[int(3*len(data.Yend)/4)])**2

    pl.semilogx(timesteps, result[:,0], '.:')
    pl.xlabel('ODE timestep, [us]')
    pl.ylabel('State population')
    pl.title('convergence of ODE solver vs timestep')
#    pl.legend(['p(DD,0)', 'p(SS,0)'], loc=2)

    return timesteps, result

