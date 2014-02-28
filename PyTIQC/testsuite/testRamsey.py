''' Ramsey sequence with/without decoherence '''

import numpy as np
import PyTIQC.core.simtools as sim
import PyTIQC.core.qctools as qc
import PyTIQC.tools.progressbar as progbar
from scipy.optimize import leastsq
import matplotlib.pyplot as pl

doRun = True

pi = np.pi

NumberOfIons = 1
NumberofPhonons = 2
hspace = sim.hspace(NumberOfIons,2,NumberofPhonons,0)
params = sim.parameters(hspace)
params.coherenceTime = 5000
params.correlationTime = 5200
dec = sim.decoherence(params)

dec.doRandNtimes = 8
dec.dict['dephase'] = True
dec.progbar = False

params.stepsize = 1
params.ODEtimestep = 1
params.detuning = 2*pi*0.002
params.progbar = False
dec.doPP = True
params.use_servers( ['local'])
params.pplog = False
dec.doPPprintstats = False

# calcpulseparams here if coherencetime is changed

def doRamseyRun(hspace, params, dec):
    tdelay = np.linspace(0,5000,500)
    YR = np.zeros([len(tdelay),hspace.dimensions], np.complex64)
    pop = np.zeros([len(tdelay),hspace.dimensions], np.complex64)

    widgets = [progbar.Percentage(), ' ', progbar.Bar(),' ', progbar.ETA()]
    pbar = progbar.ProgressBar(widgets=widgets).start()

    for i in range(len(tdelay)):
        pulseseq = sim.PulseSequence( [ \
            sim.Rcar(params, pi/2, 0), \
            sim.Delay(params, tdelay[i]), \
            sim.Rcar(params, pi/2, pi/2) \
            ] )

        data = qc.simulateevolution(pulseseq, params, dec)

        YR[i,:] = data.YR[-1,:]

        pbar.update(int(1.*(i+1)*100/(len(tdelay))))

    data1 = sim.database(tdelay, YR, hspace, pulseseq=None, statesvalid = False)
    data1.tracedpopulation(0)

    [fitfun, par] = doRamseyFit(data1)

    return data1, fitfun, par

def doRamseyFit(data1):
    p0 = [0.5, params.detuning, pi/2, 0.5, 2000]
    fitfun = lambda p, x: p[0] * np.cos(p[1]*x+p[2]) * np.exp(-np.log(2)*(x/p[4])**2) + p[3]
    errfun = lambda p, x, y: y-fitfun(p,x)
    x = data1.T
    y = np.sum(data1.YtrI, axis=1) / data1.hspace.nuions
    par, covmat, infodict, mesg, ier = leastsq(errfun, p0, args=(x,y), full_output = True)
    return fitfun, par

##################
if doRun:
    data1, fitfun, par = doRamseyRun(hspace, params, dec)
    print "coherence time (us) = ", par[-1]
    pl.plot(data1.T, data1.YtrI)
    pl.plot(data1.T, fitfun(par, data1.T))
    pl.show()
