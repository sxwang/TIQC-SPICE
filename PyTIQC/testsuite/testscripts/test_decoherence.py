
#!/usr/bin/env python
# -*- mode: Python; coding: latin-1 -*-
# Time-stamp: "2011-05-30 23:02:32 tomi"

#  file       test_evaluations.py
#  author     Thomas Monz 2011

import unittest
import PyTIQC.core.qctools as qc
import PyTIQC.core.simtools as sim
import numpy as np

import PyTIQC.evaluation.InvestigatePulseSeq as IPS
import PyTIQC.evaluation.densitymatrixreconstruction as dmr
import PyTIQC.evaluation.processtomography.proctom as proctom

pi = np.pi

from scipy.optimize import leastsq
import PyTIQC.tools.progressbar as progbar

class TestUserFunction(unittest.TestCase):

    plot = False
    previous_tests = 0
    figNum = 0

    def getFigNum(self):
	figNum = 0
        if(self.plot):
            TestUserFunction.figNum += 1
	    figNum = TestUserFunction.figNum + TestUserFunction.previous_tests
            print "(figure",figNum,")"
        return figNum

    def test_state_initialisation(self):
        """
        do Rabi with initerror and check for correct
        reduction in contrast
        """
        assert(test_state_initialisation_detailed(self.getFigNum()))
        #assert(True)

    def test_dephasing(self):
        """
        do Ramsey and look for decay
        """
        assert(test_dephasing_detailed(self.getFigNum()))

    def test_spontdecay(self):
        """
        wait for a Delay pulse and watch the D state decay to S
        """
        assert(test_spontdecay_detailed(self.getFigNum()))

    def test_heating(self):
        """
        wait for a Delay pulse and watch the phonon number increase
        """
        assert(test_heating_detailed(self.getFigNum()))


# Test class with plotting enabled

class TestUserFunctionPlots(TestUserFunction):
        plot = True


# Collect all test suites for running

all_suites = unittest.TestSuite((
  unittest.makeSuite(TestUserFunction)
  ))

plot_suites = unittest.TestSuite((
  unittest.makeSuite(TestUserFunctionPlots)
  ))


def run():
    unittest.TextTestRunner(verbosity=2).run(all_suites)


# ---------------------------

###########################################
# test_state_initialisation_detailed
###########################################
def test_state_initialisation_detailed(figNum):
    NumberOfIons = 1
    NumberOfPhonons = 1

    hspace = sim.hspace(NumberOfIons,2,NumberOfPhonons,0)
    params = sim.parameters(hspace)
    params.stateiniterr = 0.2

    dec = sim.decoherence(params)

    dec.doRandNtimes = 1000
    dec.dict['initerr'] = True

    params.y0[qc.stateToIndex(NumberOfIons*'S'+',0', hspace)] = 1

    pulseseq = sim.PulseSequence( [ \
        sim.Rcar(params, 5*pi, 0),
        ] )

    data = qc.simulateevolution(pulseseq, params, dec)
    data.tracedpopulation(figNum)

    # start with fit here
    x = data.T
    y = data.YtrN.transpose()[-1] # this is the d-state population

    # p[0] ... amplitude, should be 0.8
    # p[1] ... freq, should be params.omc
    # p[2] ... phase, should be 0
    # p[3] ... offset, should be 0

    startparams = np.array([0.8, params.omc, 0, 0])

    # 1-data ... to get the D-state population
    fitfun = lambda p, x: p[0]*np.sin(p[1]*x/2+p[2])**2 + p[3]
    errfun = lambda p, x, y: y-fitfun(p,x)

    par, covmat, infodict, mesg, ier = leastsq(errfun, startparams, args=(x,y), full_output = True)
#    print par

    # even for the 1000 realisations, allow for a 3% error
    epsilon = 0.03
    return abs(abs(par[0])-startparams[0]) < epsilon

###########################################
# test_dephasing_detailed
###########################################
def test_dephasing_detailed(figNum):
    NumberOfIons = 1
    NumberofPhonons = 0
    hspace = sim.hspace(NumberOfIons,2,NumberofPhonons,0)
    params = sim.parameters(hspace)
    dec = sim.decoherence(params)

    params.coherenceTime = 1000
    params.correlationTime = 0.1
    dec.doRandNtimes = 20
    dec.dict['dephase'] = True
    dec.progbar = False

    params.stepsize = 10
    params.ODEtimestep = 1
    params.detuning = 2*pi*0.01  #kHz
    params.progbar = False
    #dec.doPP = True
    #params.use_servers(['local'])
    #dec.doPPprintstats = False

    tdelay = np.linspace(0,1000,100)
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
    data1.tracedpopulation(figNum)

    # fitting part
    p0 = [0.5, params.detuning, pi/2, 0.5, params.coherenceTime]
    fitfun = lambda p, x: p[0] * np.cos(p[1]*x+p[2]) * np.exp(-np.log(2)*(x/p[4])**2) + p[3]
    errfun = lambda p, x, y: y-fitfun(p,x)
    x = data1.T
    y = data1.YR[:,0]
    par, covmat, infodict, mesg, ier = leastsq(errfun, p0, args=(x,y), full_output = True)
    
    epsilon = 100 # with 20 iterations allow 100us offset in coherence time
    #print par
    return data1
    #return np.abs(par[-1] - params.coherenceTime) < epsilon


###########################################
# test_spontdecay_detailed
###########################################
def test_spontdecay_detailed(figNum):
    NumberOfIons = 1
    NumberofPhonons = 1
    hspace = sim.hspace(NumberOfIons,2,NumberofPhonons,0)
    params = sim.parameters(hspace)
    dec = sim.decoherence(params)

    dec.doRandNtimes = 100
    dec.dict['spontdecay'] = True
    dec.doPP = True
    dec.doPPprintstats = False
    dec.progbar = False
    # for the test we set it to 300 mus, instead of 1.168 s
    params.lifetime = 300

    params.y0[qc.stateToIndex(NumberOfIons*'D'+',0', hspace)] = 1
    params.y0[qc.stateToIndex(NumberOfIons*'S'+',0', hspace)] = 0

    pulseseq = sim.PulseSequence( [ \
        sim.Delay(params, 1000), \
        ] )

    data = qc.simulateevolution(pulseseq, params, dec)
    data.tracedpopulation(figNum)

    # fitting part
    p0 = [1, params.lifetime]
    fitfun = lambda p, x: p[0] * np.exp(-x / float(p[1]))
    errfun = lambda p, x, y: y-fitfun(p,x)
    x = data.T
    y = data.YtrN[:,0]
    par, covmat, infodict, mesg, ier = leastsq(errfun, p0, args=(x,y), full_output = True)
    
    epsilon = 50 # with 100 iterations allow 50us offset in decay time

#    print np.abs(par[-1] - params.lifetime)
    return np.abs(par[-1] - params.lifetime) < epsilon

###########################################
# test_heating_detailed
###########################################
def test_heating_detailed(figNum):
    NumberOfIons = 1
    NumberofPhonons = 10
    hspace = sim.hspace(NumberOfIons,2,NumberofPhonons,0)
    params = sim.parameters(hspace)
    dec = sim.decoherence(params)

    # since data.RhoPNAll doesn't save phonons, we'll have to simulate sequentially
    numRuns = 5
    phonons = np.zeros(numRuns)

    dec.doRandNtimes = 1
    dec.dict['heating'] = True
    dec.progbar = False
    params.progbar = False
    params.heatingrate = 250
    maxDelay = 1000

    pulseseq = sim.PulseSequence( [ \
        sim.Delay(params, maxDelay), \
        ] )
    widgets = [progbar.Percentage(), ' ', progbar.Bar(),' ', progbar.ETA()]
    pbar = progbar.ProgressBar(widgets=widgets).start()

    for i in range(numRuns): 
        data = qc.simulateevolution(pulseseq, params, dec)
        phonons[i] = qc.indexToState(np.nonzero(data.Yend)[0][0], hspace)[0]
        pbar.update(int(1.*(i+1)*100/numRuns))
    
    epsilon = np.abs(maxDelay/params.heatingrate - np.mean(phonons))

    return epsilon < 2*np.std(phonons)  # allow 2 sigma deviation
    #return data


###################################################
if __name__ == "__main__":
    #print test_state_initialisation_detailed()
    #data1= test_dephasing_detailed(0)
    #print test_spontdecay_detailed(1)
    data = test_heating_detailed(0)

# test_evaluations.py ends here
