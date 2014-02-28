#!/usr/bin/env python
# -*- mode: Python; coding: latin-1 -*-
# Time-stamp: "2011-05-27 15:17:46 c704252"

#  file       test_carrier.py
#  author     Thomas Monz 2011

import unittest
import PyTIQC.core.qctools as qc
import PyTIQC.core.simtools as sim
import numpy as np

import pylab as pl

from scipy.optimize import leastsq
import PyTIQC.tools.progressbar as progbar

pi = np.pi


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

    def test_Rabi_carrier(self):
        """
        add functionality to do carrier rabi flops
        if the frequency is correct _and_ we have
        full contrast, let it pass
        """
        assert(test_Rabi_carrier_detailed(self.getFigNum()))

    def test_Rabi_detuned(self):
        """
        add functionality to do detuned carrier rabi flops
        if the frequency is correct (sqrt(omega^2 + delta^2))
        _and_ we have the correct amplitude (omega^2/(omega^2+delta^2))
        let it pass
        """
        assert(True)

    def test_Ramsey_carrier(self):
        """
        check that (on the carrier):
          Rcar(0.5,0,1)
          Rcar(0.5,phi,1)

        gives expected oscillation
        """
        assert(test_Ramsey_carrier_detailed(self.getFigNum()))

    def test_Ramsey_detuned(self):
        """
        check that a detuning of delta:
          Rcar(0.5,0,1, carrier)
          waittime
          Rcar(0.5,0.5,1, detuned)
        gives oscillations with the detuning
        frequency
        """
        assert(True)

    def test_ACStark(self):
        """
        with default detuning (-80MHz) check that 
          Rcar(0.5, 0)
          Rac(time, 0, 0)
          Rcar(0.5, 0)
        gives oscillations with the expected frequency.
        """
        assert(test_ACStark_detailed(self.getFigNum()))

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

def test_Rabi_carrier_detailed(figNum):
    #print TestUserFUnction.figNum, TestUserFunction.figNum_start, "<<<<"
    NumberOfIons = 1
    PhononOverhead = 2

    hspace = sim.hspace(NumberOfIons,2,NumberOfIons+PhononOverhead,0)
    params = sim.parameters(hspace)
    params.stepsize = 1
    dec = sim.decoherence(params)

    params.y0[qc.stateToIndex('S,0', hspace)] = 1
    params.printpulse = False # don't print pulse details

    pulseseq = sim.PulseSequence( [
        sim.Rcar(params, 10*pi, 0, -1)
        ] )

    data = qc.simulateevolution(pulseseq, params, dec)

    data.tracedpopulation(figNum)

    # start with fit here
    x = data.T
    y = data.YtrN.transpose()[0] # this is the s-state population

    # p[0] ... amplitude, should be 1
    # p[1] ... freq, should be params.omc
    # p[2] ... phase, should be 0
    # p[3] ... offset, should be 0

    startparams = np.array([1, params.omc, 0, 0])

    # 1-data ... to get the D-state population
    fitfun = lambda p, x: 1-p[0]*np.sin(p[1]*x/2+p[2])**2 + p[3]
    errfun = lambda p, x, y: y-fitfun(p,x)

    par, covmat, infodict, mesg, ier = leastsq(errfun, startparams, args=(x,y), full_output = True)

    #print startparams
    #print par

    #print startparams-par

    epsilon = 10**-5
    if par[0] - startparams[0] > epsilon:
        print "amplitude of Rabi oscillations wrong"
    if par[1] - startparams[1] > epsilon:
        print "frequency of Rabi oscillations wrong"
    if par[2] - startparams[2] > epsilon:
        print "phase of Rabi oscillations wrong"
    if par[3] - startparams[3] > epsilon:
        print "offset of Rabi oscillations wrong"


    return np.all(par-startparams < epsilon)

def test_Ramsey_carrier_detailed(figNum):
    NumberOfIons = 1
    PhononOverhead = 2

    hspace = sim.hspace(NumberOfIons,2,NumberOfIons+PhononOverhead,0)
    params = sim.parameters(hspace)
    dec = sim.decoherence(params)

    params.y0[qc.stateToIndex('S,0', hspace)] = 1
    params.stepsize = 1
    params.printpulse = False # don't print pulse details

    numberofpoints = 20
    phi = np.linspace(0, 2*pi, numberofpoints)
    ex = np.zeros(numberofpoints)

    for i in range(numberofpoints):
        pulseseq = sim.PulseSequence( [
            sim.Rcar(params, pi/2, 0, -1),
#            sim.Delay(params, tdelay[i]),
            sim.Rcar(params, pi/2, phi[i], -1)
            ])

        data = qc.simulateevolution(pulseseq, params, dec)
        data.tracedpopulation(figNum)
        ex[i] = data.YtrN[-1,0]

    # fig1 = pl.figure(1)
    # ax1 = fig1.add_subplot(111)
    # ax1.plot(phi, ex)
    # fig1.show()


    # p[0] ... amplitude, should be 1
    # p[1] ... because phase is in units of pi -> 1
    # p[2] ... phase, should be 0
    # p[3] ... offset, should be 0.5

    startparams = np.array([1, 1, 0, 0.5])

    # 1-data ... to get the D-state population
    fitfun = lambda p, x: p[0]/2*np.cos(p[1]*x+p[2]) + p[3]
    errfun = lambda p, x, y: y-fitfun(p,x)

    par, covmat, infodict, mesg, ier = leastsq(errfun, startparams, args=(phi,ex), full_output = True)

    #print startparams
    #print par

    #print startparams-par

    epsilon = 10**-5

    if par[0] - startparams[0] > epsilon:
        print "amplitude of Ramsey experiment wrong"
    if par[1] - startparams[1] > epsilon:
        print "frequency of Ramsey experiment wrong"
    if par[2] - startparams[2] > epsilon:
        print "phase of Ramsey experiment wrong"
    if par[3] - startparams[3] > epsilon:
        print "offset of Ramsey experiment wrong"


    return np.all(par-startparams < epsilon)


    # estimates on the fit parameters

def test_ACStark_detailed(figNum):
    NumberOfIons = 1
    hspace = sim.hspace(NumberOfIons,2,0,0)
    params = sim.parameters(hspace)
    dec = sim.decoherence(params)

    params.progbar = False

    pulseseq = sim.PulseSequence( [ \
        sim.Rcar(params, pi/2, 0),
        sim.Rac(params, 0, 0, 0),
        sim.Rcar(params, pi/2, 0)
        ] )

    ACtime = np.linspace(0,10*pi,100)
    realtime = np.zeros_like(ACtime)

    YR = np.zeros([len(ACtime), hspace.dimensions], dtype='complex')

    for i in range(len(ACtime)):
        pulseseq[1] = sim.Rac(params, ACtime[i], 0, 0)
        realtime[i] = pulseseq[1].duration
        data = qc.simulateevolution(pulseseq, params, dec)
        YR[i,:] = data.YR[-1,:]


    data1 = sim.database(realtime, YR, hspace, statesvalid = False)
    data1.tracedpopulation(figNum)

    # start with fit here
    x = data1.T
    y = data1.YtrN.transpose()[0] # this is the s-state population

    # p[0] ... amplitude, should be 1
    # p[1] ... freq, should be params.omc
    # p[2] ... phase, should be 0
    # p[3] ... offset, should be 0

    startparams = np.array([1, params.omac, 0, 0])

    # 1-data ... to get the D-state population
    fitfun = lambda p, x: 1-p[0]*np.sin(p[1]*x/2+p[2])**2 + p[3]
    errfun = lambda p, x, y: y-fitfun(p,x)

    par, covmat, infodict, mesg, ier = leastsq(errfun, startparams, args=(x,y), full_output = True)

    epsilon = 10**-3
    if par[0] - startparams[0] > epsilon:
        print "amplitude of AC oscillations wrong"
    if par[1] - startparams[1] > epsilon:
        print "frequency of AC oscillations wrong"
    if par[2] - startparams[2] > epsilon:
        print "phase of AC oscillations wrong"
    if par[3] - startparams[3] > epsilon:
        print "offset of AC oscillations wrong"

    return np.all(par-startparams < epsilon), data1, fitfun, par, startparams

if __name__ == "__main__":
    #print test_Rabi_carrier_detailed()
    #print test_Ramsey_carrier_detailed()
    #print test_ACStark_detailed()
    [val, data1, fitfun, par, startparams] = test_ACStark_detailed()
    print val
    

# test_carrier.py ends here
