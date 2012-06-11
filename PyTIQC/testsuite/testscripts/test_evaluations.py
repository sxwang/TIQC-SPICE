#!/usr/bin/env python
# -*- mode: Python; coding: latin-1 -*-
# Time-stamp: "2011-05-27 15:26:40 c704252"

#  file       test_evaluations.py
#  author     Thomas Monz 2011

import unittest
import PyTIQC.core.qctools as qc
import PyTIQC.core.simtools as sim
import numpy as np

import PyTIQC.evaluation.InvestigatePulseSeq as IPS
import PyTIQC.evaluation.densitymatrixreconstruction as dmr
import PyTIQC.evaluation.processtomography.proctom as proctom
import PyTIQC.evaluation.processtomography.quantumprocess as qproc

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

    def test_StateTomo(self):
        """
        do statetomo of SS and check whether we get the
        correct result
        """
        assert(test_StateTomo_detailed())

    def test_ProcTomo(self):
        """
        do statetomo of SS and check whether we get the
        correct result
        """
        assert(test_ProcTomo_detailed())



# Collect all test suites for running

all_suites = unittest.TestSuite((
  unittest.makeSuite(TestUserFunction)
  ))

plot_suites = all_suites

def run():
    unittest.TextTestRunner(verbosity=2).run(all_suites)


# ---------------------------

def test_StateTomo_detailed():
    NumberOfIons = 2
    PhononOverhead = 1

    hspace = sim.hspace(NumberOfIons,2,NumberOfIons+PhononOverhead,0)
    params = sim.parameters(hspace)
    dec = sim.decoherence(params)

    params.y0[qc.stateToIndex('SS,0', hspace)] = 1

    pulseseq = sim.PulseSequence( [ \
        sim.Rcar(params, 2*pi, 0, 0), \
        ] )

    ScanObj = IPS.ScanParameter_in_Sequence(pulseseq, params, dec, np.arange(3**NumberOfIons), type = 'StateTomo')

    ScanObj.runScan()
    data_dmr = ScanObj.output_dict['qstates_camera']
    rho = dmr.IterML.iterfun(data_dmr, 100)
    #if np.real(rho[3,3]) > 0.99:
    #    print 'statetomo passed'

    return np.real(rho[3,3]) > 0.99

def test_ProcTomo_detailed():
    NumberOfIons = 1
    PhononOverhead = 7

    hspace = sim.hspace(NumberOfIons,2,NumberOfIons+PhononOverhead,0)
    params = sim.parameters(hspace)
    dec = sim.decoherence(params)

    params.y0[qc.stateToIndex('S,0', hspace)] = 1

    pulseseq = sim.PulseSequence( [ \
        sim.Rcar(params, pi/2,0),
        ] )

    ScanObj = IPS.ScanParameter_in_Sequence(pulseseq, params, dec, np.arange(12**NumberOfIons), type = 'ProcTomo')

    ScanObj.runScan()
    data_proctom = ScanObj.output_dict['qstates_camera']
    chi = proctom.proctomo(data_proctom, 100)
    #if np.real(chi[0,0]) > 0.99:
    #    print 'proctomo passed'

    chiId = qproc.Unitary2Chi(pulseseq[0].Uidtr.conjugate())
    return np.max(abs(chi - chiId)) < 0.001


if __name__ == "__main__":
    print test_StateTomo_detailed()
    print test_ProcTomo_detailed()


# test_evaluations.py ends here
