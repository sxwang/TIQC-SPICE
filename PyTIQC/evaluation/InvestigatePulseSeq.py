#!/usr/bin/env python
# -*- mode: Python; coding: latin-1 -*-
# Time-stamp: "2011-08-30 20:36:33 c704252"

#  file       InvestigatePulseSeq.py
#  author     Thomas Monz

"""
general outline of this lib is to extend a
given pulse sequence (e.g. a CZ-CNOT) and
automatically add the required pulses
in the front/end to perform different forms
of tomographies (state, process, etc)
"""

import PyTIQC.core.qctools as qc
import PyTIQC.core.simtools as sim
import densitymatrixreconstruction as dmr
import processtomography.proctom as proctom

import numpy as np
import pp
mod = np.mod
pi = np.pi

import copy


def MSparity(phase, simparams, use_ideal=False):
    """ add a pi/2 pulse to MS gate, varying phase """
    total_append =[ sim.Rcar(simparams, pi/2, phase, -1)]
    total_append[0].use_ideal = use_ideal
    return total_append

# something seems to be screwed up with the ordering
# of ion counting ... left and right
# so i had to do that 'nuions-x' trick
def StateTomo(parameter, simparams, use_ideal=False):
    """
    pulses for a state tomography using addressed pulses.
    working with up to 4 qubits for the moment
    """
    # 0 - Z
    # 1 - X
    # 2 - Y
    # counting from the right, in principle ... :(

    total_append = []
    par = int(parameter)
    nuions = simparams.hspace.nuions
    # first qubit
    if (mod(par/3**0, 3) == 1): total_append.append(sim.Rcar(simparams, pi/2, 1.5*pi, 0))
    if (mod(par/3**0, 3) == 2): total_append.append(sim.Rcar(simparams, pi/2, pi, 0))

    # second qubit
    if (mod(par/3**1, 3) == 1): total_append.append(sim.Rcar(simparams, pi/2, 1.5*pi, 1))
    if (mod(par/3**1, 3) == 2): total_append.append(sim.Rcar(simparams, pi/2, pi, 1))

    # third qubit
    if (mod(par/3**2, 3) == 1): total_append.append(sim.Rcar(simparams, pi/2, 1.5*pi, 2))
    if (mod(par/3**2, 3) == 2): total_append.append(sim.Rcar(simparams, pi/2, pi, 2))

    # fourth qubit
    if (mod(par/3**3, 3) == 1): total_append.append(sim.Rcar(simparams, pi/2, 1.5*pi, 3))
    if (mod(par/3**3, 3) == 2): total_append.append(sim.Rcar(simparams, pi/2, pi, 3))

    if use_ideal:
        for pulse in total_append:
            pulse.use_ideal = True

    return total_append


def ProcTomoPrepare(parameter, simparams, use_ideal=False):
    """
    pulses for a proctomo (preparation part) using addressed pulses.
    working with up to 3 qubits for the moment
    """

    nuions = simparams.hspace.nuions
    total_prepend = []
    par = int(parameter)
    # first qubit
    if (mod(par/(4**0 * 3**nuions), 4) == 1): total_prepend.append(sim.Rcar(simparams, pi/2 , 0, nuions-1))
    if (mod(par/(4**0 * 3**nuions), 4) == 2): total_prepend.append(sim.Rcar(simparams, pi/2, pi*0.5, nuions-1))
    if (mod(par/(4**0 * 3**nuions), 4) == 3): total_prepend.append(sim.Rcar(simparams, pi, 0, nuions-1))

    # second qubit
    if (mod(par/(4**1 * 3**nuions), 4) == 1) and nuions >= 2 : total_prepend.append(sim.Rcar(simparams, pi/2 , 0, nuions-2))
    if (mod(par/(4**1 * 3**nuions), 4) == 2) and nuions >= 2 : total_prepend.append(sim.Rcar(simparams, pi/2, pi*0.5, nuions-2))
    if (mod(par/(4**1 * 3**nuions), 4) == 3) and nuions >= 2 : total_prepend.append(sim.Rcar(simparams, pi, 0, nuions-2))

    # third qubit
    if (mod(par/(4**2 * 3**nuions), 4) == 1) and nuions >= 3 : total_prepend.append(sim.Rcar(simparams, pi/2 , 0, nuions-3))
    if (mod(par/(4**2 * 3**nuions), 4) == 2) and nuions >= 3 : total_prepend.append(sim.Rcar(simparams, pi/2, pi*0.5, nuions-3))
    if (mod(par/(4**2 * 3**nuions), 4) == 3) and nuions >= 3 : total_prepend.append(sim.Rcar(simparams, pi, 0, nuions-3))

    if use_ideal:
        for pulse in total_prepend:
            pulse.use_ideal = True

    return total_prepend


def ProcTomoAnalyse(parameter, simparams, use_ideal=False):
    """
    pulses for proctomo (analysis part) using addressed pulses.
    working with up to 3 qubits for the moment
    """
    # 0 - Z
    # 1 - X
    # 2 - Y
    # counting from the right
    # i have to use strings here because neither sim nor params is
    # known here.
    nuions = simparams.hspace.nuions
    total_append = []
    par = int(parameter)
    # first qubit
    if (mod(par/3**0, 3) == 1): total_append.append(sim.Rcar(simparams, pi/2, pi*1.5, 0))
    if (mod(par/3**0, 3) == 2): total_append.append(sim.Rcar(simparams, pi/2, pi*1, 0))

    # second qubit
    if (mod(par/3**1, 3) == 1) and nuions >= 2: total_append.append(sim.Rcar(simparams, pi/2, pi*1.5, 1))
    if (mod(par/3**1, 3) == 2) and nuions >= 2: total_append.append(sim.Rcar(simparams, pi/2, pi*1, 1))

    # third qubit
    if (mod(par/3**2, 3) == 1) and nuions >= 3: total_append.append(sim.Rcar(simparams, pi/2, pi*1.5, 2))
    if (mod(par/3**2, 3) == 2) and nuions >= 3: total_append.append(sim.Rcar(simparams, pi/2, pi*1, 2))

    if use_ideal:
        for pulse in total_append:
            pulse.use_ideal = True

    return total_append


def simulateevolution(index, pulseseq, simparams, simdec):
    """ just a wrapper for qctools.simulateevolution, for pp. """
    data = qctools.simulateevolution(pulseseq, simparams, simdec)
    return data, index


class ScanParameter_in_Sequence:
    """
    | sequence ... input sequence that we want to investigate
    | simparams ... parameters of the simulator, such as hspace etc
    | parameter ... a parameter that we scan. it will affect what prepend/append will do to the sequence
    | prepend ... something prepended to the sequence, such as pulses for the process tomography preparation
    | append ... additional pulses, for instance to perform a state tomography
    """

    def __init__(self, sequence, simparams, simdec, parameter, type=None, verbose=False, save_all_data=True, doPP=False, use_ideal=False, pbar=None):
        self.simparams = simparams
        self.sequence = sequence
        self.simdec = simdec
        self.parameter = parameter
        self.type = type
        self.doPP = doPP
        self.verbose = verbose
        self.save_all_data = save_all_data
        self.use_ideal = use_ideal
        self.pbar = pbar

        self.type_dict = {
            'StateTomo': [None, StateTomo],
            'ProcTomo': [ProcTomoPrepare, ProcTomoAnalyse],
            'MSparity': [None, MSparity]
            }

        self.output_dict = {
            'qstates_camera': np.zeros((len(self.parameter),1+2**self.simparams.hspace.nuions)),
            'pmt_excitations': np.zeros((len(self.parameter),1+self.simparams.hspace.nuions+1)),
            'pmt_at_pulse': np.zeros((len(self.parameter), 1+self.simparams.hspace.nuions+1))
            }
        self.output_qstates_cam_all = \
              np.zeros((max([1,self.simdec.doRandNtimes]), len(self.parameter),1+2**self.simparams.hspace.nuions))
 

        self.sampledata = None
        self.alldata = []

        if doPP:
            self.simdec.doPP = False
            self.simdec.doPPprintstats = False
            self.simdec.progbar = False
        
    def get_full_sequence(self, index):
        """ prepend/append pulses to the given pulse sequence. """
        seq = copy.deepcopy(self.sequence)
        if self.type == None:
            return seq
        else:
            if not self.type_dict[self.type][0] == None:
                prepend_pulses = self.type_dict[self.type][0](index, self.simparams, self.use_ideal)
                seq.prepend(prepend_pulses)
            if not self.type_dict[self.type][1] == None:
                append_pulses = self.type_dict[self.type][1](index, self.simparams, self.use_ideal)
                seq.append(append_pulses)
            return seq


    def runSingle(self, index, param):
        """ run the simulation for one auto-generated pulse sequence """
        if self.verbose:
            print 'running parameter: '+str(param)
        self.simdec.progbar = False
        self.simdec.doPPprintstats = False
        # create the current seq
        pulseseq = self.get_full_sequence(param)

        dataobj = qc.simulateevolution(pulseseq, self.simparams, self.simdec)

        self.calc_output(dataobj, index)
        self.sampledata = dataobj
        if self.save_all_data:
            self.alldata.append(dataobj)

    def calc_output(self, data, index):
        """ process the data to obtain camera and PMT measurement results. """
        data.endpopulation(printstates=False)
        data.displayPMTpopulations(0)

        self.output_dict['qstates_camera'][index,0] = self.parameter[index]
        self.output_dict['qstates_camera'][index,1:] = data.YRend

        self.output_dict['pmt_excitations'][index,0] = self.parameter[index]
        self.output_dict['pmt_excitations'][index,1:] = data.P_PMT_end

        # does not work if pulseseq is empty ... remove for now
        if len(self.sequence) > 0:
            self.output_dict['pmt_at_pulse'][index,0] = self.parameter[index]
            self.output_dict['pmt_at_pulse'][index,1:] = data.P_PMT_P[-2]

        # store all runs in MonteCarlo 
        nRuns = np.shape(data.RhoPNAll)[0]
        if nRuns > 1:
            for i in range(nRuns):
                self.output_qstates_cam_all[i,index,0] = self.parameter[index]
                self.output_qstates_cam_all[i,index,1:] = np.diag(data.RhoPNAll[i,-1])
        else:
            self.output_qstates_cam_all[0] = self.output_dict['qstates_camera']

    def runBatch(self, parlist, pulseseq_list, params_list, dec_list):
        """ parallelize tomographies by sending simulations to server in batches """

        job_server = pp.Server( \
            ncpus = self.simparams.ppcpus,
            ppservers = self.simparams.ppservers,
            secret = self.simparams.ppsecret)

        jobs = [job_server.submit(simulateevolution, \
            args=(parlist[i], pulseseq_list[i], params_list[i], dec_list[i]), \
            depfuncs=(), \
            modules=('numpy','scipy', 'PyTIQC.core.simtools', 'PyTIQC.core.qmtools','PyTIQC.tools.progressbar','PyTIQC.core.qctools') ) for i in xrange(len(parlist)) ]

        for job in jobs:
            [dataobj, index] = job()
            if dataobj == None:
                print "simulationCore failed, exiting"
                return None
            self.calc_output(dataobj, index)
            self.sampledata = dataobj

        job_server.destroy()

    def runScan(self, batchsize=8):
      """ iterate over the parameters and run a simulation for each. """
      if not self.doPP:  
        for k in range(len(self.parameter)):
            self.runSingle(k, self.parameter[k])
      else:
        for k in range( int(np.ceil(len(self.parameter)/float(batchsize))) ):
            if k < len(self.parameter):
                parlist = self.parameter[k*batchsize:k*batchsize+batchsize]
            else:
                parlist = self.parameter[k*batchsize:]

            if self.pbar:
                self.pbar.update(int(min([k*batchsize, max(self.parameter)])*100./len(self.parameter)))
            if self.verbose:
                print 'running parameter: ', np.min(parlist), '-', np.max(parlist)

            params_list = []
            pulseseq_list = []
            dec_list = []
            data = None

            for p in parlist:
                # create the current seq
                pulseseq = self.get_full_sequence(p)

                pulseseq_list.append(copy.deepcopy(pulseseq))
                params_list.append(copy.deepcopy(self.simparams))
                dec_list.append(copy.deepcopy(self.simdec))

            self.runBatch(parlist, pulseseq_list, params_list, dec_list)
        if self.pbar:
            self.pbar.finish()

######################################################################
if __name__ == "__main__":
    import numpy as np
    import PyTIQC.core.simtools as sim
    import PyTIQC.core.qctools as qc

    pi = np.pi

    runStateTomo = False
    runStateTomoPP = False
    runProcTomo = False
    runProcTomoPP = True

    if runStateTomo or runStateTomoPP:
        NumberOfIons = 2
        PhononOverhead = 1

        hspace = sim.hspace(NumberOfIons,2,NumberOfIons+PhononOverhead,0)
        params = sim.parameters(hspace)
        dec = sim.decoherence(params)
        params.use_servers(['lindhard'])

        pulseseq = sim.PulseSequence( [ \
          sim.Rcar(params, 2*pi, 0, 0),
        ] )

        if runStateTomo:
            ScanObj = ScanParameter_in_Sequence(pulseseq, params, dec, np.arange(3**NumberOfIons), type='StateTomo', verbose=True)
        elif runStateTomoPP:
            ScanObj = ScanParameter_in_Sequence(pulseseq, params, dec, np.arange(3**NumberOfIons), type='StateTomo', verbose=True, doPP=True, use_ideal=True)

        ScanObj.runScan()

        data_dmr = ScanObj.output_dict['qstates_camera']
        rho = dmr.IterML.iterfun(data_dmr, 100)
        if np.real(rho[-1,-1]) > 0.99:
            print 'statetomo passed'
        else:
            print 'statetomo failed'

    if runProcTomo or runProcTomoPP:
        NumberOfIons = 2
        PhononOverhead = 1

        hspace = sim.hspace(NumberOfIons,2,NumberOfIons+PhononOverhead,0)
        params = sim.parameters(hspace)
        dec = sim.decoherence(params)
        params.use_servers(['local'])
        
        pulseseq = sim.PulseSequence( [ \
          sim.Rcar(params, 2*pi, 0, 0),
        ] )

        if runProcTomo:
            ScanObj = ScanParameter_in_Sequence(pulseseq, params, dec, np.arange(12**NumberOfIons), type='ProcTomo', verbose=True)
            ScanObj.runScan()
        elif runProcTomoPP:
            ScanObj = ScanParameter_in_Sequence(pulseseq, params, dec, np.arange(12**NumberOfIons), type='ProcTomo', verbose=True, doPP=True, use_ideal=True)
            ScanObj.runScan(batchsize=40)            

        data_proctom = ScanObj.output_dict['qstates_camera']
        chi = proctom.proctomo(data_proctom, 100)
        if np.real(chi[0,0]) > 0.99:
            print 'proctomo passed'
        else:
            print 'proctomo failed'
