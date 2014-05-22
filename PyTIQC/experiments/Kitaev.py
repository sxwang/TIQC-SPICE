import numpy as np
import PyTIQC.core.simtools as sim
import PyTIQC.core.qctools as qc
import PyTIQC.core.gates as U
import matplotlib.pyplot as pl
import copy, pickle
import pp

pi = np.pi

class Kitaev:
    """ performing Kitaev versions of QFT, order-finding, and Shor.
    currently restricted to 3 (virtual) ions only.
    | GeneratePulseSeq - generate group of pulseseq (4 total) given permutations
    | getQFT - process result of 4 sims to get 8-bit QFT result
    | simulateevolution - wrapper to repeat qc.simulateevolution 4 times
    """

    def __init__(self):
        self.data_group = []

    def GeneratePulseSeq(self, params, Perms, perm2=None, perm3=None):

        #Legacy support for old-style calls:
        if perm2 and perm3:
            Perms = [Perms,perm2,perm3]
        ''' generate pulse sequences given all permutations '''\

        self.nOps = len(Perms)
        controlRots = [ \
             [sim.Delay(params, sim.Rac(params, pi/2**(i+1), 0, 0).duration),
              sim.Rac(params, pi/2**(i+1), 0, 0)]
                    for i in range(self.nOps)]

        Hadamard0 = sim.PulseSequence( [ \
                sim.Rcar(params, pi/4, pi),
                sim.Rac(params, pi, 0),
                sim.Rcar(params, pi/4, 0) ])

        def PulseSeqWithControls(ctl):
            # 'ctl' is a bit array of control qubits
            ctlstr = np.binary_repr(ctl).zfill(self.nOps-1)
            pulseseq = sim.PulseSequence([])
            
            for k in range(self.nOps):
                pulseseq.append(Hadamard0)
                pulseseq.append(Perms[k])
                pulseseq += [copy.deepcopy(controlRots[i] \
                              [int(ctlstr[self.nOps-i-k-1])]) for i in range(k)]
                pulseseq.append(Hadamard0)

                if k == self.nOps-1:
                    pulseseq += [ \
                      sim.Hide(params, 1, True),
                      sim.Hide(params, 2, True),
                      sim.MeasInit(params, 0, incl_hidingerr=False) ] 
                else:
                    pulseseq += [ \
                      sim.Hide(params, 1, True),
                      sim.Hide(params, 2, True),
                      sim.MeasInit(params, 0),
                      sim.Hide(params, 1, False),
                      sim.Hide(params, 2, False) ]
            return pulseseq
        pulseseq_group = []
        for ctl in range(2**(self.nOps-1)):
            pulseseq_group.append( PulseSeqWithControls(ctl) )

        return pulseseq_group

    def getQFT(self, data_group=None, ind=None):
        ''' for 3 qubits, give output of QFT based on 3 classical measurements '''
        if data_group == None:
            data_group = self.data_group
        result = np.zeros(2**self.nOps)

        for k in range(2**self.nOps):
            if ind == None: reg = data_group[k%(2**(self.nOps-1))].register
            result[k] = 1
            bitstr = np.binary_repr(k).zfill(self.nOps)
            for i in range(self.nOps):
                result[k] *= reg[i][int(bitstr[self.nOps-i-1])]

        self.result = result
        return result


    def getQFTall(self, data_group=None):
        ''' get the QFT outputs for a series of simulations '''
        if data_group == None:
            data_group = self.data_group

        numruns = np.min([ \
                np.shape(data_group[k].registerAll)[0]
                for k in range(2**(self.nOps-1)) ])
        resultAll = []

        for run in range(numruns):
            result_single = self.getQFT(data_group, ind=run)
            result_single = result_single / sum(result_single)
            resultAll.append( result_single )

        resultAll = np.array(resultAll)

        self.numruns = numruns
        self.resultAll = resultAll
        self.result_mean = np.mean(resultAll, axis=0)
        self.result_std = np.std(resultAll, axis=0)

    def simulateevolution(self, pulseseq, params, dec, doPP=False):
        if not doPP:
            for ctl in range(2**(self.nOps-1)):
                print ctl
                data = qc.simulateevolution(pulseseq[ctl], params, dec)
                self.data_group.append(data)
                print np.around(data.register,3)
            result =  self.getQFT(self.data_group)

        else:
            job_server = pp.Server( \
                ncpus = 1, 
                ppservers = params.ppservers, 
                secret = params.ppsecret)

            self.data_group = [0 for _ in range(2**(self.nOps-1))]
            dec.doPPprintstats = True
            def jobfunc(x, pulseseq, params, dec):
                return PyTIQC.core.qctools.simulateevolution(pulseseq[x], params, dec)
            controls = range(2**(self.nOps-1))
            jobs = [(ctl, job_server.submit(jobfunc, \
                    args=(ctl, pulseseq, params, dec), \
                    modules=('numpy','scipy', 'PyTIQC.core.simtools', \
                             'PyTIQC.core.qctools', 'Kitaev',
                             'PyTIQC.core.qmtools', 'PyTIQC.core.sequel', \
                             'PyTIQC.tools.progressbar') ) )\
                         for ctl in controls ]
            for ctl, job in jobs:
                data1 = job()
                self.data_group[ctl] = data1

            result = self.getQFT(self.data_group)

        return self.result

    def datagroupsave(self, filename, RhoOnly=False):
        f = open(filename, 'wb')
        if RhoOnly:
            datagr = []
            for i in range(len(self.data_group)):
                datagr.append( self.data_group[i].RhoPNAll )
            pickle.dump(datagr, f)
        else:
            pickle.dump(self.data_group, f)
        f.close()
            

    def resultSave(self, filename):
        np.save(filename, self.resultAll)
