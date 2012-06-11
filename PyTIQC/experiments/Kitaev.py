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

    def GeneratePulseSeq(self, params, Perm4, Perm2, Perm):
        ''' generate 4 pulse sequences given 3 permutations '''
        control90 = [sim.Delay(params, sim.Rac(params, pi/2, 0, 0).duration), \
                     sim.Rac(params, pi/2, 0, 0)]
        control45 = [sim.Delay(params, sim.Rac(params, pi/4, 0, 0).duration), \
                     sim.Rac(params, pi/4, 0, 0)]

        Hadamard0 = sim.PulseSequence( [ \
                sim.Rcar(params, pi/4, pi),
                sim.Rac(params, pi, 0),
                sim.Rcar(params, pi/4, 0) ])

        def PulseSeqWithControls(ctl):
            # 'ctl' is a bit array of control qubits, 2 in this case
            ctlstr = np.binary_repr(ctl).zfill(2)
            pulseseq = sim.PulseSequence([ \
                # 1st qubit
                Hadamard0,
                Perm4,
                Hadamard0,
                sim.Hide(params, 1, True),
                sim.Hide(params, 2, True),
                sim.MeasInit(params, 0),
                sim.Hide(params, 1, False),
                sim.Hide(params, 2, False),
                # 2nd qubit
                Hadamard0,
                Perm2, 
                copy.deepcopy(control90[ int(ctlstr[1]) ]),
                Hadamard0,
                sim.Hide(params, 1, True),
                sim.Hide(params, 2, True),
                sim.MeasInit(params, 0),
                sim.Hide(params, 1, False),
                sim.Hide(params, 2, False),
                # 3rd qubit
                Hadamard0,
                Perm, 
                copy.deepcopy(control90[ int(ctlstr[0]) ]),
                copy.deepcopy(control45[ int(ctlstr[1]) ]),
                Hadamard0,
                sim.Hide(params, 1, True),
                sim.Hide(params, 2, True),
                sim.MeasInit(params, 0, incl_hidingerr=False),                
                ])
            return pulseseq
        pulseseq_group = []
        for ctl in range(4):
            pulseseq_group.append( PulseSeqWithControls(ctl) )

        return pulseseq_group

    def getQFT(self, data_group=None, ind=None):
        ''' for 3 qubits, give output of QFT based on 3 classical measurements '''
        if data_group == None:
            data_group = self.data_group
        result = np.zeros(8)
        if ind == None: reg = data_group[0].register
        else: reg = data_group[0].registerAll[ind]
        result[0] = reg[0][0] * reg[1][0] * reg[2][0]
        result[4] = reg[0][0] * reg[1][0] * reg[2][1]
        if ind == None: reg = data_group[1].register
        else: reg = data_group[1].registerAll[ind]
        result[1] = reg[0][1] * reg[1][0] * reg[2][0]
        result[5] = reg[0][1] * reg[1][0] * reg[2][1]
        if ind == None: reg = data_group[2].register
        else: reg = data_group[2].registerAll[ind]
        result[2] = reg[0][0] * reg[1][1] * reg[2][0]
        result[6] = reg[0][0] * reg[1][1] * reg[2][1]
        if ind == None: reg = data_group[3].register
        else: reg = data_group[3].registerAll[ind]
        result[3] = reg[0][1] * reg[1][1] * reg[2][0]
        result[7] = reg[0][1] * reg[1][1] * reg[2][1]
        #print np.around(result,3)
        self.result = result

        return result

    def getQFTall(self, data_group=None):
        ''' get the QFT outputs for a series of simulations '''
        if data_group == None:
            data_group = self.data_group

        numruns = np.min([
                np.shape(data_group[0].registerAll)[0],
                np.shape(data_group[1].registerAll)[0],
                np.shape(data_group[2].registerAll)[0],
                np.shape(data_group[3].registerAll)[0] ])
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
            for ctl in range(4):
                print ctl
                data = qc.simulateevolution(pulseseq[ctl], params, dec)
                self.data_group.append(data)
                print np.around(data.register,3)
            result =  self.getQFT(self.data_group)

        else:
            job_server = pp.Server( \
                ncpus = 4, 
                ppservers = params.ppservers, 
                secret = params.ppsecret)

            self.data_group = [0,0,0,0]
            dec.doPPprintstats = True
            def jobfunc(x, pulseseq, params, dec):
                return qctools.simulateevolution(pulseseq[x], params, dec)
            controls = [0,1,2,3]
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
