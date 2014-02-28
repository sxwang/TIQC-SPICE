""" In-circuit fidelity for full simulation. """

import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import FixedFormatter, FixedLocator
from scipy.optimize import leastsq, brentq
from scipy.interpolate import interp1d
import scipy.stats as spst
import scipy.linalg as splg
import os, sys, shelve, time, copy, shelve

import readdata as rd
import PyTIQC.tools.progressbar as progressbar
import PyTIQC.core.gates as U
import PyTIQC.core.simtools as sim
import PyTIQC.core.qctools as qc
import InvestigatePulseSeq as ips
import EvaluateData as evd
import processtomography.quantumprocess as qproc
import processtomography.proctom as proctom

pi = np.pi

############################################
# class InCircuitFidelity
############################################

class InCircuitFidelity:
    """ Evaluate the in-circuit fidelity of a gate within a pulse sequence given experimental data. """

    def __init__(self, pulseseq, params, evalobj=None, numtrials=16, doPP=True, verbose=True):
        self.pulseseq = pulseseq
        self.params = params
        self.evalobj = evalobj

        self.params.pplog = False  # avoid the "I/O operation on closed file" due to collisions?
        self.dephasetheoryfile = 'evaluation/dephasefidconv.shlv'
        self.dephasesinglefile = None

        # redefine dec to make sure all noise sources are initially off
        self.dec = sim.decoherence(params)
        self.dec.doRandNtimes = numtrials
        self.dec.doPP = doPP
        self.dec.doPPprintstats = False

        # set all pulses to be ideal, and later define non-ideal pulses
        for pulse in pulseseq:
            pulse.use_ideal = True

        self.verbose = verbose

        if not verbose:
            self.dec.progbar = False

        self.expfid = None

        self.dict_params = { 
          'addressing': [self._set_addressing, [0.01, 0.5], lambda x: [x, 0.8*x]],
          'dephase': [self._set_dephase, [0.00001, 0.14], \
                          lambda x: [0.00001, self.params.dephaseConv/x] ],
          'heating': [self._set_heatingrate,   [1, -5]],
          'intensity': [self._set_intensityfluct,  [-2, 1]],
          'initerr': [self._set_stateiniterr,      [-2,1]],
          'specmode': [self._set_specmodecoupling, [-2,1]],
          'spontdecay': [self._set_lifetime,       [1, -5]],
          'phaseoffset': [self._set_phaseOffset,   [0.0001, 1],   lambda x: [x]]
          }

    ### a set of functions to change the error parametrization
    def _set_addressing(self, new_params):
        self.params.addressingerr = new_params[0]
        self.params.addressingerr_global = new_params[1]
    def _set_dephase(self, new_params):
        self.params.correlationTime = new_params[0]
        self.params.coherenceTime = new_params[1]
    def _set_heatingrate(self, new_params):
        self.params.heatingrate = new_params[0]
    def _set_intensityfluct(self, new_params):
        self.params.intensityfluct = new_params[0]
    def _set_stateiniterr(self, new_params):
        self.params.stateiniterr = new_params[0]
    def _set_specmodecoupling(self, new_params):
        self.params.specmodecoupling = new_params[0]
    def _set_lifetime(self, new_params):
        self.params.lifetime = new_params[0]
    def _set_phaseOffset(self, new_params):
        self.params.phaseOffset = new_params[0]

    def use_measure(self, fidelity, fidtype):
        ''' define the fidelity measure to use for all error evaluations '''
        self.fidelity = fidelity
        self.fidtype = fidtype
        if self.verbose:
            print "using fidelity measure =", self.fidtype

    def use_error_metric(self, errormetric="last"):
        if errormetric == "mean" or "last" or "rel" or "one" or "wtail" or "wtail_rel":
            self.errormetric = errormetric
            if self.verbose: 
                print "using error metric = ", errormetric
        else:
            print "WARNING: invalid error metric: must be one of 'mean', 'last', 'rel', 'one', 'wtail' or 'wtail_rel'. defaulting to 'last'"
            self.errormetric = "last"

    def loadDataNPY(self, files_dict):
        ''' load experimental data and calculate state fidelities and total error '''
        self.files_dict = files_dict
        self.evalobj = evd.EvaluateData(self.files_dict, verbose=self.verbose)
        self.evalobj.calculate_exp_fidelities(self.fidelity, self.fidtype)

    def listGatetypes(self, printit=False):
        ''' list all gate types in a sequence.  
        Create a dictionary with gate tuple (type, theta, phi) as keys 
        and position in sequence as values.'''
        self.gates_dict = {}
        
        for i,pulse in enumerate(self.pulseseq.seq):
            key = (pulse.type, pulse.theta, pulse.phase)
            if key not in self.gates_dict:
                self.gates_dict[key] = [i]
            else:
                self.gates_dict[key].append(i)
        
        if printit:
            print "Available gates in sequence and their positions:"
            for (key, val) in self.gates_dict.iteritems():
                print "  ", self._printGatetype(key), ":", val

    def _printGatetype(self, gatetype):
        ''' convert the pi's in a gatetype tuple to human-readable format '''
        gatestr = '(' + gatetype[0] + ', '
        if int(pi/gatetype[1]) == pi/gatetype[1]:
            gatestr = gatestr + 'pi/' + str(int(pi/gatetype[1])) + ', '
        else:
            gatestr = gatestr + 'pi/' + str(np.around(pi/gatetype[1], 2)) + ', '
        if len(gatetype)==3 and gatetype[2] != 0:
            gatestr = gatestr + 'pi/' + str(int(pi/gatetype[2])) + ')'
        else:
            gatestr = gatestr + '0' + ')'
        return gatestr

    def setGatetype(self, gatetype, pos=-1):
        ''' set the gates of interest to have errors '''
        self.gatetype = gatetype
        self.gates = []
        for i,pulse in enumerate(self.pulseseq.seq):
            pulse.use_ideal = True
            # only match phase if given in the argument
            if len(gatetype) == 3:
                if pulse.type == gatetype[0] and \
                   pulse.theta == gatetype[1] and \
                   pulse.phase == gatetype[2]:
                   self.gates.append(i)
            else:
                if pulse.type == gatetype[0] and \
                   pulse.theta == gatetype[1]:
                   self.gates.append(i)


        if pos >= len(self.gates):
            print "warning: less than", pos+1, "gates in sequence. default to first appearance."
            pos = 0
            self.gatepos = self.gates[pos]
        elif pos >= 0:
            self.pulseseq.seq[self.gates[pos]].use_ideal = False
            self.gatepos = self.gates[pos]
        elif pos == -1:
            for ind in self.gates:
                self.pulseseq.seq[ind].use_ideal = False
            self.gatepos = self.gates[0]

        if self.verbose:
            print "matching gate", self._printGatetype(self.gatetype), "at pos", self.gatepos, "appearing", len(self.gates), "times at positions", self.gates

        if len(self.gates) == 0:
            print "warning: no gate of type", gatetype, "found"

    def setNoisetype(self, noisetype, indrange=None):
        ''' turn on the desired decoherence source '''
        self.noisetype = noisetype
        if noisetype in self.dec.dict.keys():
            self.dec.dict[noisetype] = True
            if self.verbose:
                print "matching noise:", self.noisetype
        else:
            print "noise type", noisetype, "not found in dec.dict, skipping"

        self.indrange = indrange

    def setNoiseAmp(self, noisetype, ind):
        ''' change the noise parametrization '''

        # user-defined values or use default
        if self.indrange == None:
            start = self.dict_params[noisetype][1][0]
            end = self.dict_params[noisetype][1][1]
        else:
            start = self.indrange[0]
            end = self.indrange[1]
        # convert ind to a parametrization factor by rescaling the interval.
        # this is used instead of np.linspace b/c ind can be non-integer.
        val = float(ind)/np.max([self.numwidths-1,1]) * (end-start) + start
        self.errorparam[int(ind)] = val
        
        new_param = []
        for j in range(len(self.dec.dict_params_static[noisetype])):
            new_param.append( self.dict_params[noisetype][2](val)[j] )
            
        #if self.noisetype == 'dephase':
        #    new_param[-1] = new_param[-1] * self.pulseseq[0].duration /(self.pulseseq[0].theta/2)
        self.dict_params[noisetype][0](new_param)

        if self.verbose:
            print 'set dec.dict_params[', noisetype, '] =', \
               np.around(new_param, 3)

    def randomizeevolution(self, numwidths, convtype=None, std=False, PTerrorbars=True):
        ''' parametrize noise by width '''

        self.numwidths = numwidths
        self.simfid = np.zeros(self.numwidths)
        self.simfiderr = np.zeros(self.numwidths)
        self.gatefid = np.zeros(self.numwidths)
        self.gatefid75 = np.zeros(self.numwidths)
        self.gatefid25 = np.zeros(self.numwidths)
        self.errorparam = np.zeros(self.numwidths) # converted to physical units

        tic = time.time()

        # go through the noise dictionary and tweak decoherence params
        for ind in range(self.numwidths):
            if self.noisetype == 'all':
                for noisetype in self.dec.dict_params_static.keys():
                    self.setNoiseAmp(noisetype, ind)
            else:
                self.setNoiseAmp(self.noisetype, ind)

            # re-initialize the pulseseq. BUG: not sure why it has to be done though; self.pulseseq.UtrAll gets set to [1] otherwise and pulse.UtrAll.append fails in simtools.
            #for pulse in self.pulseseq:
            #    pulse.UtrAll = []

            # run 
            [simerror, simerror_std] = self.run()
            self.simfid[ind] = simerror
            self.simfiderr[ind] = simerror_std

            # error to gate fidelity conversion
            if convtype == "analytic":
                [gatefidmed, gatefid75, gatefid25] = \
                    self.calculateGateFidelityAnalytic(ind=ind)
            elif convtype == "PT":
                [gatefidmed, gatefid75, gatefid25] = \
                    self.calculateGateFidelity(std=std)  # TEMP
                [gatefidmed, gatefid75, gatefid25] = \
                    self.calculateGateFidelityProcTomo(ErrorBars=PTerrorbars, ind=ind)
            else:
                [gatefidmed, gatefid75, gatefid25] = \
                    self.calculateGateFidelity(std=std)

            self.gatefid[ind] = gatefidmed
            self.gatefid75[ind] = gatefid75
            self.gatefid25[ind] = gatefid25

            if self.verbose:
                print "simerror = %0.4f +/- %0.4f" % (simerror, simerror_std)
                print "gate fidelity = %0.4f + %0.4f - %0.4f" % (gatefidmed, gatefid75, gatefid25)

        toc = time.time()
        if self.verbose:
            print "runtime:", toc-tic, "seconds"

    def randomizeOnce(self, indreal):
        ''' run once with the matching gate error '''
        if self.noisetype == 'all':
            for noisetype in self.dec.dict_params_static.keys():
                self.setNoiseAmp(noisetype, indreal)
        else:
            self.setNoiseAmp(self.noisetype, indreal)   
        [matcherror, matcherror_std] = self.run()
        self.matchfid = matcherror
        self.matchfiderr = matcherror_std

    def run(self):
        ''' run the simulation and calculate final error '''

        # get the partial pulseseq and starting state
        if self.errormetric == "one":
            self.pulseseq_part = sim.PulseSequence([self.pulseseq.seq[self.gatepos]])
        else:
            self.pulseseq_part = sim.PulseSequence(self.pulseseq.seq[self.gatepos:])
        # if using population measures, estimate input densmat
        if self.fidtype == 'jozsafid' or self.fidtype == 'tracedist-rho':
            self.params.use_rho0(self.evalobj.rhoexp[self.gatepos])
        else:
            pop = np.diag(self.evalobj.rhoexp[self.gatepos])
            rhoIn = np.diag(pop)
            for ii in range(len(pop)):
                for jj in range(len(pop)):
                    if ii != jj:
                        rhoIn[ii, jj] = np.sqrt(pop[ii]*pop[jj]) * \
                            np.exp(1j * np.angle(self.evalobj.rho[self.gatepos][ii,jj]))
            self.params.use_rho0(rhoIn)

        # run it
        pulseseq_inst = copy.deepcopy(self.pulseseq_part)
        self.data = qc.simulateevolution(pulseseq_inst, self.params, self.dec)
        self.data.RhoPNAll = np.array(self.data.RhoPNAll)

        self.evalobj.data = self.data
        self.evalobj.calculate_sim_fidelities(self.fidelity, self.fidtype, \
                                                  startat=self.gatepos, \
                                                  one = (self.errormetric == "one"))

        simerror, simerror_std = self.getICFerror()

        return simerror, simerror_std

    def getICFerror(self):
        ''' calculate total ICF error for sim and exp '''

        if self.expfid == None:
            # exp error list can be missing the first value
            if len(self.evalobj.expfid_dict[self.fidtype]['fid']) == \
               len(self.evalobj.expfid_dict[self.fidtype]['fiderr'])+1:
                self.evalobj.expfid_dict[self.fidtype]['fiderr'] = \
                    np.insert(self.evalobj.expfid_dict[self.fidtype]['fiderr'],0,0)

            if self.errormetric == "one":
                fidpart = self.evalobj.expfid_dict[self.fidtype]['fid']\
                                  [self.gatepos:self.gatepos+2]
                fidparterr = self.evalobj.expfid_dict[self.fidtype]['fiderr'] \
                                  [self.gatepos:self.gatepos+2]
            else:
                fidpart = self.evalobj.expfid_dict[self.fidtype]['fid'][self.gatepos:]
                fidparterr = self.evalobj.expfid_dict[self.fidtype]['fiderr'] \
                                  [self.gatepos:]

            if self.errormetric == "mean":
                self.expfid = np.sum(fidpart[0]-fidpart)/np.sum(len(fidpart))
                self.expfiderr = np.sqrt(np.sum(fidparterr**2))/np.sum(len(fidpart))
            elif self.errormetric == "last" or self.errormetric == "one":
                self.expfid = 1-fidpart[-1]
                self.expfiderr = fidparterr[-1]
            elif self.errormetric == "rel":
                self.expfid = 0
                self.expfiderr = 0
            elif self.errormetric == "wtail":
                weights = np.exp(-np.array(range(len(fidpart))))
                self.expfid = np.sum( (1-fidpart) * weights)
                self.expfiderr = np.sqrt(np.sum( fidparterr**2 * weights**2))
            elif self.errormetric == "wtail_rel":
                weights = np.exp(-np.array(range(len(fidpart))))
                self.expfid = np.sum( (fidpart[0]-fidpart) * weights)
                self.expfiderr = np.sqrt(np.sum( fidparterr**2 * weights**2))
            elif self.errormetric == "raw":
                self.expfid = 0
                self.expfiderr = 0
        
        if self.errormetric == "one":
            fidpart = self.evalobj.fid_dict[self.fidtype]['mean']\
                              [self.gatepos:self.gatepos+2]
            fidparterr = self.evalobj.fid_dict[self.fidtype]['err'] \
                              [self.gatepos:self.gatepos+2]
        else:
            fidpart = self.evalobj.fid_dict[self.fidtype]['mean'][self.gatepos:]
            fidparterr = self.evalobj.fid_dict[self.fidtype]['err'] \
                              [self.gatepos:]
        if self.errormetric == "mean":
            simerror = np.sum(fidpart[0]-fidpart)/np.sum(len(fidpart))
            simerror_std = np.sqrt(np.sum(fidparterr**2))/np.sum(len(fidpart))
        elif self.errormetric == "last" or self.errormetric == "one":
            simerror = 1-fidpart[-1]
            simerror_std = fidparterr[-1]
        elif self.errormetric == "rel":
            relfidlist = []
            for rhosim in self.evalobj.data.RhoPNAll[:, -1]:
                if self.errormetric == "rel":
                    relfidlist.append(self.fidelity(self.evalobj.rhoexp[-1], rhosim))
                else:
                    relfidlist.append(self.fidelity(self.evalobj.rhoexp[self.gatepos+1], rhosim))
            simerror = 1-np.mean(relfidlist)
            simerror_std = np.std(relfidlist)
        elif self.errormetric == "wtail":
            weights = np.exp(-np.array(range(len(fidpart))))
            simerror = np.sum( (1-fidpart) * weights)
            simerror_std = np.sqrt(np.sum( fidparterr**2 * weights**2))
        elif self.errormetric == "wtail_rel":
            weights = np.exp(-np.array(range(len(fidpart))))
            simerror = np.sum( (fidpart[0]-fidpart) * weights)
            simerror_std = np.sqrt(np.sum( fidparterr**2 * weights**2))
        elif self.errormetric == "raw":
            simerror = U.fidelities_dict[self.fidtype](self.evalobj.rhoexp[self.gatepos+1], self.evalobj.rhosim[self.gatepos+1])
            simerror_std = 0
            
        return simerror, simerror_std

    def calculateGateFidelityAnalytic(self, ind):
        ''' extract gate fidelity from analytical formula if available. 
            use process fid instead of average fid'''
        x = self.errorparam[ind]
        nuions = self.params.hspace.nuions
        strgatetype = str((self.gatetype[0], self.gatetype[1]))

        if self.noisetype != 'dephase':
            print "no analytical formula defined for ", self.noisetype
            return self.calculateGateFidelity(std=False)

        if self.dephasesinglefile != None:
            dephasecsv = np.loadtxt(self.dephasesinglefile, delimiter=',')
        else:
            try:
                d = shelve.open(self.dephasetheoryfile)
                dephasecsv = d[str(nuions)+strgatetype]
            except KeyError, e:
                print "Exception occurred in gate fidelity conversion: ", e
                return self.calculateGateFidelity(std=False)

        fidphase = lambda x: np.interp(x, dephasecsv[:,0], dephasecsv[:,1])
                
        return fidphase(x), fidphase(x), fidphase(x) # no errorbars yet

    def calculateGateFidelity(self, std):
        ''' extract gate fidelity from the *current* stored unitaries '''
        # note!!! self.data.pulseseq not self.pulseseq
        gatefid = []
        for pulse in self.data.pulseseq: 
            gatefid_pulse = []
            if not pulse.use_ideal:
              for i in range(self.data.numruns):
                chi_id = qproc.Unitary2Chi(pulse.Uidtr)
                chi_sim = qproc.Unitary2Chi(pulse.UtrAll[i])
                #print "Unitaries:"
                #print np.around(chi_id,3) # TEMP
                #print np.around(chi_sim,3)
                gatefid_pulse.append(U.fidelity(chi_id, chi_sim))
              gatefid.append(gatefid_pulse)

        gatefid = np.array(gatefid).flatten()
        self.gatefidsample = gatefid

        fid_mean = np.mean(gatefid)
        fid_std = np.std(gatefid)
        fid_med = np.median(gatefid)
        fid_75qtl = spst.scoreatpercentile(gatefid, 75)
        fid_25qtl = spst.scoreatpercentile(gatefid, 25)

        if std:
            return fid_mean, fid_mean+fid_std, fid_mean-fid_std
        else:
            return fid_med, fid_75qtl, fid_25qtl

    def calculateGateFidelityProcTomo(self, ErrorBars=False, ind=None):
        ''' use proctomo to convert error of gate to fidelity '''
        #print "Using ProcTomo to calculate gate fidelity. This may take a while ..."

        if ind != None:
            indreal = [ind]
            if ErrorBars:
                indreal = [ind, ind, ind]
        elif not ErrorBars:
            indreal = [self.error0]
        else:
            indreal = [self.error0, self.error0max, self.error0min]

        result = []
        for ind in indreal:
            if self.noisetype == 'all':
                for noisetype in self.dec.dict_params_static.keys():
                    self.setNoiseAmp(noisetype, ind)
            else:
                self.setNoiseAmp(self.noisetype, ind)   

            if not self.verbose:
                widgets = [progressbar.Percentage(), ' ', progressbar.Bar(),' ', progressbar.ETA()]
                pbar = progressbar.ProgressBar(widgets=widgets).start()        
            else:
                pbar = None

            pos = self.gates[0] # only do it once, for now
            newpulseseq = sim.PulseSequence([ copy.deepcopy(self.pulseseq[pos]) ])
            self.pulseseq[pos].UtrAll = []
            newpulseseq[0].use_ideal=False

            ScanObj = ips.ScanParameter_in_Sequence(newpulseseq, self.params, self.dec, np.arange(12**self.params.hspace.nuions),type='ProcTomo', verbose=self.verbose, doPP=True, use_ideal=True, pbar=pbar)
            ScanObj.runScan(batchsize=40)
            data_proctom = ScanObj.output_dict['qstates_camera']
            chi = proctom.proctomo(data_proctom, 100)
            chiId = qproc.Unitary2Chi(newpulseseq[0].Uidtr.conjugate())
            result.append(U.jozsafid(chiId, chi))
            #print "PT:"
            #print np.around(chiId,3)
            #print np.around(chi,3) # TEMP

        if ind != None and ErrorBars:
            result.sort()
            result.append(result[0])
            result.remove(result[0])

        if not ErrorBars:
            print self._printGatetype(self.gatetype),":", np.around(result[0],4)
            return [result[0], result[0], result[0]]
        else:
            print self._printGatetype(self.gatetype) + " : %0.4f + %0.4f - %0.4f" % (np.around(result[0],4), np.around(result[1]-result[0],4), np.around(result[0]-result[2],4) )
            return result



    def fitICFerror(self, displ=1):
        ''' match simulated error to experiment '''
        # TODO: for errormetric="rel", need to find max, not implemented yet

        if len(self.gates) == 0: return

        x = self.errorparam
        simfid = self.simfid
        simfiderr = self.simfiderr

        xx = np.linspace(np.min(x), np.max(x))
        p0 = np.polyfit(x, simfid, 3)
        fittedfunc = lambda x: np.polyval(p0, x)
        fittedfunc_0 = lambda x: fittedfunc(x) - self.expfid

        # if simfiderr are small then use these
        simerrormin = simfid - simfiderr
        simerrormax = simfid + simfiderr
        p1 = np.polyfit(x, simerrormin, 3)
        fittedfunc_min = lambda x: np.polyval(p1,x)
        p2 = np.polyfit(x, simerrormax, 3)
        fittedfunc_max = lambda x: np.polyval(p2,x)
        # fittedfunc_min/max or fittedfunc
        fittedfunc_min_0 = lambda x: fittedfunc(x) - (self.expfid+self.expfiderr)
        fittedfunc_max_0 = lambda x: fittedfunc(x) - (self.expfid-self.expfiderr)

        try: 
            self.error0 = brentq(fittedfunc_0, np.min(x), np.max(x))
            try:
                self.error0min = brentq(fittedfunc_min_0, np.min(x), np.max(x))
            except ValueError, e:
                self.error0min = self.error0
            try:
                self.error0max = brentq(fittedfunc_max_0, np.min(x), np.max(x))
            except ValueError, e:
                self.error0max = self.error0

            fiderrorbars = (self.error0max - self.error0min) /2
            if self.verbose:
                print 
                print "FIT SUMMARY:"
                print "exp error = %0.4f +/- %0.4f" % (self.expfid, self.expfiderr)
                print "parametrized matching error of gate = %0.4f (+ %0.4f | - %0.4f)" % (self.error0, self.error0max-self.error0, self.error0-self.error0min)
        except ValueError, e:
            print "experror: %0.4f" % self.expfid
            print "ERROR: fidelity matching failed:", e
            self.error0 = -1
            self.error0min = -1
            self.error0max = -1

        if displ:
            fig = pl.figure(displ)
            fig.clf()
            ax4 = fig.add_subplot(111)
            # fit to simulation, upper & lower bound
            ax4.plot(xx, fittedfunc(xx), 'g-', linewidth=2)
            #ax4.plot(xx, fittedfunc_max(xx), 'g-', linewidth=1)
            #ax4.plot(xx, fittedfunc_min(xx), 'g-', linewidth=1)
            # simulated values & errorbars
            ax4.plot(self.errorparam, self.simfid, 'b.')
            ax4.errorbar(self.errorparam, self.simfid, self.simfiderr, fmt=None, ecolor='b')
            # experiment, upper & lower bound
            ax4.plot([0,np.max(x)], [self.expfid, self.expfid], 'r', linewidth=2)
            ax4.plot([0,np.max(x)], [self.expfid+self.expfiderr, self.expfid+self.expfiderr], 'r', linewidth=1)
            ax4.plot([0,np.max(x)], [self.expfid-self.expfiderr, self.expfid-self.expfiderr], 'r', linewidth=1)
            # labels and such
            ax4.set_ylim([0,np.max(self.simfid)+0.05])
            ax4.set_xlabel('Error parametrization')
            ax4.set_ylabel('total error')
            ax4.set_title('ICF: error vs param')
            fig.show()

    def convertErrorToGateFid(self, displ=2):
        ''' fit the simulated gate fidelities to error param and find matching gate fidelity '''

        if len(self.gates) == 0: return -1, -1, -1
        if self.error0 == -1: return -1, -1, -1

        x = self.errorparam
        xx = np.linspace(np.min(x), np.max(x))
        p = np.polyfit(self.errorparam, self.gatefid, 3)
        self.gatefidfunc = lambda x: np.polyval(p, x)

        ylimmin = np.min(self.gatefid25)-0.02

        try:
            gatefidmid = self.gatefidfunc(self.error0)
            gatefidmin = self.gatefidfunc(self.error0min)
            gatefidmax = self.gatefidfunc(self.error0max)
            gatefiderr = (gatefidmax-gatefidmin)/2.

            if gatefidmax == gatefidmid: gatefidmax = 1. # upper bound on error
            print "gate fidelity obtained from fitting to error param:"
            print self._printGatetype(self.gatetype), \
                ": %0.4f (+ %0.4f | - %0.4f)" \
                %(gatefidmid, gatefidmax-gatefidmid, gatefidmid-gatefidmin)

        except ValueError, e:
            print "ERROR: function evaluation failed: ", e
            gatefidmid = np.min(self.gatefidfunc(x))
            gatefiderr = 0
            print self._printGatetype(self.gatetype), ": (lower bound) %0.4f" %(gatefidmid)

        if displ:
            fig = pl.figure(displ)
            fig.clf()
            ax5 = fig.add_subplot(111)
            # fit to simulation
            ax5.plot(xx, self.gatefidfunc(xx), 'g-', linewidth=2)
            # simulated values
            ax5.errorbar(self.errorparam, self.gatefid, yerr=[ \
                    self.gatefid-self.gatefid25, \
                    self.gatefid75-self.gatefid], fmt='b.')
            # experiment, upper & lower bount
            ax5.plot([self.error0,self.error0], [ylimmin,1],'r',linewidth=2)
            ax5.plot([self.error0min,self.error0min], [ylimmin,1],'r',linewidth=1)
            ax5.plot([self.error0max,self.error0max], [ylimmin,1],'r',linewidth=1)
            ax5.set_ylim([ylimmin, 1])
            ax5.set_xlabel('gate error parametrization')
            ax5.set_ylabel('gate fidelity')
            ax5.set_title('Fidelity vs error parametrization')

            fig.show()

        return gatefidmid, gatefidmax-gatefidmid, gatefidmid-gatefidmin

    def plotSimvsGatefid(self, displ=3):

        p = np.polyfit(self.gatefid, self.simfid, 3)
        gatefidfunc = lambda x: np.polyval(p, x)
        xx = np.linspace(np.min(self.gatefid), np.max(self.gatefid))
        gatefidfunc_0 = lambda x: gatefidfunc(x) - self.expfid
        gatefidfunc_min_0 = lambda x: gatefidfunc(x) - (self.expfid + self.expfiderr)
        gatefidfunc_max_0 = lambda x: gatefidfunc(x) - (self.expfid - self.expfiderr)

        try:
            self.error1 = brentq(gatefidfunc_0, np.min(self.gatefid), 
                                 np.max(self.gatefid))
            try:
                self.error1min = brentq(gatefidfunc_min_0, 
                                 np.min(self.gatefid), np.max(self.gatefid) )
            except ValueError, e:
                self.error1min = self.error1
            try:
                self.error1max = brentq(gatefidfunc_max_0, 
                                 np.min(self.gatefid), np.max(self.gatefid) )
            except ValueError, e:
                self.error1max = self.error1

            fiderrorbars = (self.error1max - self.error1min) /2
            if self.error1max == self.error1: 
                self.error1max = 1  # upper bound on fidelity

            print "gate fidelity obtained from fitting after gate fidelity conversion:"
            print self._printGatetype(self.gatetype), \
                ": %0.4f (+ %0.4f | - %0.4f)" \
                %(self.error1, self.error1max-self.error1, self.error1-self.error1min)
        except Exception, e:
            print "ERROR: fitting failed: ", e
            self.error1 = -1
            self.error1min = -1
            self.error1max = -1

        if displ:
            fig = pl.figure(displ)
            fig.clf()
            ax = fig.add_subplot(111)
            ax.errorbar(self.gatefid, self.simfid, self.simfiderr, fmt='b.',
                        label="simulation error")
            ax.plot(xx, gatefidfunc(xx), 'g-', linewidth=2, 
                        label="simulation fit")
            ax.plot([np.min(self.gatefid),1], [self.expfid, self.expfid], 'r', linewidth=2, label="experiment error")
            ax.plot([np.min(self.gatefid),1], [self.expfid+self.expfiderr, self.expfid+self.expfiderr], 'r', linewidth=1)
            ax.plot([np.min(self.gatefid),1], [self.expfid-self.expfiderr, self.expfid-self.expfiderr], 'r', linewidth=1)
            ax.set_xlabel('gate fidelity')
            ax.set_ylabel('error metric of pulse sequence')
            #ax.set_title('ICF for: ' + self.noisetype)
            ax.legend()
            ax.set_xlim([np.min(self.gatefid)-0.02, 1])
            ax.set_ylim([0, np.max(self.simfid)+np.max(self.simfiderr)+0.02])
            fig.show()

        return self.error1, self.error1max-self.error1, self.error1-self.error1min


    def plotevolution(self, legends=False, displ=4, randomize=True):
        ''' make plots of populations at end of each pulse '''
        # if randomized, generate data from exp-sim matched errors. Otherwise use the most recent simulation results.
        if randomize:
            self.randomizeOnce(self.error0)

        fig = pl.figure(displ)
        fig.clf()

        legendstr = []
        for i in range(2**self.params.hspace.nuions):
            legendstr.append(np.binary_repr(i).zfill(self.params.hspace.nuions))

        YRideal = np.zeros_like(self.data.YRPN)
        YRexp = np.zeros_like(self.data.YRPN)
        for i in range(len(self.pulseseq)+1):
            YRideal[i] = np.diag(self.evalobj.rho[i])
            YRexp[i] = np.diag(self.evalobj.rhoexp[i])

        def drawlines(ax, ymax):
          for i in range(4):
            i += 1
            ax.plot([i,i],[0, ymax],'0.9')

        ax1 = fig.add_subplot(311)
        ax1.plot(YRideal)
        if legends: ax1.legend(legendstr)
        #drawlines(ax1, 0.5)
        ax1.set_ylim([0,np.max(YRideal)+0.02])
        ax1.set_ylabel('state populations')
        ax1.set_title('Ideal state evolution')

        ax2 = fig.add_subplot(312)
        ax2.plot(YRexp)
        if legends: ax2.legend(legendstr)
        #drawlines(ax2, 0.5)
        ax2.set_ylim([0,np.max(YRideal)+0.02])
        ax2.set_ylabel('state populations')
        ax2.set_title('Experimental state evolution')

        ax3 = fig.add_subplot(313)
        ax3.plot(self.data.YRPN)
        if legends: ax3.legend(legendstr)
        #drawlines(ax3, 0.5)
        ax3.set_ylim([0,np.max(YRideal)+0.02])
        ax3.set_xlabel('pulse number')
        ax3.set_ylabel('state populations')
        ax3.set_title('Simulated state evolution')

        fig.show()
