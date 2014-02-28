''' 
Calculate and compare fidelities for ideal, simulated, and experimental data.
'''

import numpy as np
import matplotlib.pyplot as pl
import scipy.stats as spst
import os, sys, shelve, time, copy

import PyTIQC.core.gates as U
import PyTIQC.core.simtools as sim
import readdata as rd
import densitymatrixreconstruction as dmr

pi = np.pi

###############################################################
# class EvaluateData
###############################################################

class EvaluateData:
    """ collect, analyze, compare data from experiment and simulation.
    Input: a dict of file names for experiment/simulation data. Files contain numpy arrays of density matrix vs pulse."
    """
    def __init__(self, files_dict=None, verbose=True):
        # indicators
        self.rho = None
        self.exploaded = False
        self.experrloaded = False
        self.simloaded = False

        if files_dict:
            self.files_dict = files_dict
            self.loadDataNPY(self.files_dict)
        self.verbose=verbose

        # 2 nested dictionaries to store all the fidelity information
        self.fid_dict = {}
        self.expfid_dict = {}

    def loadidealdata(self, rhoideal):
        ''' set the ideal densmat '''
        self.rho = rhoideal

    def loadexpdata(self, filepath, timetags, is_range, is_pop, print_fidelities = True):
        ''' load experimental data and optionally calculate state fidelities (jozsafid) '''
        auspath = rd.PathObject(filepath)
        self.YRexp = []
        self.rhoexp = []

        # load populations
        if is_pop:
            auspath = rd.PathObject(filepath)
            if is_range:
                dataobj_list = rd.ReadDataMultiple(timetags,is_range=True)
                dataobj = dataobj_list.pop()
                for item in dataobj_list:
                    dataobj += item
                print dataobj.parameters['cycles']
            else:
                dataobj = rd.ReadData(timetags[0])

            self.YRexp = dataobj.data_dict['cprb'][:,1:]
            self.YRexp_err = np.sqrt(self.YRexp*(1-self.YRexp) / dataobj.parameters['cycles'] )
            self.rhoexp = None

        # else load densmats  
        else:  
            data_list = rd.ReadDataMultiple(timetags, is_range=is_range)
            # if len(data) = len(ideal)-1, assume the missing one is the prep'd state
            if self.rho != None:
                if len(data_list) == np.shape(self.rho)[0]-1:
                    self.rhoexp.append(self.rho[0])
                elif len(data_list) != np.shape(self.rho)[0]:
                    print "length of data_list and ideal results don't match. Adding initial density matrix."
                    self.rhoexp.append(self.rho[0])
            else:
                print "ideal densmat not loaded - check len(rhoexp)."

            for dat in data_list: #range(len(expfiles)):
                rhoexp1 = dmr.IterML.iterfun(dat, 100)
                YRexp1 = np.diag(rhoexp1)
                self.rhoexp.append(rhoexp1)
                self.YRexp.append(YRexp1)

            self.YRexp = np.array(self.YRexp)
            self.rhoexp = np.array(self.rhoexp)
            if print_fidelities:
                for i in range(np.min([len(self.rhoexp), len(self.rho)])):
                    print i, ": ", U.jozsafid(self.rho[i], self.rhoexp[i])

        self.exploaded = True

    def loadsimdata(self, data, print_fidelities=False):
        ''' load simulation data directly from variables '''
        try: # data is a database object
            self.data = data
            self.data.RhoPNAll = np.array(self.data.RhoPNAll)
            self.rhosim = self.data.RhoPN
        except AttributeError: # or it's just the densmats
            # make a dummy data object
            self.data = sim.database(np.zeros(1), np.zeros(1), sim.hspace(1,2,0,0))
            self.data.RhoPNAll = data

        self.simloaded = True
        if print_fidelities:
            for i in range(np.min([len(self.rhosim), len(self.rho)])):
                print i, ": ", U.jozsafid(self.rho[i], self.rhosim[i])

    def loadexpdataNPY(self, files_dict):
        ''' load experiment data from NPY file '''
        if 'expdatafile' in files_dict.keys():
            try:
                self.rhoexp = np.load(files_dict['expdatafile'])
                self.exploaded=True
            except Exception, e:
                print "experimental data not loaded:", files_dict['expdatafile']
                print e    
        if 'expdataerrfile' in files_dict.keys(): 
            try:
                d = shelve.open(files_dict['expdataerrfile'])
                self.experr_dict = d['experr_dict']
                d.close()
                self.experrloaded=True
            except Exception, e:
                print "experimental error data not loaded:", files_dict['expdataerrfile']
                print e  


    def loadDataNPY(self, files_dict):
        ''' load ideal, exp, and sim data. Alternatively input directly the variables for the simulation: data, pulseseq, params. '''
        self.loadexpdataNPY(files_dict)

        if 'idealdatafile' in files_dict.keys():
            try:
                self.rho = np.load(files_dict['idealdatafile'])
            except Exception, e:
                print "ideal densmat not loaded:", files_dict['idealdatafile']
                print e 
        if 'simdatafile' in files_dict.keys():
            try:
                d = shelve.open(files_dict['simdatafile'])
                self.data = d['data']
                pulseseq = d['pulseseq']
                params = d['params']
                if 'dec' in d.keys():
                    self.dec = d['dec']
                self.data.RhoPNAll = np.array(self.data.RhoPNAll)
                #self.data.RhoPN = self.data.RhoPN/self.data.numruns
                self.rhosim = self.data.RhoPN
                d.close()
                self.simloaded = True
            except Exception, e:
                print "simulation data not loaded:", files_dict['simdatafile']
                print e       

    def getPopFromRho(self):
        ''' extract populations from density matrix '''
        try:
            self.YR = np.zeros((np.shape(self.rho)[0], np.shape(self.rho)[1]))
            for i in range(np.shape(self.rho)[0]):
                self.YR[i] = np.diag(self.rho[i])
        except Exception, e:
            print "self.YR not defined: ", e

        try:
            self.YRexp = np.zeros((np.shape(self.rhoexp)[0], np.shape(self.rhoexp)[1]))
            for i in range(np.shape(self.rhoexp)[0]):
                self.YRexp[i] = np.diag(self.rhoexp[i])
        except Exception, e:
            print "self.YRexp not defined: ", e
            
        try:
            self.YRsim = np.zeros((np.shape(self.rhosim)[0], np.shape(self.rhosim)[1]))
            for i in range(np.shape(self.rhosim)[0]):
                self.YRsim[i] = np.diag(self.rhosim[i])
        except Exception, e:
            print "self.YRsim not defined: ", e



    def calculate_sim_fidelities(self, fidelity, fidtype, startat=0, one=False, 
                                 grouperrors=False):
        ''' calculate fidelities from a simulation data set including error bars '''
        if not self.simloaded:
            print "ERROR: must load simulation data first"
            return 

        if one:
            rho = self.rho[startat:startat+2]
        else:
            rho = self.rho[startat:]
        data = self.data
        # create a dictionary as an entry in the top-level dict 
        self.fid_dict[fidtype] = {}

        fid_end_list = []
        for rhosim in data.RhoPNAll[:, -1]:
            fid_end_list.append( fidelity(rho[-1], rhosim) )
        self.fid_end_list = np.array(fid_end_list)
        fidmean = np.mean(fid_end_list)
        fiderr = np.std(fid_end_list)
        if self.verbose:
            print fidtype, np.around(fidmean,3), "+/-", np.around(fiderr,3)

        # fidelity vs pulse number
        fid_list = []
        fid_med = []   # median
        fid_75qtl = [] # 75 percent quartile
        fid_25qtl = [] # 25 percent quartile
        for i in range(np.shape(rho)[0]):
            fid_pulse = []
            if not grouperrors:
                for rhosim in data.RhoPNAll[:,i]:
                    fid_pulse.append( fidelity(rho[i], rhosim) )
            else:
                groups = int(np.ceil(np.shape(data.RhoPNAll)[0]/float(grouperrors)))
                for j in range(groups):
                    minind = j*grouperrors
                    maxind = np.min( [(j+1)*grouperrors, np.shape(data.RhoPNAll)[0]] )
                    rhoavg = np.mean(data.RhoPNAll[minind:maxind,i], axis=0)
                    fid_pulse.append( fidelity(rho[i], rhoavg) )
            fid_pulse = np.array(fid_pulse)

            fid_list.append(fid_pulse)
            fid_med.append(np.median(fid_pulse))
            fid_75qtl.append(spst.scoreatpercentile(fid_pulse, 75))
            fid_25qtl.append(spst.scoreatpercentile(fid_pulse, 25))
                
        self.fid_dict[fidtype]['list'] = np.array(fid_list)
        self.fid_dict[fidtype]['mean'] = np.mean(fid_list, axis=1)
        self.fid_dict[fidtype]['err'] = np.std(fid_list, axis=1)
        self.fid_dict[fidtype]['med'] = np.array(fid_med)
        self.fid_dict[fidtype]['75qtl'] = np.array(fid_75qtl)
        self.fid_dict[fidtype]['25qtl'] = np.array(fid_25qtl)

        # pad with 1's in front, if startat != 0
        fid1 = np.ones(np.shape(self.rho)[0] - np.shape(rho)[0])
        self.fid_dict[fidtype]['mean'] = \
            np.hstack([fid1, self.fid_dict[fidtype]['mean'] ])
        self.fid_dict[fidtype]['err'] = \
            np.hstack([0*fid1, self.fid_dict[fidtype]['err'] ])
        self.fid_dict[fidtype]['med'] = \
            np.hstack([fid1, self.fid_dict[fidtype]['med'] ])
        self.fid_dict[fidtype]['75qtl'] = \
            np.hstack([0*fid1, self.fid_dict[fidtype]['75qtl'] ])
        self.fid_dict[fidtype]['25qtl'] = \
            np.hstack([0*fid1, self.fid_dict[fidtype]['25qtl'] ])

        #if self.verbose:
        #    print fidtype,  "%0.3f + %0.3f - %0.3f" % ( fid_med[-1], fid_75qtl[-1]-fid_med[-1], fid_med[-1]-fid_25qtl[-1] )

    def calculate_exp_fidelities(self, fidelity, fidtype, startat=0):
        ''' calculate fidelities from an experimental data set '''
        if not self.exploaded:
            print "ERROR: must load experimental data first"
            return 

        rho = self.rho
        rhoexp = self.rhoexp

        self.expfid_dict[fidtype] = {}
        fid_exp = []
        for i in range(len(rho)):
            fid_exp.append( fidelity(rho[i], rhoexp[i]) )
        self.expfid_dict[fidtype]['fid'] = np.array(fid_exp)

        if self.experrloaded:
            self.expfid_dict[fidtype]['fiderr'] = self.experr_dict[fidtype]

    def calculatePlotFidelities(self, displ=1, plottitle='', startat=0, grouperrors=False, single=False, fidelity=None, fidtype=None):
        ''' plot sim/exp comparison for all fidelity measures.
        | displ: figure #, plottitle: title for plot
        | startat: the starting pulse in simulation; for ICF: ICFobj.gatepos
        | grouperrors: get stdev averaged over 100 points (set to 100) or False
        | single: make plot for a single fidelity metric
        | fidelity, fidtype: must be set to not None if single=True '''

        rho = self.rho[startat:]
        pslen = np.shape(rho)[0]

        # if no plotting, just calculate the fidelities and exit
        if displ == 0:
          for fidtype,fidelity in U.fidelities_dict.items():
            if self.simloaded:
                self.calculate_sim_fidelities(fidelity, fidtype, startat, 
                                              grouperrors=grouperrors)
            if self.exploaded:
                self.calculate_exp_fidelities(fidelity, fidtype, startat)
          return

        # dictionary to set plot order and titles
        fidplot_dict = { \
            'tracedist-rho': (1, "trace distance"),
            'jozsafid': (3, "fidelity"),
            'tracedist-pop': (2, "classical trace distance"),
            'sso': (4, "squared statistical overlap") }

        fig = pl.figure(displ)
        fig.clf()

        def makeSubplot(ax, fidelity, fidtype, title=True, xlabel=True, ylabel=True):
            ax.hold(True)
            # get fidelities of sim
            if self.simloaded:
                self.calculate_sim_fidelities(fidelity, fidtype, startat, 
                                              grouperrors=grouperrors)
                ax.plot(self.fid_dict[fidtype]['med'], 'b.-', label='sim')
                ax.errorbar(range(len(self.fid_dict[fidtype]['med'])),\
                    self.fid_dict[fidtype]['med'], \
                    yerr=[ \
                    self.fid_dict[fidtype]['med']-self.fid_dict[fidtype]['25qtl'], \
                    self.fid_dict[fidtype]['75qtl']-self.fid_dict[fidtype]['med'] ], \
                    fmt=None, color='b')
            
            # get fidelities of exp
            if self.exploaded:
                self.calculate_exp_fidelities(fidelity, fidtype, startat)
                fid_exp = self.expfid_dict[fidtype]['fid']
                if self.experrloaded:
                    ax.plot(fid_exp, 'r.-', label='exp')
                    experr = self.expfid_dict[fidtype]['fiderr']
                    if len(experr) == len(fid_exp)-1:
                        experr = np.hstack((0,experr))
                    ax.errorbar(range(len(fid_exp)), fid_exp, experr, fmt='r.-')
                else:
                    ax.plot(fid_exp, 'r.-', label='exp')
            
            # make the plot look nice
            # this block is just to set ylim nicely
            if self.simloaded and self.exploaded:
                ylimmin = np.min(np.hstack([self.fid_dict[fidtype]['med'],fid_exp]))
            elif self.simloaded:
                ylimmin = np.min(self.fid_dict[fidtype]['med'])
            elif self.exploaded:
                ylimmin = np.min(fid_exp)
            ax.set_ylim([ylimmin-0.05, 1])
            ax.set_xlim([0, np.shape(self.rho)[0]])
            if xlabel:
                ax.set_xlabel('pulse #')
            if ylabel:
                ax.set_ylabel('state fidelity')
            if title:
                ax.set_title(plottitle+' '+fidplot_dict[fidtype][1])
            #ax.legend(['sim', 'exp'], loc='lower left')
            ax.legend(loc='lower left')
            ax.hold(False)

        if single:
            if fidelity == None or fidtype == None:
                print "Error: must set fidelity and fidtype for single=True"
                return
            ax = fig.add_subplot(111)
            makeSubplot(ax, fidelity, fidtype, title=False)
        else:
            for fidtype,fidelity in U.fidelities_dict.items():
                subpl = fidplot_dict[fidtype][0]
                if subpl == 2 or subpl == 4: ylabel = False 
                else: ylabel=True
                if subpl == 1 or subpl == 2: xlabel = False 
                else: xlabel = True
                ax = fig.add_subplot(2,2,fidplot_dict[fidtype][0])
                makeSubplot(ax, fidelity, fidtype, True, xlabel,ylabel)

        fig.show()


    def calcAllFid(rho, rhoexp, displ=1):
        ''' calculate final fidelity using different measures '''
        fidelity_exp = np.zeros(len(rho))
        sso_exp = np.zeros(len(rho))
        popdist_exp = np.zeros(len(rho))

        for i in range(len(rho)):
            fidelity_exp[i] = U.fidelity(rho[i], rhoexp[i])
            sso_exp[i] = U.sso(rho[i], rhoexp[i])
            popdist_exp[i] = U.popdist(rho[i], rhoexp[i])

        fig = pl.figure(displ)
        fig.clf()
        ax = fig.add_subplot(111)
        ax.plot(fidelity_exp, '.-')
        ax.plot(sso_exp, '.-')
        ax.plot(popdist_exp, '.-')
        ax.set_ylim([np.min(np.hstack([fidelity_exp, sso_exp, popdist_exp]))-0.05, 1])
        ax.set_xlabel('pulse #')
        ax.set_ylabel('state fidelity')
        ax.legend(['fidelity', 'sso', 'popdist'], loc='lower left')
        fig.show()

        return fidelity_exp, sso_exp, popdist_exp

#########################################################
# functions to load raw cprb data from IBK experiments.
# possibly superceded by PyTIQC.evaluation.readdata.py
#########################################################

def loaddatacprb(datafile, NumberOfIterations = 100, populations=False):
    ''' load and calculate densitymatrix from cprb file '''

    file_obj = open(datafile)
    file_string = file_obj.read()
    file_obj.close()

    file_string = file_string.replace("\r\n", "\n")
    file_string = file_string.replace("\r", "\n")
    ftemp = open('tempcprb', 'wb')
    ftemp.write(file_string)
    ftemp.close()

    data = np.loadtxt('tempcprb', skiprows=1)
    os.remove('tempcprb')

    if populations:
        x = data[:,0]
        y = data[1:,0]
        yerr = np.sqrt(y*(1-y))/NumberOfIterations
        return x,y,yerr
    else:
        densitymatrix = dmr.IterML.iterfun(data, NumberOfIterations)
        return densitymatrix

def loadevaldata(filename, YRideal, Nqubits, NumIterations=100):
    ''' load population file and compare to ideal calculations '''
    data = np.loadtxt(filename, skiprows=1)
    popdata = np.zeros(2**Nqubits)

    y0 = np.zeros(2**Nqubits+1)
    y0[-1] = 1
    data = np.insert(data, 0, y0, axis=0)
    data = data[:,1:]

    experr = np.zeros(np.shape(data)[0])
    for k in range(len(data)):
        fid = np.sum(abs(YRideal[k,:] - data[k,:]) )
        #print ' at pulse', k+1, ' : ', str(np.around(fid,3))
        experr[k] = fid

    return data, experr


