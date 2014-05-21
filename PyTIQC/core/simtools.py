''' this module contains all the key class definitions for the simulation: hspace, decoherence, parameters, pulse, database. '''

import numpy as np
import numpy.matlib as npml
import scipy.misc as smisc
import scipy.linalg as splg
import scipy.constants as constants
import pylab as pl
from mpl_toolkits.mplot3d import Axes3D

import qmtools as qm
import gates as U

import PyTIQC.core
import copy
import cPickle as pickle
import time
import os

def nchoosek(n,k):
    return smisc.comb(n,k)

pi = np.pi

############################################################
# class hspace
############################################################
class hspace:
  """ hilbert space. a struct containing:

  | nuions -- number of ions (default = 1)
  | levels -- ion levels (default = 2, can be 3)
  | maxphonons -- default = 0
  | densmat -- use density matrix formalism. not implemented yet
  """

  def __init__(self, nuions=1, levels=2, maxphonons=0, densmat=0):
    if levels != 2 and levels != 3:
        raise ValueError('number of levels must be either 2 or 3')

    self.nuions = nuions 
    self.levels = levels
    self.maxphonons = maxphonons
    self.densmat = 0
    self.dimensions = (self.maxphonons+1) * self.levels ** self.nuions
    self.visible = np.ones(self.dimensions)
    self.calculate_operators()
    self.initial_state()

  def initial_state(self, state="ground", nBar=0, qstate=None):
    """ define an initial state """

    self.y0 = np.zeros(self.dimensions, np.complex128)

    # default: ground state, p(S...S,0) = 1, all else = 0
    if state == "ground":
        self.y0[qm.stateToIndex(self.nuions*'S'+',0', self)] = 1

    # thermal state: normalized poisson distribution
    if state == "thermal":
        # distribution
        def pThermal(n, nBar):
            return float(nBar**n) / (nBar + 1)**(n+1)
        # set the population of each state
        for n in range(self.maxphonons+1):
            self.y0[qm.stateToIndex(self.nuions*'S'+',0',self)+n] = \
                pThermal(n, nBar)
        # normalize, then convert to amplitudes
        self.y0 = self.y0 / np.sum(self.y0)
        self.y0 = np.sqrt(self.y0)

    # quantum state: specify a state string
    if state == "quantum":
        if qstate==None or len(qstate) != self.nuions:
            raise ValueError('Invalid value for qstate: must be string of length=nuions')
        self.y0[qm.stateToIndex(qstate+',0', self)] = 1

  def calculate_operators(self):
    """ define a set of standard operators (a, a-dagger, sig_z, sig_m, sig_p etc) useful later for constructing hamiltonians. """
    self.operator_dict = {}

    # zero matrix
    self.operator_dict['zero'] = np.zeros((self.dimensions,self.dimensions))

    # projection matrix (for spontaneous decay)
    if self.levels == 2:
        proj = np.array([[0,0],[1,1]])
    elif self.levels == 3:
        proj = np.array([[0,0,0], [1,1,0], [0,0,1]]) 
    self.operator_dict['proj'] = proj

    # projection matrix into state S
    if self.levels == 2:
        projs = np.array([[0,0],[0,1]])
    elif self.levels == 3:
        projs = np.array([[0,0,0], [0,0,0], [0,0,1]]) 
    self.operator_dict['projs'] = projs

    # creation/annihilation operator for phonons
    a = np.zeros((self.maxphonons+1, self.maxphonons+1))
    for k in xrange(self.maxphonons):
      a[k,k+1] = k+1
    a = np.sqrt(a)
    a_dagger = a.transpose()

    self.operator_dict['a'] = a
    self.operator_dict['a_dag'] = a_dagger
    self.operator_dict['id_a'] = np.eye(self.maxphonons+1)

    # raising/lowering for qubit
    id = np.diag(np.ones(self.levels))
    if self.levels == 2:
        sigm = np.array([[0,1],[0,0]])
    elif self.levels == 3:
        sigm = np.array([[0,1,0],[0,0,0],[0,0,0]])
    # let's store lowering as a matrix of the form [:,:,qubit]
    # add one last element, which is the sum over all operators
    # which might be useful for one or the other thing
    lev = self.levels
    sigm_matrix = np.zeros((lev**self.nuions, lev**self.nuions, self.nuions+1))
    sigp_matrix = np.zeros((lev**self.nuions, lev**self.nuions, self.nuions+1))
    for k in xrange(self.nuions):
      a = 1
      for l in xrange(self.nuions):
        if l==k:
          a = npml.kron(a,sigm)
        else:
          a = npml.kron(a, id)
      sigm_matrix[:,:,self.nuions-1-k] = a
      sigp_matrix[:,:,self.nuions-1-k] = a.transpose()

    sigm_matrix[:,:,-1] = np.sum(sigm_matrix, axis = 2)
    sigp_matrix[:,:,-1] = np.sum(sigp_matrix, axis = 2)

    self.operator_dict['lowering'] = sigm_matrix
    self.operator_dict['raising'] = sigp_matrix

    # raising/lowering on the auxiliary level for qubit - same as above
    id = np.diag(np.ones(self.levels))
    if self.levels == 3:
        sigmx = np.array([[0,0,0],[0,0,1],[0,0,0]])
        lev = self.levels
        sigmx_matrix = np.zeros((lev**self.nuions, lev**self.nuions, self.nuions+1))
        sigpx_matrix = np.zeros((lev**self.nuions, lev**self.nuions, self.nuions+1))
        for k in xrange(self.nuions):
          a = 1
          for l in xrange(self.nuions):
            if l==k:
              a = npml.kron(a,sigmx)
            else:
              a = npml.kron(a, id)
          sigmx_matrix[:,:,self.nuions-1-k] = a
          sigpx_matrix[:,:,self.nuions-1-k] = a.transpose()

        sigmx_matrix[:,:,-1] = np.sum(sigmx_matrix, axis = 2)
        sigpx_matrix[:,:,-1] = np.sum(sigpx_matrix, axis = 2)

        self.operator_dict['lowering_aux'] = sigmx_matrix
        self.operator_dict['raising_aux'] = sigpx_matrix

    # now, finally, sigz
    if self.levels == 2:
        sigz = np.array([[1,0],[0,-1]])
    elif  self.levels == 3:
        sigz = np.array([[1,0,0],[0,-1,0],[0,0,0]])
    sigz_matrix = np.zeros((lev**self.nuions, lev**self.nuions, self.nuions+1))
    for k in xrange(self.nuions):
      a = 1
      for l in xrange(self.nuions):
        if l==k:
          a = npml.kron(a,sigz)
        else:
          a = npml.kron(a, id)

      sigz_matrix[:,:,self.nuions - 1 -k] = a

    sigz_matrix[:,:,-1] = np.sum(sigz_matrix, axis = 2)

    self.operator_dict['sigz'] = sigz_matrix


############################################################
# class decoherence
############################################################
class decoherence:
  """ decoherence effects to be included in Hamiltonian."""

  def __init__(self, params):
    # general decoherence parameters
    self.doRandNtimes = 0 # Monte-Carlo: average over this many instances; if not zero, start averaging
    self.doRandBatch = 40 # save the data to file after every batch runs
    self.doPP = False
    self.doPPprintstats = True
    self.doSQL = False
    self.doSQLname = None

    self.params = params

    self.progbar = True
    self.intensfluctV = []
    self.spontdecayV = []
    self.heatingV = []

    self.dict = {
      'none': True,
      'all': False, # turn everything on
      'addressing': False, # addressing error - constant during all evolution
      'dephase': False, # dephasing, should be replaced by dephasing_laser
#      'dephase_Bfield': False, # dephasing due to magnetic field fluctuations
      'heating': False, # heating
      'hiding': False, # hiding error
      'intensity': False, # intensity fluctuations
      'initerr': False, # qubit initialisation error, assumes that you accidently prepared S'
         # will roll dice and, if true, set the coupling to all lasers to zero
         # via a modified addressing error matric
      'specmode': False, # spectator mode coupling, modelled as an initialized intensity fluctuation (changes coupling globally for a single run)
      'spontdecay': False, # spont. emission from the qubit
      'switchtime': False, # switching time between pulses
      'phaseoffset': False # constant phase offset on each pulse
      }

    # mapping between dec.dict and parameters
    # the pair of numbers is the range for np.logspace(base=2) to change the default params for ICF
    # [start, end]; 2**start is smaller error than 2**end
    self.dict_params_static = { 
      'addressing': [self.params.addressingerr, self.params.addressingerr_global],
      'dephase': [self.params.correlationTime, self.params.coherenceTime],
      'heating': [self.params.heatingrate],
      'intensity': [self.params.intensityfluct],
      'initerr': [self.params.stateiniterr],
      'specmode': [self.params.specmodecoupling],
      'spontdecay': [self.params.lifetime],
      'hiding': [self.params.hidingMeaserr, self.params.hidingerr]
      }

  def calcDephasing(self, T, stepsize):
    """ construct a time-correlated gaussian series as self.dephase.
        Tau and timeperpulse has units of real time (us). 
        For details see "How to generate exponentially correlated Gaussian
        random numbers" by Markus Deserno (UCLA 2002) """

    # get the two time parameters
    # note!!! the factor params.dephaseConv is measured "experimentally" 
    # (i.e. by running simulations)
    # the analytical result is 2log(2) -  see Chwalla's thesis eqn.2.62 or 2.59
    tau = self.params.correlationTime
    try:
        dE = self.params.dE
    except AttributeError:
        dE = self.params.dephaseConv / self.params.coherenceTime
    self.dE = dE

    # generate gaussian random numbers
    corrg = np.zeros_like(T)
    gauss = np.random.normal(size=len(T))
    gauss -= np.mean(gauss)  # make sure mean=0

    # if tau is too small then turn off correlation.
    if tau < 1:
        corrg = gauss
    else:
        # construct correlation with iterative algorithm
        f = np.exp(-1.*stepsize/tau)
        corrg[0] = gauss[0]
        for k in xrange(len(corrg)-1):
            corrg[k+1] = f*corrg[k] + np.sqrt(1-f**2)*gauss[k+1]
            corrg -= np.mean(corrg)

    self.decT = T
    self.dephaseV = corrg*dE   # scale to proper variance (corollary 1)
    #self.dephaseV = np.ones_like(corrg)*dE # for testing: remove randomization


  def plotcorrelation(self):
    """ plot the sequence and correlation function; for testing """
    acor = np.correlate(self.dephaseV, self.dephaseV, 'full')
    cor_fig = pl.figure()
    ax1 = cor_fig.add_subplot(211)
    ax1.plot(self.decT, self.dephaseV)
    ax1.set_title('dec.dephaseV')
    ax1.set_ylabel('phase')
    ax2 = cor_fig.add_subplot(212)
    ax2.plot(self.decT, acor[-len(self.dephaseV):])
    ax2.set_xlabel('time [us]')
    ax2.set_ylabel('correlation')
    self.acor = acor[-len(self.dephaseV):]

  def calcPoisson(self, T, stepsize, lifetime, ions=1):
    ''' generate poisson process for spontaneous decay and heating'''
    # poisson-distribution based on expected number of decays
    numjumps = ions*np.sum(np.random.poisson(float(stepsize) / lifetime, len(T)))
    jumptimes = np.random.uniform(0, T[-1], numjumps)

    self.decT = T
    # decT[nonzero(jumpVector)] is where jumps will occur
    jumpVector = np.zeros_like(T)
    for i in range(len(jumptimes)):
      jumpVector[T.searchsorted(jumptimes[i])] = 1

    return jumpVector

  def calcGaussian(self, T, stepsize, width, ions=1):
    ''' generate normally distributed noise for intensity fluct '''
    self.decT = T
    noiseVector = np.sqrt( 1+np.random.normal(scale=width, size=len(self.decT)) ) 
    return noiseVector

  def calcSpontaneousDecay(self, T, stepsize, lifetime, ions=1):
    ''' generate poisson process for spontaneous emission '''
    self.spontdecayV = self.calcPoisson(T, stepsize, lifetime, ions)

  def calcHeating(self, T, stepsize, rate):
    ''' generate poisson process for heating '''
    self.heatingV = self.calcPoisson(T, stepsize, rate)

  def calcIntensFluct(self, T, stepsize, width):
    ''' generate vector for intensity fluctuations '''
    self.intensfluctV = self.calcGaussian(T, stepsize, width)
    self.intensfluctV[self.intensfluctV<0] = 0

############################################################
# class parameters
############################################################
class parameters:
  """ all parameters for one simulation """

  def __init__(self, hspace, lab="UIBK"):
      self.hspace = hspace

      if lab=="UIBK":   # standard Innsbruck parameters
          self.omz = 2*pi*1.2  # secular frequency
          self.omc = 2*pi*0.023  # carrier freq in MHz
          self.omsb = 2*pi*0.023   # sideband freq in MHz
          self.pitime = pi / self.omc      # carrier pi-time
          self.omx = 2*pi*0.01  # corresponds to time for hiding/unhiding
          self.recoilangle_addr = 68  # addressed beam
          self.recoilangle_glob = 22 # global beam
          self.wavelength = 729e-9
          self.ionmass = 40 * 1.67e-27
      else:             # MIT parameters
          # TODO: read these from file
          self.omz = 2*pi*0.8
          self.pitime = 4.7
          self.omc = pi / self.pitime      # carrier frequency
          self.omsb = 1.1*pi / self.pitime  # drive sideband a little harder (no eta)
          self.recoilangle_addr = 0
          self.recoilangle_glob = 0
          self.wavelength = 674e-9
          self.ionmass = 88 * 1.67e-27

      # switching time between global/indiv pulses in us
      self.sameionswitchtime = 7
      self.diffionswitchtime = 46

      # spontaneous decay: lifetime of ca40
      self.lifetime = 1168000

      # intensity fluctuations
      self.intensityfluct = 0.02 # percentage fluctuations

      # coupling to spectator modes
      self.specmodecoupling = 0.02 # percentage to modify addressing matrix

      # initialization error
      # this is the probability of ending up in S'(m=+1/2)
      # because it's still 'bright' during the detection
      # we handle that like a non-coupling S state
      self.stateiniterr = 0.003

      # dephasing: UIBK parameters in us
      self.correlationTime = 333
      self.coherenceTime = 5000

      # dephasing-pulse: a constant phase offset on each pulse
      self.phaseOffset = 0.01

      # heating: 1 quanta per "params.heatingrate" us
      self.heatingrate = 141000 

      # addressing error in percentage (intensity)
      self.addressingerr = 0.035
      self.addressingerr_global = 0.02

      # detection settings in mus (like the lifetime)
      self.detection_time_ccd = 8000
      self.detection_time_pmt = 3000

      # hiding error - net population after each hide/unhide cycle
      self.hidingerr = 0.99
      self.hidingMeaserr = 0.005 # chance of ending up in D after each hide

      # pp servers - read from file
      self.server_dict = {}
      try:
          rootpath = PyTIQC.core.__path__
          f = open(rootpath[0]+'/server_dict', 'r')
          source = (f.readlines())
          f.close()
          for line in source:
              line.strip()
              if line[0] == '#': continue
              line = line.split(':')
              key = line[0]
              val = line[1].split('#')[0].strip()
              self.server_dict[key] = val
      except Exception, e:
          print "server_dict file parsing failed, setting server_dict to empty"
          print e
          self.server_dict = {}



      # for full hamiltonian, want to use lamb-dicke approx?
      self.LDapproximation = True

      # parameter for dephasing conversion - obtained by matching Ramsey decay
      self.dephaseConv = 2*pi*0.1551

      # parameters for AC stark pulse
      self.ACdelta = -2*pi*20 #MHz
      self.ACintenFac = 6.3 # faster omrabi to compensate for detuning
      self.omac = 2*pi/80.*2 # 2pi/80us, freq on the AC stark pulse
      # this one needs to be found experimentally for a given omac, here for 80us.
      self.ACcorr = 0.0927 # factor to add to HT for ACstark from S-P. 

      # parameter for MS pulse
      self.MSdelta = 2*pi*0.028 # detuning in MHz
      self.shortestMS = 2  # shortest MS pulse e.g. 16 for MS(pi/16) for QFT
      self.ODEtimestep = 0.01
      self.shape = 5 # pulseshaping time for MS gate
      self.doMSshapingCorr = True
      # these params need to be found for every detuning. Here they are for 28kHz.
      self.MScorr = {
          'duration': 1.1642,  #for shape=3: 1.0973
          'omrabi': 1.02672,
          'shortPulseFactorPos': {2: 1.000, 4: 1.0106, 8: 1.0158, 16: 1.01835},
          'shortPulseFactorNeg': {2: 1.010, 4: 1.0215, 8: 1.0271, 16: 1.02986}
          }
      # set of corrections to omc_ms due to hiding ions
      # row (1st index): number of ions participating in MS
      # column (2nd index): total number of ions
      self.doMShidecorr = False
      self.MShidecorr = np.array([ \
          [ -1, -1, -1, -1, -1, -1],
          [ -1, -1, -1, -1, -1, -1],
          [ -1, -1, -1, 1.0112, 1.0231, -1],
          [ -1, -1, -1, -1, 1.023, 1.0363]])

      # choose step size for storing solution
      self.stepsize = 1

      # default starting state
      self.y0 = self.hspace.y0
      self.y0_dict = None
 
      self.includelightshift = True

      # we're counting ions from right to left
      tmp_addressing = np.fliplr(np.eye(hspace.nuions))
      self.addressing = np.vstack([tmp_addressing, np.ones(hspace.nuions)]) # last line for global addressing

      self.detuning = 0
      self.dotimedep = False

      # choose an ODE solver for the MS gate
      self.solver = "ZVODE"

      # times are in us, freqs are in MHz. these factors help convert if needed
      # not used so far for anything yet.
      self.freqscale = 1e6
      self.timescale = 1e-6

      # print each pulse or progress bar displays for ODE solver
      self.printpulse = False
      self.progbar = True

      # pp settings. default is use all local cores
      self.ppcpus = 'autodetect' # local # of cores to use
      self.ppsecret = "tiqc_cluster1"
      self.ppservers = ['local']
      self.pplog = True  # do logging to file 'pp.log'

      # save data with timestamp of start time
      self.savedataname = time.strftime("%Y%m%d-%H%M%S") + '-data-autosave.shlv'
      self.savedata = False

      # option to not save intermediate data points
      self.saveallpoints = True
 
      self.calcEta()
      self.calcPulseParams()

  def initial_state(self, state="ground", nBar=0, qstate=None):
    """ define an initial state; wrapper to hspace's function """
    self.hspace.initial_state(state, nBar, qstate)
    self.y0 = self.hspace.y0

  def use_servers(self, servers):
      ''' select servers to use for pp '''
      if 'all' in servers or servers == 'all':
          self.ppcpus = 'autodetect' # local # of cores to use
          self.ppsecret = "tiqc_cluster1"
          self.ppservers = []
          for k in self.server_dict:
              self.ppservers.append(self.server_dict[k])
          self.ppservers = tuple(self.ppservers)
          if len(self.ppservers) > 0:
              self.ppcpus = 0
      else:
          for serv in servers:
              if serv != 'local':
                try:
                    self.ppservers.append(self.server_dict[serv])
                except KeyError:
                    print "server", serv, "not found in server_dict, skipping"
              else:
                  self.ppcpus = 'autodetect'
          self.ppservers = tuple(self.ppservers)
          if 'local' not in servers:
              self.ppcpus = 0

      if len(self.ppservers)==0 and self.ppcpus==0:
          print "no server defined, using local by default"
          self.ppcpus = 'autodetect'

  def calcEta(self):
    """ calculate the Lamb-Dicke parameter eta. It should have length = number of ions+1 (last one for global) """
    self.eta = np.zeros(self.hspace.nuions)
    self.eta = \
        np.ones(self.hspace.nuions) / \
        np.linalg.norm(np.ones(self.hspace.nuions)) * \
        self.lambdicke(self.wavelength, self.ionmass, self.omz/2/pi, self.recoilangle_addr, self.freqscale)
    self.eta = np.append(self.eta,  
        1/np.linalg.norm(np.ones(self.hspace.nuions)) * \
        self.lambdicke(self.wavelength, self.ionmass, self.omz/2/pi, self.recoilangle_glob, self.freqscale) )

  def setShortestMS(self, timefac):
    """ given an MS gate duration, calculate and reset self.shortestMS """
    for i in range(10): # assuming nothing shorter than pi/2**10 = pi/1024
        if timefac * 2**i == int(timefac * 2**i):
            self.shortestMS = self.shortestMS * 2**i 
            break
    print 'set shortestMS to', self.shortestMS
    self.calcPulseParams()

  def calcPulseParams(self):
    """ calculate rabi frequencies on various detunings """

    # lightshift correction
    if self.includelightshift:
        self.lsfactor = self.omc**2 / 2

    # coherence time - correct for # of ions
    #comment this out for now so it doesn't accidentally get called twice
    #self.coherenceTime = self.coherenceTime / float(self.hspace.nuions**2)

    # omc_ms is the carrier freq of the MS gate
    # set this according to the smallest MS gate we want to use
    self.omc_ms = abs(self.MSdelta) / self.eta[-1]
    frac = self.shortestMS / 2
    self.omc_ms = self.omc_ms / np.sqrt(frac)
    self.MSintenFac = self.omc_ms / self.omc

  def lambdicke(self, l, m, om, w, fs):
    """ calculate lamb-dicke factor for wavelength l, trap frequency om, recoil angle w, frequency scale fs """
    hq = constants.hbar
    y = 2*pi / l * np.sqrt(hq/(2*m*2*pi*om*fs))* abs(np.cos(w/360.*(2*pi)))
    return y

  def set_addressing(self):
    ''' modify the addressing matrix. Used by decoherence.'''
    tmp_addressing = np.zeros((self.hspace.nuions, self.hspace.nuions))
    for k in xrange(self.hspace.nuions):
      tmp_addressing[k,self.hspace.nuions-k-1] = 1
      # addressing error for addressed beam
      if self.hspace.nuions-k-1 -1 >= 0:
          tmp_addressing[k,self.hspace.nuions-k-1 -1] = self.addressingerr
      if self.hspace.nuions-k-1 +1 < self.hspace.nuions:
          tmp_addressing[k,self.hspace.nuions-k-1 +1] = self.addressingerr

    self.addressing = np.vstack([tmp_addressing, np.ones(self.hspace.nuions)]) # last line for global addressing
    N = self.hspace.nuions
    for k in xrange(N):
        middleion = [int(N)/2, int(N)/2 - (1 if N%2==0 else 0)]
        self.addressing[-1, k] = (1-self.addressingerr_global)**min([abs(k-middleion[0]), abs(k-middleion[1])])

  def use_rho0(self, rho0):
    ''' given a densmat as a starting state, decompose it into eigenvectors '''
    [e, v] = splg.eig(rho0)
    e = np.real(e)
    self.y0_dict = {}
    # first calculate the normalization factor
    etotal = 0
    for i in range(len(e)):
        if e[i] > 0.01:
            etotal += e[i]
    # then put the ones with eigenval > 0.01 into a dict
    for i in range(len(e)):
        if e[i] > 0.01:
            y0 = np.zeros_like(self.y0)
            y0[::self.hspace.maxphonons+1] = np.conj(v[:,i])
            self.y0_dict[e[i]/etotal] = y0

############################################################
# class pulse
############################################################
class pulse:
  """ base class for actual laser pulses (Rcar, Rblue, Rac, RMS).
    (Do not instantiate this class directly.)
    | Note: put pi's in the function call! helps to remember which parameter is which.
  """

  def __init__(self, params, theta, phi, ion, use_ideal=False):
    self.theta = theta
    self.phase = phi
    self.starttime = -1
    self.endtime = -1
    self.ion = ion # saved to allow for changes later
    self.U = 1 #unitary for the pulse

    # option to do timedep HT (e.g. for MS) regardless of param.dotimedepending
    self.dotimedepPulse = False
    self.dobichro = False

    # select ion
    self.targetion = np.copy(params.addressing[ion, :])

    # option to use ideal unitary (from ideal/gates.py) to do evolution
    self.use_ideal = use_ideal

    # a list to keep all unitaries generated in every simulation
    self.UtrAll = []

  def calculateIdealUnitary(self, params, hiddenions=[]):
    """ calculate the ideal unitary for this pulse """
    nuions = params.hspace.nuions
    op_ida = params.hspace.operator_dict['id_a']
    import PyTIQC.core.gates as U # not sure why I need this now ... rev 255

    if self.type == "R":
        if self.ion==-1:
            self.Uidtr = U.Rg(self.theta, self.phase, nuions, hiddenions)
        else:
            self.Uidtr = U.Ri(self.theta, self.phase, nuions, self.ion, hiddenions)
    elif self.type == "Z":
        self.Uidtr = U.Uz(self.theta, nuions, self.ion, hiddenions)
    elif self.type == "M":
        self.Uidtr = U.MS(np.sign(params.MSdelta)*self.theta, self.phase, nuions, hiddenions)
    elif self.type == "D" or self.type == "H" or self.type == "I":
        self.Uidtr = U.G(U.I, nuions)/nuions

    elif self.type == "Rb":
        # no phonons defined in core.gates, so directly use the Hamiltonian
        h = qm.Hamilton()
        HT, lqc = h.Hamiltonian(self, params)
        tlen = self.theta / self.omrabi / params.eta[self.ion]
        Ugate = splg.expm2(-1j * tlen * HT)
        Ulqc1 = np.diag(np.exp(-1j * 0 * lqc))
        Ulqc2 = np.diag(np.exp(-1j * (-tlen) * lqc))
        U = np.mat(Ulqc2) * np.mat(Ugate) * np.mat(Ulqc1)        
        #U = splg.expm2(-1j * tlen * HT)
        self.Uid = np.asarray(U)
        self.Uidtr = np.diag(np.ones(len(np.diag(self.Uid))))
    else:
        raise ValueError("Ideal unitary for pulse not defined")

    if self.type != "Rb":
        self.Uid = np.kron(self.Uidtr, op_ida)


  def traceM(self, hspace):
    """ calculate the unitary, ignore the motional qubit """
    if isinstance(self.U, int):
        self.Utr = self.U
    else:
      nuions = hspace.nuions
      levels = hspace.levels
      phonons = hspace.maxphonons+1
      self.Utr = np.zeros([levels**nuions, levels**nuions], np.complex128)
      for i in range(levels**nuions):
        for j in range(levels**nuions):
          qi = i*phonons
          qj = j*phonons
          self.Utr[i, j] = self.U[qi, qj]

  def maketimedep(self, pulseshape=None, shapingCorr=None):
    """ to be called by simulateevolution """
    self.omrabi_t = lambda t: self.omrabi

  def __str__(self):
      return self.type + '(' + str(np.around(self.theta,3)) + ', ' + str(np.around(self.phase,3)) + ', ' + str(self.ion) + ')'


########
# Rcar
########
class Rcar(pulse):
  """ carrier pulse.
      Use detuning correction for ignorelightshift only if
      detuning > 2*omc, so the taylor expansion is valid. """
  def __init__(self, params, theta, phi, ion=-1, use_ideal=False):
    pulse.__init__(self, params, theta, phi, ion, use_ideal)
    self.type = 'R'
    self.detuning = params.detuning
    if not params.includelightshift and self.detuning > 2*params.omc:
      self.detuning = self.detuning - params.lsfactor / self.detuning

    self.duration = theta / params.omc
    self.omrabi = params.omc

    self.calculateIdealUnitary(params)

########
# Rblue
########
class Rblue(pulse):
  """ pulse on the blue sideband. Always do detuning correction if ignorelightshift is on. """
  def __init__(self, params, theta, phi, ion=-1, use_ideal=False):
    pulse.__init__(self, params, theta, phi, ion, use_ideal)
    self.type = 'Rb'
    self.detuning = params.detuning + params.omz
    if not params.includelightshift:
      self.detuning = self.detuning - params.lsfactor / self.detuning

    self.duration = theta / params.omsb / params.eta[ion]
    self.omrabi = params.omsb
    self.phase = self.phase - pi/2 # to match hartmut's

    if self.use_ideal:
        print "warning: ideal U for Rblue not defined"

########
# Rac
########
class Rac(pulse):
  """ same as Rcar, with detuning = params.ACdelta (20 MHz) and intensity adjusted """
  def __init__(self, params, theta, phi, ion=0, use_ideal=False):
    pulse.__init__(self, params, theta, phi, ion, use_ideal)
    self.type = 'Z'
    self.detuning = params.detuning + params.ACdelta
    if not params.includelightshift:
      self.detuning = self.detuning - params.lsfactor / self.detuning

    self.duration = abs(theta) / params.omac
    self.omrabi = params.omc * params.ACintenFac

    self.calculateIdealUnitary(params)

########
# RMS
########
class RMS(pulse):
  """ MS gate: bichromatic pulse with timedep always on """
  def __init__(self, params, theta, phi, ion=-1, use_ideal=False):
    # optimize for one gate and adjust intensity accordingly
    self.type = 'M'

    self.MSdelta = params.MSdelta if theta > 0 else -params.MSdelta
    self.timefactor = abs(theta) / (pi/params.shortestMS)
    if abs(self.timefactor) < 1:
        params.setShortestMS(self.timefactor)
        self.timefactor = theta / (pi/params.shortestMS)

    if np.sign(self.MSdelta) > 0:
        shortPulseFactor = params.MScorr['shortPulseFactorPos']
    else:
        shortPulseFactor = params.MScorr['shortPulseFactorNeg']
    if params.shortestMS in shortPulseFactor:
        self.spf = shortPulseFactor[params.shortestMS]
    else:
        self.spf = 1

    pulse.__init__(self, params, theta, phi, ion, use_ideal)

    # intermediate variables for MS params
    car_det = params.omz + self.MSdelta
    omrabi = params.omc_ms

    # here are the parameters that exist for every pulse
    self.detuning = 0  # center-line detuning?
    self.detuning_b = car_det
    self.detuning_r = -car_det
    self.omrabi = omrabi
    self.omrabi_r = omrabi
    self.omrabi_b = omrabi

    # phase of the beat signal
    self.phase_rb = 0
    # phase of the light field. should define x/y
    self.phase_light = pi-phi

    # bichromatic part
    self.dotimedepPulse = True
    self.dobichro = True

    # the factor of 1.04 in the pulse duration has been figured
    # out via simulation and _only_ works for shape=3 and shortestMS up to 16
    # other shape times will likely need a different factor
    self.duration_fac = params.MScorr['duration']
    self.duration_base = 2*pi / abs(params.MSdelta)
    if params.shape == None:
      self.duration =  self.timefactor * self.duration_base
      #  self.duration = theta
    else:
      self.duration = self.duration_fac * self.timefactor * self.duration_base

    self.omrabi_corr = params.MScorr['omrabi']

    self.calculateIdealUnitary(params)

  def maketimedep(self, pulseshape, shapingCorr):
    """ to be called by simulateevolution, and also calculate pulse shaping """
    # with optional pulse shaping
    # not sure where this factor comes from yet ...
    corr = self.omrabi_corr
    if pulseshape == None:

      self.omrabi_rt = lambda t: corr * self.omrabi_r / 2
      self.omrabi_bt = lambda t: corr * self.omrabi_b / 2
    else:
      if shapingCorr:
        self.omrabi_rt = lambda t: corr * self.omrabi_r / 2 * self.BlackmanShape(t, self.starttime, self.endtime, pulseshape, self.duration_fac*self.duration_base) / self.spf
        self.omrabi_bt = lambda t: corr * self.omrabi_b / 2 * self.BlackmanShape(t, self.starttime, self.endtime, pulseshape, self.duration_fac*self.duration_base) / self.spf
      else:
        self.omrabi_rt = lambda t: corr * self.omrabi_r / 2 * self.BlackmanShape(t, self.starttime, self.endtime, pulseshape) / self.spf
        self.omrabi_bt = lambda t: corr * self.omrabi_b / 2 * self.BlackmanShape(t, self.starttime, self.endtime, pulseshape) / self.spf


  def BlackmanShape(self, t, starttime, endtime, ramptime=3, rampevery=None):
    ''' generate the waveform for pulse shaping.
        returns the field? intensity? factor at time t '''
    # this is a very ugly cheat but it works. It sets the duration 
    # of a multiple-of-shortestMS pulse the same as repeating the shortestMS pulse multiple times. 
    # then only a small set (1/shortestMS) of power correction factors need to be defined.

    # make the interval markers
    if rampevery:
        tposmarkers = np.arange(starttime, endtime, rampevery)
        tposramp = tposmarkers + ramptime
        tnegmarkers = np.arange(starttime, endtime, rampevery) + rampevery
        tnegramp = tnegmarkers - ramptime

    # check if t is in the positive-ramp region
    if rampevery:
        tt = None
        for i in range(len(tposmarkers)):
            if t > tposmarkers[i] and t < tposramp[i]:
                tt = t - tposmarkers[i]

        for i in range(len(tnegmarkers)):
            if t > tnegramp[i] and t < tnegmarkers[i]:
                tt = tnegmarkers[i] - t

        if tt:
            return 1./2 * (0.84 - np.cos(1.*tt/ramptime * np.pi) + 0.16 * np.cos(2.*tt/ramptime * np.pi))
        else:
            return 1.
    else:
        if t > starttime+ramptime and t < (endtime - ramptime):
            return 1.
        else:
            if starttime < t and t < starttime+ramptime:
                tt = t - starttime
            else:
                tt = endtime - t
            return 1./2 * (0.84 - np.cos(1.*tt/ramptime * np.pi) + 0.16 * np.cos(2.*tt/ramptime * np.pi))


########
# Hide
########
class Hide(pulse):
  """ pulse to hide/unhide ions. like Delay, except modifies the addressing matrix """
  def __init__(self, params, ion, hide, use_ideal=False):
    pulse.__init__(self, params, pi, 0, ion, use_ideal)
    self.type = 'H'
    self.hide = hide # set hide=0 to unhide
    self.duration = pi / params.omx

    self.phase = 0
    self.detuning = 0
    self.omrabi = 0
    self.targetion = np.zeros(params.hspace.nuions)
    self.dotimedepPulse = False

    self.calculateIdealUnitary(params)

########
# Delay
########
class Delay(pulse):
  """ do nothing. used to keep track of times. """
  def __init__(self, params, tdelay):
    pulse.__init__(self,params, 0, 0, ion=-1, use_ideal=False)
    self.type = 'D'
    self.phase = 0
    self.duration = tdelay
    self.detuning = 0
    self.omrabi = 0
    self.targetion = np.zeros(params.hspace.nuions)
    self.dotimedepPulse = False

    self.calculateIdealUnitary(params)

##########
# MeasInit
##########
class MeasInit(pulse):
  """ a special instruction to measure and initialize one qubit.
  | experiment: hide other qubits, measure (scatter/detect), apply sigma beam (to re-init to S), raman cool (re-init phonons), unhide other qubits.
  | simulation (ideal): trace over other qubits, store pop of qubit in register, project qubit to S. """
  def __init__(self, params, ion, incl_hidingerr=True):
      pulse.__init__(self, params, 0, 0, ion, use_ideal=False)
      self.type = "I"
      self.hspace = params.hspace

      self.duration = 600 # time in us for measurement & recooling
      self.phase = 0
      self.detuning = 0
      self.omrabi = 0
      self.dotimedepPulse = False

      self.incl_hidingerr = incl_hidingerr

      self.calculateIdealUnitary(params)

  def measure(self, y):
      ''' returns the traced population of the desired qubit'''
      reg = np.zeros(2)
      for (ind, pop) in enumerate(abs(y**2)):
          state = qm.indexToState(ind, self.hspace)
          if state[1][self.hspace.nuions-self.ion-1] == 1: # S
              reg[0] += pop
          else: # = 0, D
              reg[1] += pop

      return reg

  def Uinit(self):
      ''' returns projection matrix to re-initialize qubit to 0 '''
      projs = self.hspace.operator_dict['projs']
      projmtx = npml.kron( U.R(projs, self.hspace.nuions, self.ion), self.hspace.operator_dict['id_a'] )
      return projmtx
      

############################################################
# class PulseSequences
############################################################
class PulseSequence():
    """ a container class for a pulse sequence. """

    def __init__(self, seq, addswitchtime=False, params=None):
        """ initialize the sequence """
        if type(seq) != list:
            raise Exception("Pulseseq should be a list of pulses")

        self.seq = []
        for item in seq:
            if hasattr(item, 'seq'):
                for pulse in item.seq:
                    self.seq.append(copy.deepcopy(pulse))
            else:
                self.seq.append(item)

        self.seqqc = self.seq # this list optionally has delays inserted, to be used in qctools
        self.makepulsesequence()

        if addswitchtime:
            if params==None:
                raise Exception("must specify params to add switch times")
            self.addDelays(params)

    def makepulsesequence(self):
        """ set the times properly in a pulse sequence (starttime, endtime).
        This avoids declaring time as a global variable.
        """
        time = 0
        for pulse in self.seq:
            pulse.starttime = time
            time = time + pulse.duration
            pulse.endtime = time
        self.totaltime = time

    def __len__(self):
        """ get length of the sequence """
        return len(self.seq)

    def __getitem__(self, i):
        """ indexing operator: get the ith pulse """
        return self.seq[i]

    def __setitem__(self, i, pulse):
        """ set the ith pulse """
        self.seq[i] = pulse
        self.makepulsesequence()

    def __contains__(self, b):
        """ 'in' operator: outcome of test b in self.seq."""
        return b in self.seq

    def prepend(self, prepend_pulses):
        """ prepend the sequence with another sequence """
        for k in xrange(len(prepend_pulses)):
            self.seq.insert(0, prepend_pulses[len(prepend_pulses)-1 -k])
        self.makepulsesequence()
        self.seqqc = self.seq

    def append(self, append_pulses):
        """ append the sequence with another sequence """
        for k in append_pulses:
            self.seq.append(k)
        self.makepulsesequence()
        self.seqqc = self.seq

    def extend(self, pulseseq):
        """ like append, but input is a pulse sequence object """
        self.seq.extend(pulseseq.seq)
        self.makepulsesequence()
        self.seqqc = self.seq

    def __add__(self, other):
        ''' concatenation operator '''
        self.append(other)
        return self

    def changeions(self, params, ionindices):
        ''' for a given pulse sequence, input is a list or tuple 
        indicating which ions to apply the Z gates to.
        NOTE: this assumes the default sequence uses ions in sequence 0,1,2... '''
        for pulse in self.seq:
            if pulse.type == "Z":
                pulse.ion = ionindices[pulse.ion]
                pulse.targetion = np.copy(params.addressing[pulse.ion, :])
            elif pulse.type == "H":
                pulse.ion = ionindices[pulse.ion]

    def addDelays(self, params):
        ''' add appropriate delays to pulseseq to account for ion switching'''

        if len(self) == 0:
            self.seqqc = []
            return 

        # first make a new sequence and insert blanks where Delays go
        self.seqqc = []
        for pulse in self.seq:
            self.seqqc.append(pulse)
            self.seqqc.append(None)

        # insert the proper delays
        k = 0
        while k+2 < len(self.seqqc):
            if self.seqqc[k].ion == self.seqqc[k+2].ion:
                self.seqqc[k+1] = Delay(params, params.sameionswitchtime)
            else:
                self.seqqc[k+1] = Delay(params, params.diffionswitchtime)

            # inherit use_ideal from previous pulse
            self.seqqc[k+1].use_ideal = self.seqqc[k].use_ideal
            k+=2
        self.seqqc.pop()  # remove the last None

        # set the times properly
        time = 0
        for pulse in self.seqqc:
            pulse.starttime = time
            time = time + pulse.duration
            pulse.endtime = time    
        self.totaltime = time

        # now go through the original seq and change the times to match.
        for [i, pulse] in enumerate(self.seq):
            pulse.starttime = self.seqqc[2*i].starttime
            pulse.endtime = self.seqqc[2*i].endtime

    def getdigitalrepresentation(self):
        ''' return a time vector with 1's where a pulse is on and 0's where not, to only us precision. '''
        t = np.zeros(self.totaltime)
        for pulse in self.seq:
            t[int(pulse.starttime) : int(pulse.endtime)] = 1
        return t

    def __str__(self):
        ''' print the entire sequence in format similar to input. '''
        outstr = ''
        for n, pulse in enumerate(self.seq):
            outstr+= str(n) + ': ' + pulse.type + '(' + str(np.around(pulse.theta,3)) + ', ' + str(np.around(pulse.phase,3)) + ', ' + str(pulse.ion) + ')'
            outstr+= '\n'
        return outstr

    def counttypes(self):
        ''' count the number of types of each pulse '''
        types_dict = {}
        for pulse in self.seq:
            if pulse.type not in types_dict.keys():
                types_dict[pulse.type] = 1
            else:
                types_dict[pulse.type] += 1

        shortform = ''
        for [ptype, count] in types_dict.items():
            print count, ptype, "pulses"
            shortform += str(count) + ptype + ","
        print shortform

    def totalUnitary(self):
        ''' return the ideal unitary of the whole sequence '''
        uni = 1
        for pulse in self.seq:
            uni = np.dot(pulse.Uidtr, uni)

        return uni

############################################################
# function plotBar3D
############################################################
def plotBar3D(rho, displ=1):
  """ generic function to make a 3d plot of density matrices """
  j = np.seterr(divide='ignore')

  dim = len(rho)

  fig = pl.figure(displ)

  xpos, ypos = np.meshgrid(np.arange(dim), np.arange(dim))
  xpos = xpos.flatten()
  ypos = ypos.flatten()
  zpos = np.zeros(len(xpos.flatten()))
  dx = 0.5*np.ones_like(xpos)
  dy = dx.copy()
  
  rect1 = fig.add_subplot(1,3,1).get_position() 
  ax1 = Axes3D(fig, rect1)
  ax1.bar3d(xpos, ypos, zpos, dx, dy, abs(rho).flatten())
  ax1.set_zlim3d([0,1])

  rect2 = fig.add_subplot(1,3,2).get_position() 
  ax2 = Axes3D(fig, rect2)
  realrho = np.real(rho)
  ax2.bar3d(xpos, ypos, zpos, dx, dy, realrho.flatten())
  ax2.set_zlim3d([-0.5,0.5])

  rect3 = fig.add_subplot(1,3,3).get_position() 
  ax3 = Axes3D(fig, rect3)
  imagrho = np.imag(rho)
  ax3.bar3d(xpos, ypos, zpos, dx, dy, imagrho.flatten())
  ax3.set_zlim3d([-0.5,0.5])

  pl.show()

############################################################
# function loaddata
# parallel to data.save()
############################################################
def loaddata(filename):
    """ parallel to data.save() """
    infile = open(filename, 'rb')
    data = pickle.load(infile)
    infile.close()

    if data.averaged == False:
        data.mean(data.numruns)

    return data

############################################################
# class database
############################################################
class database:
  """ tools to manipulate data. class members are various views to the data:

  | T, Y = time, full states
  | YR = all state populations
  | Yend = final state
  | Pend = final populations
  | YRend = final state populations traced over phonons
  | RhoEnd = final density matrix traced over phonons
  | YtrI = all state populations traced over ions
  | YtrN = all state populations traced over phonons
  | TP,YP = time and states at end of each pulse (or next timestep)
  | YRP = populations at end of each pulse (or next timestep)
  | YRPN = YRP traced over motional states
  | RhoPN = density matrices traced over motional states at each pulse
  | RhoPNAll = All density matrices for all MC instances. Everything else can be extracted from RhoPNAll.
  """

  def __init__(self, T, Y, hspace, pulseseq=None, register=[], statesvalid=True):
    self.T = T
    self.Y = Y
    self.hspace = hspace
    self.statesvalid = statesvalid # if averaged then set statesvalid to false

    self.phonons = hspace.maxphonons+1
    self.nuions = hspace.nuions
    self.levels = hspace.levels

    self.numruns = 1
    self.numstates = 1
    self.averaged = True # to keep track of averaging/saving
    if not statesvalid:
      self.YR = np.real(Y)
      self.Pend = self.YR[-1]
      self.Yend = None
    else:
      self.YR = abs(Y)**2
      self.Yend = Y[-1]
      self.Pend = abs(self.Yend)**2

    if pulseseq != None:
        self.pulseseq = pulseseq
        self.statesatpulse()
        # save all densmats at all pulses for calculating errors later
        self.RhoPNAll = []
        self.RhoPNAll.append(self.RhoPN)
        # state and densmat at the end
        self.YRend = self.YRPN[-1]
        self.RhoEnd = self.RhoPN[-1]
        # classical registers
        self.registerAll = [] # all classical registers
        self.registerAll.append(register)
        self.register = register

        for pulse in self.pulseseq:
            pulse.traceM(hspace)
            try:
                pulse.UtrAll = []
                pulse.UtrAll.append(pulse.Utr)
            except Exception, e:
                print "Exception in database creation:", e

  def __iadd__(self, other): 
    self.YR = self.YR + other.YR
    # add the states and densmats at each pulse
    self.YRPN = self.YRPN + other.YRPN
    self.RhoPN = self.RhoPN + other.RhoPN
    self.YRend = self.YRend + other.YRend
    self.RhoEnd = self.RhoEnd + other.RhoEnd
    self.RhoPNAll.append(other.RhoPN)
    self.statesvalid = False
    self.numruns += 1
    self.averaged = False

    if self.pulseseq and other.pulseseq:
        for i in range(len(self.pulseseq)):
            self.pulseseq[i].UtrAll.append(other.pulseseq[i].Utr)

    # add and store the classical registers
    self.registerAll.append(other.register)
    for n in range(len(self.register)):
        self.register[n] = self.register[n] + other.register[n]

    return self

  def mean(self, k):
    """ Calculate average of populations and density matrices """
    self.YR = self.YR/k
    self.YRend = self.YRend/k
    self.RhoEnd = self.RhoEnd/k
    self.YRPN = self.YRPN/k
    self.RhoPN = self.RhoPN/k

    self.statesvalid = False
    self.RhoPNAll = np.array(self.RhoPNAll)
    self.averaged = True

    if self.pulseseq:
        for pulse in self.pulseseq:
            pulse.UtrAll = np.array(pulse.UtrAll)

    for n in range(len(self.register)):
        self.register[n] = self.register[n]/k

  def times(self, k):  # shouldn't __mul__ work here?
    self.YRPN = self.YRPN*k
    self.RhoPN = self.RhoPN*k
    self.YRend = self.YRend*k
    self.RhoEnd = self.RhoEnd*k
    self.RhoPNAll = np.array(self.RhoPNAll)*k
    self.registerAll = np.array(self.registerAll)*k
    for n in range(len(self.register)):
        self.register[n] = self.register[n]*k

    return self

  def addstate(self, other): # only slightly different from __iadd__
    self.YRPN = self.YRPN + other.YRPN
    self.RhoPN = self.RhoPN + other.RhoPN
    self.YRend = self.YRend + other.YRend
    self.RhoEnd = self.RhoEnd + other.RhoEnd
    self.RhoPNAll = self.RhoPNAll + other.RhoPNAll
    self.statesvalid = False
    self.numstates += 1
    self.registerAll = self.registerAll + other.registerAll
    for n in range(len(self.register)):
        self.register[n] = self.register[n] + other.register[n]

    return self

  def save(self,filename=None):
    """ save the object to a file; default filename is timestamped """
    if not filename:
      filename = time.strftime("%Y%m%d-%H%M") + '-data.pk'
    output = open(filename, 'wb')
    pickle.dump(self, output)
    output.close()

  def statesatpulse(self, displ=0):
    ''' save the states at the end of pulse times '''
    if self.pulseseq == None:
        print "Error: self.pulseseq not defined"
        return

    self.TP = np.zeros(len(self.pulseseq)+1)
    for i in range(len(self.pulseseq)):
      self.TP[i+1] = self.pulseseq[i].endtime
    indP = self.T.searchsorted(self.TP)
    if self.statesvalid:
      self.YP = self.Y[indP, :]
    self.YRP = self.YR[indP, :]

    # trace over motional states
    self.YRPN = np.zeros((np.size(indP), self.levels**self.nuions))
    if self.statesvalid:
      RhoP = np.zeros((np.size(indP), np.shape(self.Y)[1], np.shape(self.Y)[1]),np.complex128)
      for i in range(len(indP)):
        RhoP[i,:,:] = np.outer(np.conj(self.Y[indP[i],:]), self.Y[indP[i],:])
        self.RhoPN = np.zeros((np.size(indP), self.levels**self.nuions, self.levels**self.nuions), np.complex128)

    for i in range(self.levels**self.nuions):
      self.YRPN[:,i] = np.sum(self.YRP[:,i*self.phonons:(i+1)*self.phonons],1)
      if self.statesvalid:
        for j in range(self.levels**self.nuions):
          self.RhoPN[:,i,j] = np.trace( RhoP[:, \
              i*self.phonons:(i+1)*self.phonons,
              j*self.phonons:(j+1)*self.phonons], \
              axis1 = 1, axis2 = 2)

    # make plot
    if not(displ): return
    else:
      pulsestate_fig = pl.figure(displ)
      pulsestate_fig.clf()
      ax = pulsestate_fig.add_subplot(111)
      ax.plot(self.YRPN)
      pl.ylim(0,1)
      pl.xlabel('Pulse number')
      pl.ylabel('State population')
      pulsestate_fig.show()

  def endpopulation(self, printstates=True):
    """ end state traced over motional qubit """

    if printstates:
      for j in range(len(self.Pend)):
        if self.statesvalid:
          print "|%s>: %1.3f, %2.3f + %1.3fi" \
            %(qm.indexToState(j, self.hspace)[2], self.Pend[j], np.real(self.Yend[j]), np.imag(self.Yend[j]))
        else:
          print "|%s>: %1.3f" \
            %(qm.indexToState(j, self.hspace)[2], self.Pend[j])

  def tracedpopulation(self, displ):
    """ trace over the motional qubit. plot D state pop for each ion. """
    # trace over motional states
    self.YtrN = np.zeros((np.size(self.T), self.levels**self.nuions))
    for i in range(self.levels**self.nuions):
      self.YtrN[:,i] = np.sum(self.YR[:,i*self.phonons:(i+1)*self.phonons],1)

    # trace over ion
    # use bitmasks to find states
    self.YtrI = np.zeros((np.size(self.T), self.nuions))
    for i in range(self.nuions):
      ionind = 1 << i
      for j in range(self.hspace.dimensions):
        s = qm.indexToState(j, self.hspace)[2].split(',')[0]
        s = s.replace('S','1').replace('D','0').replace('A','2')
        sid = int(s,self.levels)
        if sid & ionind :
          self.YtrI[:,i] = self.YtrI[:,i] + self.YR[:,j]

    # make plot
    if not(displ): return
    else:
      traced_fig = pl.figure(displ)
      traced_fig.clf()
      ax = traced_fig.add_subplot(111)
      pl.hold(True)
      for i in range(self.nuions):
        ax.plot(self.T, 1-self.YtrI[:,i], label="ion %d"%(self.nuions-i-1)) #count ions in reverse
      pl.hold(False)

      pl.ylim(0,1)
      ax.legend()
      pl.xlabel('Time [us]')
      pl.ylabel('D state population')
      traced_fig.show()

  def displaypopulation(self, displ, plotlegend=False):
    """ display final populations in text form and
        make plots of population for each ion+phonon state"""
    # calculate abs populations and end states
    self.endpopulation(printstates=False)
    # print end results to terminal
    print "Final population:"
    for j in range(len(self.Pend)):
      if self.Pend[j]>=0.001: # only print non-zero populations
        if self.statesvalid:
          print "|%s>: %1.3f, %2.3f + %1.3fi" \
            %(qm.indexToState(j, self.hspace)[2], self.Pend[j], np.real(self.Yend[j]), np.imag(self.Yend[j]))
        else:
          print "|%s>: %1.3f" \
            %(qm.indexToState(j, self.hspace)[2], self.Pend[j])

    # make plot
    if not(displ): return
    else:
      population_fig = pl.figure(displ)
      population_fig.clf()
      for j in range(self.nuions*self.levels):
        ax = population_fig.add_subplot(self.nuions*self.levels, 1, j+1)
        ax.hold(True)
        for k in range(self.phonons):
          if self.hspace.visible[j*self.phonons+k]:
            ax.plot(self.T, self.YR[:,j*self.phonons+k], \
                 label=qm.indexToState(j*self.phonons+k,self.hspace)[2])
            pl.ylim(0,1)
            pl.ylabel('population')

        ax.hold(False)
        if plotlegend:
          ax.legend()

      pl.xlabel('Time [us]')
      population_fig.show()

  def displayPMTpopulations(self, displ, plotlegend = False):
    """ make plot of populations you would observe with a PMT"""
    try: self.YtrN
    except AttributeError:
      self.tracedpopulation(0)
      # now we have self.YtrN ... already traced
      # over phonons

    # 0 to N excitations ... so it's N+1 excitations we can
    # detect on the PMT
    self.P_PMT = np.zeros((len(self.T), self.nuions+1))

    # we already traced out the phonons. make a fake hspace
    # to use standard function for obtaining the correct
    # labeling/entries
    fake_hspace = copy.deepcopy(self.hspace)
    fake_hspace.maxphonons = 0

    for k in xrange(2**self.nuions):
      # get number of excitations, add to the corresponding column
      # in the PMT populations
      NumberOfExcitations = qm.indexToExcitations(k, fake_hspace)
      # keep in mind that excitations = D = 0
      self.P_PMT[:,self.nuions-NumberOfExcitations] += np.abs(self.YtrN[:,k])

    self.P_PMT_end = self.P_PMT[-1,:]

    if self.pulseseq != None:
        self.TP = np.zeros(len(self.pulseseq)+1)
        for i in range(len(self.pulseseq)):
          self.TP[i+1] = self.pulseseq[i].endtime
        indP = self.T.searchsorted(self.TP)
        self.P_PMT_P = self.P_PMT[indP, :]

    if not(displ): return
    else:
      pmtpopulation_fig = pl.figure(displ)
      pmtpopulation_fig.clf()

      ax = pmtpopulation_fig.add_subplot(111)
      ax.hold(True)

      for k in xrange(self.nuions+1):
        ax.plot(self.T, self.P_PMT[:,k], label = str(k))
        pl.ylim(0,1)

      ax.hold(False)
      if plotlegend:
        ax.legend()

      pl.ylabel('PMT excitations')
      pl.xlabel('Time [us]')
      pmtpopulation_fig.show()

  def displayPhaseSpace(self, displ):
    """
    plot timeevolution in phasespace
    some math tricks need to be done to get the correct result
    """

    self.x_p = np.zeros((len(self.T),2))

    a_full = npml.kron(np.diag(np.ones(self.hspace.levels**self.hspace.nuions)), self.hspace.operator_dict['a'])
    a_dagger_full = npml.kron(np.diag(np.ones(self.hspace.levels**self.hspace.nuions)), self.hspace.operator_dict['a_dag'])


    x_op = a_full+a_dagger_full
    p_op = 1j*(a_full-a_dagger_full)

    for k in xrange(len(self.T)):
      Y_tmp = self.Y[k,:]

      x = np.real(np.dot(np.dot(Y_tmp.conjugate().transpose(), x_op ), Y_tmp))
      p = np.real(np.dot(np.dot(Y_tmp.conjugate().transpose(), p_op ), Y_tmp))

      self.x_p[k,0] = x
      self.x_p[k,1] = p

    if not(displ): return
    else:
      phasespace_fig = pl.figure(displ)
      phasespace_fig.clf()

      ax = phasespace_fig.add_subplot(111)

      ax.plot(self.x_p[:,0],self.x_p[:,1])

      pl.ylabel('momentum')
      pl.xlabel('location')
      phasespace_fig.show()



############################################################
# TODO: put the following functions (for spont.decay errors during measurement)
# into a proper class, or absorb into qmtools, or ...?

def DecayedPopulations_CCD(populations, params):
    """
    given a vector of populations, return populations that would
    correspond to decay during the ccd measurement time (calculated
    from the measurement time provided in params)
    """

    populations_decayed = np.zeros(len(populations))
    # make sure everything is float
    for k in xrange(len(populations)):
      populations_decayed[k] = populations[k]

    # standard labeling starting with D..D to S...S
    NumberOfIons = int(np.log2(len(populations)))

    # multiplier used later on to easily derive the
    # population index
    multiplier = 2**(np.arange(NumberOfIons)[::-1])

    fake_hspace = hspace(NumberOfIons)

    # first calculate the decayed populations,
    # we'll then sample from that
    decay_prob = 1 - np.exp(-params.detection_time_ccd*1. / params.lifetime)

    for k in xrange(len(populations)):
      # for each population
      StateArray = qm.indexToState(k, fake_hspace)[1] # with D = 0
#      print StateArray
#      print type(StateArray)
      NumberOfExcitations = qm.indexToExcitations(k, fake_hspace)
      # directly tells us the max number of decays possible
      for l in np.arange(1,NumberOfExcitations+1):
          decayed_state_list = generate_decay_list([StateArray], decays = l)
          for decayed_state in decayed_state_list:
              ind = int(np.sum(decayed_state * multiplier))
#              print ind
              shift_pop = decay_prob**l * populations[k]
              populations_decayed[k] -= shift_pop
              populations_decayed[ind] += shift_pop


#    print populations_decayed
    return populations_decayed

def generate_decay_list(statearray_list, decays=1):
    # statearray = [1,0] for D, S
    # 0 = D
    # 1 = S
    # generate all possible decays in a table like
    #    0 = input
    #    1 = all states with one decay
    #    2 = all states with two decays
    # etc

    if type(statearray_list) == np.ndarray:
        statearray_list = statearray_list.tolist()

    statearrays_new = []
    for statearray in statearray_list:
        if type(statearray) == np.ndarray:
            statearray = statearray.tolist()
        for k in xrange(len(statearray)):
            if statearray[k] == 0:
                statearray_decayed = copy.copy(statearray)
                statearray_decayed[k] = 1
                statearrays_new.append(statearray_decayed)
    # remove dublicates
    out = []
    for element in statearrays_new:
        if element not in out:
            out.append(element)

    if decays>1:
        out = generate_decay_list(out, decays = decays-1)

    return out


def DecayedExcitations_PMT(excitations, params):
    """
    given a vector of excitations, return excitations that would
    correspond to decay during the pmt measurement time (calculated
    from the measurement time provided in params)
    """
    excitations_decayed = np.zeros(len(excitations))
    # make sure everything is float
    for k in xrange(len(excitations)):
      excitations_decayed[k] = excitations[k]

    # standard labeling starting with D..D to S...S
    NumberOfIons = int(len(excitations)-1)

    # first calculate the decayed populations,
    # we'll then sample from that
    decay_prob = 1 - np.exp(-params.detection_time_pmt*1. / params.lifetime)

    # # figure out all permutations of decays
    # all_perms = 0
    # for k in np.arange(1,len(excitations)+1):
    #   all_perms += nchoosek(NumberOfIons,k)

    for k in xrange(len(excitations)):
      # for each excitation
      ex = NumberOfIons - k
      for l in xrange(ex+1):
        # l = number of decays
        if l > 0:
          ex_shifted = nchoosek(ex, l) * decay_prob**l * (1-decay_prob)**(ex-l) * excitations[k]
          excitations_decayed[k+l] += ex_shifted
          excitations_decayed[k] -= ex_shifted

    return excitations_decayed


if __name__ == "__main__":

  hsp = hspace(2,2,2,0)
  parms = parameters(hsp)

  parms.detection_time_ccd = -np.log(0.1)
  parms.detection_time_pmt = -np.log(0.1)
  parms.lifetime = 1

  test_pop = np.array([1,0])
  b = DecayedExcitations_PMT(test_pop, parms)
  print b
  print 'should be: ', [0.9, 0.1]

  # test_pop = np.array([1,0,0])
  # b = DecayedExcitations_PMT(test_pop, parms)
  # print b
  # print 'should be: ', [0.81, 0.18, 0.01]

  a = generate_decay_list([[0,0]])
  b = generate_decay_list(a)

  test_pop = [1,0,0,0]
  c = DecayedPopulations_CCD(test_pop, parms)
