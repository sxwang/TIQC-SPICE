''' classes and functions for use by qctools. Defines the Hamiltonian and Noise objects. ''' 

import numpy as np
import numpy.matlib as npml
import scipy.linalg
import scipy.linalg as splg
import simtools as sim
import SchroedingerEqSolvers as ses

#######################################
# class Hamiltonian
#######################################
class Hamilton:
    ''' functions for defining unitary or timedep hamiltonian '''
    def __init__(self):
        pass

    def couplings(self, targetion, omrabi, phase, detuning, wsec, eta, hspace):
        """
        everything is provided in units of omega, so we remove
        that part and only focus on what couples with what -
        this has the advantage that we can later on add a
        timedepending omega(t) without problems
        """
        eip = np.exp(1j*phase)
        dimensions = hspace.dimensions
        Hint = np.zeros([dimensions, dimensions], np.complex128)

        # off-diagonal elements
        for k in range(1, dimensions):
            [p1,s1, st] = indexToState(k, hspace)
            for m in range(k):
                element = 0
                [p2,s2, st] = indexToState(m, hspace)
                ds = s1-s2
                ma = np.nonzero(ds)[0]             # matlab's "find"
                if len(ma) == 1:                # electronic state diff by only one ion
                    ion = np.argmax(abs(ds))       # get the index
                    if np.sum(abs(p1-p2)) == 0 :   # carrier transition
                        element = targetion[ion]
                    elif np.sum(abs(p1-p2)) == 1:  # red/blue sideband
                        mode = 1
                        element = 1j * np.sqrt(max(p1,p2)) * \
                                  eta[ion] * \
                                  np.sign(ds[ion]) * \
                                  targetion[ion]
                Hint[k,m] = element * eip

        return Hint

    def energies(self, targetion, omrabi, phase, detuning, wsec, eta, hspace):
        """
        calculate the energies/diagonal elements. will be added to the coupling
        terms to provide the hamiltonian
        """

        dimensions = hspace.dimensions

        HTdiag = np.zeros(dimensions)
        for k in range(dimensions):
            [p,s] = indexToState(k, hspace)[0:2]
            HTdiag[k] = p*wsec + detuning*np.sum(s)
        HTdiag = np.diag(HTdiag)

        return HTdiag

    def ACshift_corr(self, pulse, params):
        ''' extra term in HT to account for stark shift on S-P for ACstark pulse '''
        hspace = params.hspace
        nuions = params.hspace.nuions
        sys = np.zeros((hspace.levels**hspace.nuions,hspace.levels**hspace.nuions))
        for k in xrange(hspace.nuions): 
            sys += hspace.operator_dict['sigz'][:,:,k] * pulse.targetion[nuions-k-1]
        hz = npml.kron( sys, hspace.operator_dict['id_a'] )
        hz_diag = np.diag(hz)

        return params.ACcorr * pulse.omrabi**2 * hz_diag

    def Hamiltonian(self, pulse, params, LDApprox=True):
        ''' Hamiltonian definition using hspace.operator_dict '''
        opdict = params.hspace.operator_dict

        targetion = pulse.targetion[::-1]
        omrabi = pulse.omrabi
        phase = pulse.phase
        detuning = pulse.detuning
        wsec = params.omz
        eta = params.eta
        hspace = params.hspace

        # prefactor including the timedepending omega and the detuning+phase
        prefac = np.exp(1j * phase) * omrabi * 1/2.

        # coupling to the sideband. argument of the exponent/ld approximation
        LDpart = 1j * eta[-1] * (opdict['a']+opdict['a_dag']) 
        # calculate the coupling based on the addressing errors

        # Lamb-Dicke approximation yes/no?
        # kron with qubits to get the hamiltonian
        # here we'll use the fact that eta may be different for each ion
        sys_ho = np.zeros(((hspace.maxphonons+1)*hspace.levels**hspace.nuions, \
                          (hspace.maxphonons+1)*hspace.levels**hspace.nuions ), \
                              np.complex128)
        if LDApprox:
            for k in xrange(hspace.nuions):
                #etak = eta[-1] if pulse.ion == -1 else eta[k]
                sys_ho += npml.kron( targetion[k]*opdict['raising'][:,:,k],\
                               np.diag(np.ones(hspace.maxphonons+1)) + \
                               1j * eta[k] * (opdict['a']+opdict['a_dag']) )
        else:
            for k in xrange(hspace.nuions):
                #etak = eta[-1] if pulse.ion == -1 else eta[k]
                sys_ho += npml.kron( targetion[k]*opdict['raising'][:,:,k],\
                               splg.expm2(1j * eta[k] * (opdict['a']+opdict['a_dag'])) )

        # multiply with rabi-frequency part
        H = prefac * sys_ho

        # diagonal terms
        sysz = np.zeros((hspace.levels**hspace.nuions,hspace.levels**hspace.nuions))
        for k in xrange(hspace.nuions):
            sysz += opdict['sigz'][:,:,k] * 1/2.
        energies = -detuning * npml.kron(sysz, opdict['id_a']) + \
                   wsec * npml.kron(np.diag(np.ones(hspace.levels**hspace.nuions)), \
                                        np.dot(opdict['a_dag'], opdict['a']) )
        # subtract a zero offset
        energies -= np.diag(energies[0,0]*np.ones_like(np.diag(energies)))

        # add h.c. part
        HT = H + H.transpose().conjugate() + energies

        # diagonal elements of matrix for basis transformation (lasertoqc)
        lqc = np.diag(HT)

        if pulse.type == 'Z':
            HT = HT+np.diag(self.ACshift_corr(pulse, params))
            # to fix discrepancy between lab/Volkmar and here
            if pulse.theta < 0:
                HT = -HT
                lqc = -lqc

        return HT, lqc



    def Hamiltonian_timedep_complete(self, targetion, omrabi, phase, detuning, wsec, eta, hspace, LDApprox = False):
        """
        everything in the full interaction frame of the ions including the
        harmonic oscillation (opposing to the original version from hartmut
        that still had the phonon energy nonzero)

        in addition, this version won't omit the exp(+- i omega_z terms t) in
        attached to creation/annihilation operators

        this is eqn 3.6 in Christian Roo's thesis.

        """

        # H = exp(i delta_c t+phase) omega(t) (sum_i address_error_i sigma_plus_i) (exp(i eta (a exp(-i omegaz t) + adag exp(i omegaz t)))) + h.c.
        # !!! note !!!
        # here eta is assumed to be the same for all qubits. this is
        # not necessarily the case.
        #a more correct version would allow for
        # individual etas for each qubit ... maybe later. not
        # really required here and would only slow down the code

        targetion = targetion[::-1]

        # prefactor including the timedepending omega and the detuning+phase
        prefac = lambda t: np.exp(1j * (detuning * t + phase)) * omrabi(t) * 1/2.

        # coupling to the sideband. argument of the exponent/ld approximation
        LDpart = lambda t:  1j * eta[-1] * (hspace.operator_dict['a'] * np.exp(-1j * wsec * t ) + hspace.operator_dict['a_dag'] * np.exp(1j * wsec * t ) )

        # calculate the coupling based on the addressing errors
        sys = np.zeros((hspace.levels**hspace.nuions,hspace.levels**hspace.nuions))
        for k in xrange(hspace.nuions):
            sys += targetion[k]*hspace.operator_dict['raising'][:,:,k]

        # Lamb-Dicke approximation yes/no?
        # kron with qubits to get the hamiltonian
        if LDApprox:
            sys_ho = lambda t: npml.kron( sys, np.diag(np.ones(hspace.maxphonons+1)) + LDpart(t) )
        else:
            sys_ho = lambda t: npml.kron( sys, splg.expm2( LDpart(t)) )

        # multiply with rabi-frequency part
        H = lambda t: prefac(t) * sys_ho(t)

        # add h.c. part
        HT = lambda t: H(t)+H(t).transpose().conjugate()

        return HT


#################################
# class Noise
#################################
class Noise:
    def __init__(self):
        pass

    def sumFunctions(self, flist):
        ''' return a lambda function as the sum of a list of lambda functions '''
        def ftotal(t):
            return np.sum([func(t) for func in flist], axis=0)
        return ftotal

    def prodFunctions(self, flist):
        ''' return a lambda function as the product of a list of lambda functions '''
        def fprod(t):
            return np.prod([func(t) for func in flist])
        return fprod

    def Noise(self, params, dec):
        ''' generate the noise_dict object for use by qctools '''
        hspace = params.hspace

        noise_dict = {}
        ZeroMtx = lambda t: params.hspace.operator_dict['zero']
        noiseOne = lambda t: 1

        # default none
        noise_dict['none'] = [noiseOne, ZeroMtx, ZeroMtx]
        noise_dict['all'] = [noiseOne, ZeroMtx, ZeroMtx]

        # additive parts
        noise_dict['dephase'] = [noiseOne, self.Noise_dephase(params, dec), ZeroMtx]

        # multiplicative parts
        noise_dict['intensity'] = [self.Noise_intensity(params, dec), ZeroMtx, ZeroMtx]

        # non-unitary part
        # the additive term is the collapse operator. 
        # doesn't work yet for any pulse type besides Delay. needs some debugging.
        # 05/21/2012
#        noise_dict['spontdecay'] = [noiseOne, self.Noise_spontdecayH(params, dec),
#                                    self.Noise_spontdecay(params, dec)]
        noise_dict['spontdecay'] = [noiseOne, ZeroMtx,
                                    self.Noise_spontdecay(params, dec)]

        noise_dict['heating'] = [noiseOne, ZeroMtx,
                                 self.Noise_heating(params, dec)]

        return noise_dict

    def Noise_dephase(self, params, dec):
        ''' generate correlated dephasing error
            due to correlation, must pre-compute using dec.calcDephasing
            and search with searchsorted '''
        hspace = params.hspace
        sys = np.zeros((hspace.levels**hspace.nuions,hspace.levels**hspace.nuions))
        for k in xrange(hspace.nuions):
            sys += hspace.operator_dict['sigz'][:,:,k]
        noise_ho = npml.kron( sys, hspace.operator_dict['id_a'] )

        noise = lambda t: dec.dephaseV[dec.decT.searchsorted(t)] * noise_ho \
                          if dec.decT.searchsorted(t) < len(dec.dephaseV) \
                          else dec.dephaseV[-1] * noise_ho

        return noise

    def Noise_intensity(self, params, dec):
        ''' generate intensity error by randomly drawing a gaussian (uncorrelated) '''
        noise = lambda t: dec.intensfluctV[dec.decT.searchsorted(t)] if dec.decT.searchsorted(t) < len(dec.intensfluctV)-1 else dec.intensfluctV[-1]
        return noise

    def Noise_spontdecay(self, params, dec):
        ''' randomly switch on projective matrix for spontaneous decay by setting HT to zero '''
        ZeroMtx = params.hspace.operator_dict['zero']

        cond = lambda t: bool(dec.spontdecayV[dec.decT == t][0]) if np.any(dec.decT==t) else False

        proj = params.hspace.operator_dict['proj']
        projmtx = lambda t: npml.kron( self.R(proj, params.hspace.nuions, np.random.random_integers(0, params.hspace.nuions-1)), params.hspace.operator_dict['id_a'] )
        
        randion = np.random.random_integers(0, params.hspace.nuions-1)
        jump = params.hspace.operator_dict['lowering'][:,:,randion]
        jumpop = npml.kron( jump, params.hspace.operator_dict['id_a'] )
        jumpH = 0.5j * np.dot(jumpop.conj().transpose(), jumpop)

        noise = lambda t: projmtx(t) if cond(t) else ZeroMtx

        return noise

    def Noise_spontdecayH(self, params, dec):
        ''' non-Hermitian Hamiltonian for spontaneous decay '''
        jump = params.hspace.operator_dict['lowering'][:,:,-1]
        jumpop = npml.kron( jump, params.hspace.operator_dict['id_a'] )
        jumpH = -0.5j * np.dot(jumpop.conj().transpose(), jumpop)

        noise = lambda t: jumpH

        return noise        

    def Noise_heating(self, params, dec):
        ''' randomly switch on projective matrix for heating by setting HT to zero '''
        ZeroMtx = params.hspace.operator_dict['zero']

        cond = lambda t: bool(dec.heatingV[dec.decT == t][0]) if np.any(dec.decT==t) else False

        proj_a = params.hspace.operator_dict['a_dag']
        proj_a[-1,-1] = 1 # make sure we don't exceed the phonon space
        projmtx = lambda t: npml.kron( np.eye(params.hspace.levels**params.hspace.nuions), proj_a )

        noise = lambda t: projmtx(t) if cond(t) else ZeroMtx

        return noise

    def Noise_heatingH(self, params, dec):
        ''' non-Hermitian Hamiltonian for heating '''

        # TODO: fix. 
        # this generates an operator with maxphonons/2*j on the diagonals
        # and causes an "overflow in exp" error in qctools
        jump = params.hspace.operator_dict['a_dag']
        iden = np.eye(2**params.hspace.nuions)
        jumpop = npml.kron( iden, jump )
        jumpH = -0.5j * np.dot(jumpop.conj().transpose(), jumpop)

        noise = lambda t: jumpH

        return noise        

        

    #####################
    ### the following two functions taken from core.gates
    def kronN(self, Us):
        n = len(Us)
        if n == 2:
            return npml.kron(Us[0], Us[1])
        elif n > 2:
            return npml.kron(Us[0], self.kronN(Us[1:]))
        else:
            return Us[0]

    def R(self, gate, n, i):

        Is = np.zeros([n,2,2], np.complex128)
        Is[:] = np.array([[1,0],[0,1]])
        Is[i] = gate

        return self.kronN(Is)
    #####################


######################################
# other useful functions
######################################
def indexToState(ind, hspace):
    """ return # of phonons and binary repr of state """
    maxphonons = hspace.maxphonons
    levels = hspace.levels
    nuions = hspace.nuions
    # phonon number
    phonon = ind % (maxphonons+1)

    def unpackdigits(numstr):
        ''' like np.unpackbits, but for base 10 '''
        numarr = []
        for c in numstr:
            numarr.append(float(c))
        return np.array(numarr)

    # state can be expressed in binary
    indstate = (ind - phonon) / (maxphonons+1)
    # statestr is a string like |SS,0>
    # statenum is an array of 0,1,2's
    statestr = np.base_repr(indstate, base=levels).zfill(nuions)
    statenum = unpackdigits(statestr)
    statestr = statestr.replace('0', 'D').replace('1', 'S').replace('2', 'A')
    statestr = statestr + ',' + str(phonon)

    return phonon, statenum, statestr

def stateToIndex(statestr, hspace):
    """ return index given string repr of state; inverse of indexToState() """
    maxphonons = hspace.maxphonons
    levels = hspace.levels
    nuions = hspace.nuions

    statestr = statestr.replace('D', '0').replace('S', '1').replace('A', '2')
    statestr = statestr.split(',')
    indstate = int(statestr[0],levels) * (maxphonons+1)
    phonon = int(statestr[1])

    index = indstate + phonon

    return index

def indexToExcitations(ind, hspace):
    """ return the number of excitations/D states for a given index"""
    StateArray = indexToState(ind, hspace)[1] # with D = 0, S=1
    # take all qubits, make them bright, subtract the '0' D states
    NumberOfExcitations = int(hspace.nuions - np.sum(StateArray))
    return NumberOfExcitations


def SEsolver(solver):
    """ select ODE solver to use. Default = ZVODE. """
    if solver == "ZVODE":
       return ses.ODE_timeevo
    elif solver == "Cheby":
       return ses.Chebyshev_timeevo
    else:
       print "**** warning: solver not selected!****"
       return None
