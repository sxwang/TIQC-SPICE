#!/usr/bin/env python
# -*- mode: Python; coding: latin-1 -*-
# Time-stamp: "2011-05-30 21:30:10 tomi"

#  file       proctom.py
#  copyright  (c) Thomas Monz 2009

# from pylab import *
# from scipy.linalg import *
#
#
# The evaluation for 1 qubit is different than the one used in matlab !!!!
# This on should be the correct new one
"""proctom.py

This is a module for Quantum Process tomography. By default it uses th cprb camera data.
following functions are available

proctomo(data, NumberOfIterations = 100, S0 = None)

proctomo_trace(data, NumberOfIterations = 100, S0 = None, trace_list = [])

Further functions not intended for use with the interactive shell:

QProcess
"""
import numpy as np
import scipy.linalg as lnlg
import numpy.matlib as npml
import time


import types

try:
    import qpython
except:
    pass


import PyTIQC

try:
    import qpython.tools as qtools
    random_matrices = qtools.random_matrices
except:
    pass
#    import random_matrices




# ----------------------------------------
# ----------------------------------------
# ----------------------------------------

def getopbase():
    X = np.array([[0,1],[1,0]])
    Y = np.array([[0,-1j],[1j,0]]) # proper hermitian basis
    Z = np.array([[1,0],[0,-1]])
    Id = np.array([[1,0],[0,1]])

    Paulis=[Id, X, Y, Z]

    return Paulis

def baseappend(opin, opappend):
    len_opin = len(opin)
    len_opappend = len(opappend)
    opout = []
    for k in xrange(len_opappend):
        for l in xrange(len_opin):
            opout.append(npml.kron(opappend[k],opin[l]))
    return opout

# ----------------------------------------
# ----------------------------------------
# ----------------------------------------


def proctomo(data, NumberOfIterations = 100, use_bell_basis=False):
    """Process tomography for an arbitrary number of ions

    Params
    data: either datafile string or a a cprb matrix
    NumberOfIterations: Number of iterations for the maximum likelihood algorithm
    """

    Paulis = getopbase()

    NumberOfIons = np.int(np.log2(data.shape[1]-1))

    RhoIn, Obs, ObsVal, NRho, NObs = LoadProctomData(data)
#  tic = time.time()
    RhoTObs = np.zeros((NRho*NObs,4**NumberOfIons,4**NumberOfIons),dtype=complex)
    TransposedRhoTObs = np.zeros((NRho*NObs,4**NumberOfIons,4**NumberOfIons),dtype=complex)
    for m in xrange(NRho):
        for n in xrange(NObs):
            tmp = npml.kron(RhoIn[m,:,:].transpose(),Obs[n,:,:])
            RhoTObs[m*NObs+n,:,:] = tmp
            TransposedRhoTObs[m*NObs+n,:,:] = tmp.transpose()
    del RhoIn, Obs, tmp

#            TransposedRhoTObs[m*NObs+n,:,:] = transpose(RhoTObs[m*NObs+n,:,:])
#  print time.time()-tic,'seconds'

    # reserving some more memory space
    QOps = np.zeros((2**NumberOfIons,2**NumberOfIons,(2**NumberOfIons)**2))
    QOps2 = np.zeros((2**NumberOfIons,2**NumberOfIons,(2**NumberOfIons)**2), dtype = complex)
    AB = np.zeros((4**NumberOfIons,4**NumberOfIons), dtype = complex)
    AA = np.zeros((4**NumberOfIons,4**NumberOfIons), dtype = complex)

    v = np.eye(2**NumberOfIons);
    # Quantenoperatoren
    for m in xrange(2**NumberOfIons):
        for n in xrange(2**NumberOfIons):
            QOps[:,:,n+2**NumberOfIons*m] = np.outer(v[:,m],v[:,n])

    for k in xrange(4**NumberOfIons):
        Op_tmp = 1
        for l in xrange(NumberOfIons):
            Op_tmp = npml.kron(Paulis[np.mod(k/4**l,4)],Op_tmp)
        QOps2[:,:,k]=Op_tmp


    for m in xrange(4**NumberOfIons):
        for n in xrange(4**NumberOfIons):
            AB[m,n] = np.sum(QOps[:,:,m]*np.transpose(QOps2[:,:,n]))
            AA[m,n] = np.sum(QOps[:,:,m]*np.transpose(QOps[:,:,n]))

    C = np.dot(lnlg.inv(AB), AA)
    del AB, AA, QOps, QOps2

    # --------------------------

    dimH = 2**NumberOfIons
    dimK = dimH
    idH = np.eye(dimH)
    idK = np.eye(dimK)

    S0 = 1.*npml.kron(idH,idK)/dimK

    Kstart = np.zeros((4**NumberOfIons,4**NumberOfIons))
    # --------------------------

    #print RhoTObs.imag.max()

    ObsValCalc = np.zeros(NRho*NObs)
    #S = np.zeros(((2**NumberOfIons)**2,(2**NumberOfIons)**2))

#    tic = time.time()
    for k in xrange(NumberOfIterations):
        ObsValCalc2 = np.real(np.sum(np.sum(npml.multiply(S0, TransposedRhoTObs), axis = 1), axis = 1))
        # for mn in xrange(NRho*NObs):
        #     ObsValCalc[mn] = np.sum(S0*np.transpose(RhoTObs[mn,:,:]))
        # alternative: this seems to be a factor of 2 faster for 2 ions, but becomes a factor of two slower for 3 ions ...
#         S0_long = tile(S0,(NRho*NObs,1,1))
#         ObsValCalc = sum(sum(S0_long*TransposedRhoTObs,axis=2),axis=1)

        # tensordot(a,b,(0,0)) does something like sum_i a[i]*b[i,:,:]
        K = np.tensordot(ObsVal/ObsValCalc2,RhoTObs,(0,0))

        # lagrange multiplication
        lamquad = __ptrace(np.dot(K,np.dot(S0,K)),dimH,dimK,2)
        # here is some complain about real/imag definitions ...
        laminv = lnlg.inv(lnlg.sqrtm(lamquad))
        Laminv = npml.kron(laminv,idK)

        # new s-matrix
        S = np.dot(Laminv,np.dot(K,np.dot(S0,np.dot(K,Laminv))))
        S0 = S

#    print time.time()-tic,'seconds'
    # calculate corresponding chi matrix
    # all the info is in the S matrix

    V = np.zeros(((2**NumberOfIons)**2,(2**NumberOfIons)**2))
    for q in xrange(2**NumberOfIons):
        for m in xrange(2**NumberOfIons):
            V[:,q+2**NumberOfIons*m] = npml.kron(v[:,q],v[:,m])

    Chi = np.zeros(((2**NumberOfIons)**2,(2**NumberOfIons)**2),dtype = complex)
    for p in xrange(4**NumberOfIons):
        for q in xrange(4**NumberOfIons):
            Chi[p,q] = np.dot(np.conjugate(np.transpose(V[:,p])),np.dot(S,V[:,q]))

    Chi_final = np.dot(C,np.dot(Chi,np.transpose(np.conjugate(C))))
    return Chi_final

# ----------------------------------------
# ----------------------------------------
# ----------------------------------------


def LoadProctomData(data):
    """
    loads data and does some calculations on input/output states, etc
    mainly a helper-function
    """

    NumberOfIons = np.int(np.log2(data.shape[1]-1))
    dat = data[:,1:]

    Id = np.array([[1,0],[0,1]])
    X = np.array([[0,1],[1,0]])
    Y = np.array([[0,-1j],[1j,0]])
    Z = np.array([[1,0],[0,-1]])

    NRho = 4**NumberOfIons # 4 output dm's per ion (S,S+D,S+iD,D)
    NObs = 6**NumberOfIons # the 6 basis we measure: z, -z, x, -x, y, -y

    # reserve some memory
    RhoIn = np.zeros((NRho,2**NumberOfIons,2**NumberOfIons),dtype=complex)
    Obs = np.zeros((NObs,2**NumberOfIons,2**NumberOfIons),dtype=complex)
    ObsVal = np.zeros(NObs*NRho)

    # now we'll smartly figure out the thetas/rhos of the rotations
    theta_proc = np.array([0, 0.5, 0.5, 1])*np.pi
    phi_proc = np.array([0, 0, 0.5, 0])*np.pi

    # a bit tricky to get the right thetas for the ions now
    theta_list = []
    for k in xrange(NumberOfIons):
        theta_tmp = 1
        for l in xrange(NumberOfIons):
            if k==l:
                theta_tmp = npml.kron(theta_proc,theta_tmp)
            else:
                theta_tmp = npml.kron(np.ones(4),theta_tmp)
        theta_list.append(theta_tmp)

    phi_list = []
    for k in xrange(NumberOfIons):
        phi_tmp = 1
        for l in xrange(NumberOfIons):
            if k == l:
                phi_tmp = npml.kron(phi_proc,phi_tmp)
            else:
                phi_tmp = npml.kron(np.ones(4),phi_tmp)
        phi_list.append(phi_tmp)

    # states and matrixes
    psistart = np.zeros(2**NumberOfIons)
    psistart[2**NumberOfIons-1] = 1

#  print (theta_list,phi_list)
#  return (theta_list,phi_list)
    for k in xrange(NRho):
        Op_tmp = np.zeros((NumberOfIons,2,2),dtype=complex)
        for l in xrange(NumberOfIons):
            Op_tmp[l,:,:] = lnlg.expm(1j*theta_list[l][k]/2*(np.cos(phi_list[l][k])*X-np.sin(phi_list[l][k])*Y))

        Op = Op_tmp[0,:,:]
        if NumberOfIons > 1:
            for l in xrange(NumberOfIons-1):
                Op = npml.kron(Op,Op_tmp[l+1,:,:])

        psi_in = np.dot(Op,psistart)
        RhoIn[k,:,:] = np.outer(psi_in,np.conjugate(psi_in))
    #print RhoIn

    # new lets start to work on teh projectors
    ez_list = []
    for k in xrange(2**NumberOfIons):
        ez_list.append(np.eye(2**NumberOfIons)[:,k])

    thetaIXY_state = np.array([0, 0.5, 0.5])*np.pi
    phiIXY_state = np.array([0, 1.5, 1.0])*np.pi

    theta_state_list = []
    for k in xrange(NumberOfIons):
        theta_tmp = 1
        for l in xrange(NumberOfIons):
            if k==l:
                theta_tmp = npml.kron(thetaIXY_state,theta_tmp)
            else:
                theta_tmp = npml.kron(np.ones(3),theta_tmp)
        theta_state_list.append(theta_tmp)
#    theta_state_list.insert(0,theta_tmp)

    phi_state_list = []
    for k in xrange(NumberOfIons):
        phi_tmp = 1
        for l in xrange(NumberOfIons):
            if k == l:
                phi_tmp = npml.kron(phiIXY_state, phi_tmp)
            else:
                phi_tmp = npml.kron(np.ones(3),phi_tmp)
        phi_state_list.append(phi_tmp)
#    phi_state_list.insert(0,phi_tmp)

    for k in xrange(np.int(NObs/2**NumberOfIons)):
        Op_tmp = np.zeros((NumberOfIons,2,2),dtype=complex)
        for l in xrange(NumberOfIons):
            Op_tmp[l,:,:]=lnlg.expm(-1j*theta_state_list[l][k]/2*(np.cos(phi_state_list[l][k])*X-np.sin(phi_state_list[l][k])*Y));

        Op = Op_tmp[0,:,:]
        if NumberOfIons > 1:
            for l in xrange(NumberOfIons-1):
#         Op = npml.kron(Op_tmp[l+1,:,:],Op)
                Op = npml.kron(Op_tmp[l+1,:,:],Op)

#        print Op
#        sqrtm(Op)

        for m in xrange(2**NumberOfIons):
            psi = np.dot(Op,np.eye(2**NumberOfIons)[:,m])
            Obs[k*(2**NumberOfIons)+m,:,:] = np.outer(psi,psi.conjugate())
#    print Obs


    # calculate obsval
    for k in xrange(NRho):
        for l in xrange(np.int(NObs/(2**NumberOfIons))):
            obsind = NObs*k + (2**NumberOfIons)*l;
            datind = NObs/(2**NumberOfIons) * k + l;
            for m in xrange(2**NumberOfIons):
                ObsVal[obsind+m]=dat[datind,m]
    return (RhoIn, Obs, ObsVal, NRho, NObs)


# ----------------------------------------
# ----------------------------------------
# ----------------------------------------

def __ptrace(C,d1,d2,gamma):
    """ helper function for proctomo - don't touch it"""
    D = d1*d2

    if gamma == 2:
        ind = np.reshape(np.conjugate(np.transpose(np.reshape(np.arange(D),(d2,d1),'FORTRAN'))),(1,D),'FORTRAN')
        C = C[np.ix_(ind.tolist()[0],ind.tolist()[0])]
    else:
        d3 = d2
        d2 = d1
        d1 = d3

    Ct = 0
    L = np.arange(d1)
    for k in xrange(d2):
        ind = L + k*d1
        indl = ind.tolist()
        Ct += C[np.ix_(indl,indl)]
    return Ct

# ----------------------------------------
# ----------------------------------------
# ----------------------------------------

def proc_channel_operators(NumberOfIons):
    """
    return the operators A that we need for our quantum channels
    """

    Paulis = getopbase()

    A = np.zeros((2**NumberOfIons, 2**NumberOfIons, 4**NumberOfIons),dtype=complex)

    # create operators
    for k in xrange(4**NumberOfIons):
        tmp = 1
        for l in xrange(NumberOfIons):
            # print np.mod(int(1.*k/4**l),4)
            tmp = npml.kron(tmp,Paulis[np.mod(int(1.*k/4**l),4)])
        A[:,:,k] = tmp

    return A

#A=proc_channel_operators(2)
# print A

def proc_channel_output(chi, input_rho):

    NumberOfIons = int(np.log(len(chi))/np.log(4))

    A = proc_channel_operators(NumberOfIons)

    if isinstance(input_rho, np.ndarray):
        NumberOfOutputRhos = 1
    if isinstance(input_rho, types.ListType):
        NumberOfOutputRhos = len(input_rho)

    if NumberOfOutputRhos == 1:
        # create ouput
        rhoout = np.zeros((2**NumberOfIons,2**NumberOfIons),dtype=complex)

        for m in xrange(4**NumberOfIons):
            for n in xrange(4**NumberOfIons):
                rhoout += chi[m,n] * np.dot(A[:,:,m], np.dot(input_rho,A[:,:,n].transpose().conjugate()))

        list_rhoout = rhoout

    else:
        list_rhoout = []
        for o in xrange(NumberOfOutputRhos):
            rho_in = input_rho[o]
            rhoout = np.zeros((2**NumberOfIons,2**NumberOfIons),dtype=complex)

            for m in xrange(4**NumberOfIons):
                for n in xrange(4**NumberOfIons):
                    rhoout += chi[m,n]* np.dot(A[:,:,m], np.dot(rho_in,A[:,:,n].transpose().conjugate()))

            list_rhoout.append(rhoout)


    return list_rhoout



# ----------------------------------------
# ----------------------------------------
# ----------------------------------------

def mean_process_fidelity(chi_exp, chi_id, NumberOfSampleStates = 100, type_rho_distribution = 'PureHaar'):
    """
    mean process fidelity: for a given number of random states,
    calculate the jodza-fidelity between the obtain and expected
    output state

    possible distributions random states can be picked from:
       'Haar':  random, pure states according to the Haar measure
       'HilSch': random states according to the Hilbert-Schmidt measure
    """
    dist_type = type_rho_distribution
    func = {}
    func['PureHaar'] = random_matrices.randomPureRho_Haar
    func['HilSch'] = random_matrices.randomRho_HilbertSchmidt
    func['ParTr_n2'] = random_matrices.randomRho_PartialTraceN2

    NumberOfIons = int(np.log(len(chi_exp))/np.log(4))

    rho_in_list = []
    for k in xrange(NumberOfSampleStates):
        rho_in_list.append(func[type_rho_distribution](2**NumberOfIons))

    rho_out_exp_list = proc_channel_output(chi_exp, rho_in_list)
    rho_out_id_list = proc_channel_output(chi_id, rho_in_list)

    fid_list = []
    for k in xrange(NumberOfSampleStates):
        tmp = PyTIQC.evaluation.densitymatrixreconstruction.densitymatrix.DensityMatrixObject(rho_out_exp_list[k])

        fid = tmp.jozsafid(rho_out_id_list[k])

        fid_list.append(fid)

    print 'Process fidelity: %0.4g' % (np.sum(chi_exp* chi_id.transpose()))
    print 'Mean fidelity ('+type_rho_distribution+'): %0.4g +- %0.4g' % (np.mean(np.array(fid_list)), np.std(np.array(fid_list)))


# ----------------------------------------
# ----------------------------------------
# ----------------------------------------



def demo():
    dat_2ions_id = np.loadtxt('2ion_proctomo_id.dat')
    dat_3ions_id = np.loadtxt('3ion_proctomo_id.dat')

    chi2exp = proctomo(dat_2ions_id, NumberOfIterations = 50)
#  chi3id = proctomo(dat_3ions_id, NumberOfIterations = 50)

    if chi2exp[0,0] > 0.95:
        print 'passed'

    #chi1id = np.zeros((4**1,4**1))
    #chi1id[0,0] = 1
    # chi1exp = np.array([
    #     [.5, 0, 0, -0.5j],
    #     [0, 0, 0, 0],
    #     [0, 0, 0, 0],
    #     [0.5j, 0, 0, 0.5],
    #     ])
    chi1exp = np.array([
    [0.25, 0,   0,  0.25],
    [0,   0.25, -0.25j, 0],
    [0,   0.25j, 0.25, 0],
    [0.25, 0,   0,  0.25],
    ])

    # chi1exp = np.array([
    # [0.5, 0,   0,  -0.5j],
    # [0,   0., 0., 0],
    # [0,   0., 0., 0],
    # [0.5j, 0,   0,  0.5],
    # ])


    # a = mean_process_fidelity(chi1exp, chi1id,type_rho_distribution = 'PureHaar',NumberOfSampleStates = 10)

    rho1 = np.array([[0.5, 0.5],[.5, .5]])
    rhoo = proc_channel_output(chi1exp, rho1)

    print rhoo

    a = getopbase()
    A = baseappend(a,a)

    print np.all(npml.kron(a[0],a[1]) == A[1])




if __name__ == "__main__":
    print
    demo()



# proctom.py ends here
