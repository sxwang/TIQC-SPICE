#!/usr/bin/env python
# -*- mode: Python; coding: latin-1 -*-
# Time-stamp: "2011-07-10 15:09:55 tomi"

#  file       quantumprocess.py
#  copyright  (c) Thomas Monz 2010

import numpy as np

import proctom
import optimizechi

import numpy.matlib as npml


#import PyTIQC.evaluation.densitymatrixreconstruction.densitymatrix as dm
import densitymatrix as dm

def Unitary2Chi(U, B = None):
    """
    map unitary process U onto the according chi matrix
    option: provide an orthogonal basis B instead of the standard
    Pauli basis for the chi matrix
    """

    NumberOfQubits  = int(np.log2(U.shape[0]))
    chidim = 4**NumberOfQubits

    a = proctom.getopbase()

    A = a

    for k in xrange(NumberOfQubits-1):
        A = proctom.baseappend(A,a)

    lam = np.zeros(chidim, dtype=complex)
    for k in xrange(chidim):
        lam[k] = np.trace(np.dot(A[k].transpose().conjugate(),U))/2**NumberOfQubits

    chi = np.zeros((chidim,chidim), dtype=complex)
    for k in xrange(chidim):
        for l in xrange(chidim):
            chi[k,l]=lam[k]*lam[l].conjugate()
    return chi

def ChoiIn(NumberOfQubits):
    # lets generate a density matrix of pairs of
    # maximally entangled state such as
    # 00+11 (1,3) \otimes 00+11 (2,4)

    dimOrig = int(2**(NumberOfQubits))
    dimChoi = int(2**(2*NumberOfQubits))
    psi = np.zeros((dimChoi,1))
    indices = np.arange(0,dimChoi,dimOrig+1)
    psi[indices] = 1



    rho = np.outer(psi,psi.conjugate())
    # renorm
    rho /= np.trace(rho)

    return rho

def OSumToChoi(chi):
    """
    map process representation in operator sum form
    into choi-jamilkowski density matrix
    representation
    """

    NumberOfQubits = int(np.log2(chi.shape[0])/2)


    choi = ChoiIn(NumberOfQubits)

    id = np.diag(np.ones(2**NumberOfQubits))

    # rewrote it to be faster using kron instead of generating a larger
    # chi matrix

    # generate the choi matrix of the doubled hilbert space
#    chi_choi = np.zeros((4**(2*NumberOfQubits),4**(2*NumberOfQubits)))
    # put the original chi matrix into the top corner so that it
    # acts onto the qubits 1->N but doesn't do anything to N+1->2N
#    chi_choi[0:4**(NumberOfQubits),0:4**(NumberOfQubits)] = chi

#    rho_choi = proctom.proc_channel_output(chi_choi, rho)

    rhoout = np.zeros((2**(2*NumberOfQubits),2**(2*NumberOfQubits)),dtype=complex)

    A = proctom.proc_channel_operators(NumberOfQubits)

    for m in xrange(4**NumberOfQubits):
        for n in xrange(4**NumberOfQubits):
            rhoout += chi[m,n] * np.dot(npml.kron(id,A[:,:,m]), np.dot(choi, npml.kron(id,A[:,:,n]).transpose().conjugate()))

    rho_choi = rhoout


    return rho_choi

def KrausToChoi(Kraus):
    # kraus input should be of the form kraus[k,:,:]
    NumberOfQubits = int(np.log2(Kraus[0,:,:].shape[1]))
    rhoout = np.zeros((2**(2*NumberOfQubits),2**(2*NumberOfQubits)),dtype=complex)
    NumOfOperators = Kraus.shape[0]

    id = np.diag(np.ones(2**NumberOfQubits))

    choi = ChoiIn(NumberOfQubits)

    for k in xrange(NumOfOperators):
        rhoout += np.dot(npml.kron(id,Kraus[k,:,:]), np.dot(choi, npml.kron(id,Kraus[k,:,:]).transpose().conjugate()))

    return rhoout


class QuantumProcess:
    def __init__(self, chi_matrix):
        self.chi_matrix = chi_matrix
        self.choi_matrix = OSumToChoi(chi_matrix)
        self.choi_matrix_obj = dm.DensityMatrixObject(self.choi_matrix)

    def fid(self, id_chi):
        return self.cj_fidelity(id_chi)

    def fidelity(self, id_chi):
        return np.trace(np.dot(self.chi_matrix, id_chi))

    def cj_fidelity(self, id_chi):
        choi_matrix_id = OSumToChoi(id_chi)
        return self.choi_matrix_obj.jozsafid(idealrho=choi_matrix_id)

    def cj_distance(self, id_chi):
        choi_matrix_id = OSumToChoi(id_chi)
        return self.choi_matrix_obj.trdistance(idealrho=choi_matrix_id)


def proctomo_obj(otherself, data, NumberOfIterations):
    """Processtomography which returns a QProcess Object for use with MonteCarlo analysis"""
    chi_data = proctom.proctomo(data, NumberOfIterations)
    my_obj = QuantumProcess(chi_data)
    return my_obj

def proctomo_trace(data, NumberOfIterations = 100, trace_list = [], id_chi=None, use_bell_basis=False):
    """Function for processtomography when you trace over ions
    if id_chi is given it tries to optimize the process tomography with local operations"""
    import qpython.tools.cameratools as ct
    dat_traced = ct.traceions(data,trace_list)
    pt = proctomo(dat_traced, NumberOfIterations=NumberOfIterations, use_bell_basis=use_bell_basis)
    if id_chi != None:
        print "unrotated fid: " + str(abs(pt[0,0]))
        optimizechi.findoptimizedrot(id_chi, pt)
    return pt



if __name__ == "__main__":
    chi = np.zeros((4**1,4**1))
    # chi[0,0] = .5
    # chi[3,3] = 0.5
    chi[2,2] = 1

    choi = OSumToChoi(chi)
#    print choi

    p = 0.5
    mykraus = np.zeros((2,2,2), dtype=complex)
    mykraus[0,:,:] = np.array([[0,np.sqrt(p)],[0,0]])
    mykraus[1,:,:] = np.array([[1,0],[0,np.sqrt(1-p)]])

    choi = KrausToChoi(mykraus)

    ew, ev = np.linalg.eig(choi)
    for k in xrange(ev.shape[0]):
        if ew[k] > 10e-5:
            print ew[k]
            print ew[k]*ev[:,k].reshape(2,2).transpose()
            # op=np.zeros((2,2), dtype = complex)
            # print 2*ew[k]*np.outer(ev[:,k],ev[:,k].transpose().conjugate())

    # chi_id = np.zeros((4**1,4**1))
    # chi_id[0,0]=.5
    # chi_id[3,3]=0.5
    # chi_id[1,1]=0.

#    choi1 = OSumToChoi(chi)
#    choi2 = OSumToChoi(chi_id)

    # import qpython.tools.quantum_tools as qtls
    # print choi1-choi2
    # print np.trace(qtls.sqrtm_dm(np.dot(np.dot(qtls.sqrtm_dm(choi2),choi1),qtls.sqrtm_dm(choi2))))**2

    # testqp = QuantumProcess(chi)
    # print testqp.cj_fidelity(chi_id)
    # print testqp.fidelity(chi_id)
    # print testqp.cj_distance(chi_id)

    # UToff_TM_normal=np.array([[1, 0, 0, 0, 0, 0, 0, 0],
    #                       [0, 1, 0, 0, 0, 0, 0, 0],
    #                       [0, 0, 1, 0, 0, 0, 0, 0],
    #                       [0, 0, 0, 0, 0, 0, 0, 1j],
    #                       [0, 0, 0, 0, np.exp(1j*(-1./np.sqrt(2))*np.pi), 0, 0, 0],
    #                       [0, 0, 0, 0, 0, np.exp(1j*(-1./np.sqrt(2))*np.pi), 0, 0],
    #                       [0, 0, 0, 0, 0, 0, np.exp(1j*(-1./np.sqrt(2))*np.pi), 0],
    #                       [0, 0, 0, np.exp(1j*(-1./2-1./np.sqrt(2))*np.pi), 0, 0, 0, 0]])

    # Unitary2Chi(UToff_TM_normal)

    # a = proctom.getopbase()
    # u = a[3]
    # chiout = Unitary2Chi(u)
    # print chiout



# quantumprocess.py ends here
