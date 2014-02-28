"""optimizechi: A class for optimizations for process matrices
"""

import numpy as np
import scipy.linalg as lnlg
import numpy.matlib as npml
import time
import scipy.optimize as optimize

import proctom

expm = lnlg.expm
from math import pi

X = np.array([[0,1],[1,0]])
Y = np.array([[0,-1j],[1j,0]])
Z = np.array([[1,0],[0,-1]])
Id = np.array([[1,0],[0,1]])

#Paulis=[Id, X, Y, Z]

Paulis = proctom.getopbase() # added 2010-07-26 by TM: to solve operator definition hickups

def findoptimizedrot(ideal_chi, my_chi, return_chi=False, use_random=True, verbose = True):
    """Finds single qubit operations that maximize the overlap
    between an experimental chi matrix and an ideal chi matrix"""
    try:
        if ideal_chi.shape != my_chi.shape:
            raise RuntimeError("Matrices don't have the same dimensions")
    except AttributeError:
        my_chi = my_chi.chi
    nr_of_ions = np.log2(my_chi.shape[0]) / 2
    # if nr_of_ions > 1:
    #     raise NotImplementedError("Can only optimize one qubit by now")
    if use_random:
        p0 = np.random.random_sample(3*nr_of_ions)
    else:
        p0 = np.zeros(3*nr_of_ions)
    errfunc = lambda p: chirotinfidel(p,ideal_chi, my_chi)
    p1, covmat, infodict, mesg, ier = optimize.leastsq(errfunc, p0, full_output=True)
    fidelity = rotated_fidelity(p1, ideal_chi, my_chi)
    if verbose:
        print "fid: " + str(fidelity)
    if return_chi:
        rot_chi = rotated_chi(p1, my_chi)
        return p1, fidelity, rot_chi
    else:
        return p1, fidelity

def findoptimized_list(ideal_chi, my_list):
    p0 =  np.random.random_sample(3)
    errfunc = lambda p: chirotinfidel_list(p,ideal_chi, my_list)
    p1, covmat, infodict, mesg, ier = optimize.leastsq(errfunc, p0, full_output=True)
    fid_list = []
    for my_chi in my_list:
        fidelity = rotated_fidelity(p1, ideal_chi, my_chi)
        print "fid: " + str(fidelity)
        fid_list.append(fidelity)
    return p1, fid_list

def chirotinfidel_list(rot_vector, ideal_chi, my_list):
    fid_list = []
    for my_chi in my_list:
        fid_list.append(chirotinfidel(rot_vector, ideal_chi, my_chi)[0])
    return [max(fid_list),1,1,1]

def chirotinfidel(rot_vector, ideal_chi, my_chi,return_chi=False):
    "The infidelity function"
    NrOfQubits = int(len(rot_vector)/3)
    # rz = expm(-1j*pi/2 * rot_vector[0] * X)
    # ry = expm(-1j*pi/2 * rot_vector[1] * Y)
    # rx = expm(-1j*pi/2 * rot_vector[2] * X)
    R = lambda rot_vector: np.dot(expm(-1j*pi/2 * rot_vector[2] * X),np.dot(expm(-1j*pi/2 * rot_vector[1] * Y),expm(-1j*pi/2 * rot_vector[0] * X)))

    Rtot = 1
    AllPaulis = [1]
    for k in xrange(NrOfQubits):
        Rtot = npml.kron(Rtot, R(rot_vector[3*k:3*k+3]))
        AllPaulis = proctom.baseappend(AllPaulis, Paulis)

    k = 0
    l = 0
    res1 = np.zeros((4**NrOfQubits,4**NrOfQubits),dtype=np.complex)
    for p_matrix in AllPaulis:
        RA = np.dot(Rtot, p_matrix)
        for p_matrix2 in AllPaulis:
            res0 = np.dot(p_matrix2.conj().T, RA) / 2**NrOfQubits
            res1[k,l] = np.trace(res0)
            l += 1
        l = 0
        k += 1

    chirot = np.dot(res1.T,np.dot(my_chi, res1.conj()))
#    IF1 = np.real(1-np.trace(np.dot(ideal_chi,chirot)))
#    IF1 = np.real(1-np.sum(ideal_chi*chirot.transpose()))
#    IF = np.real(1-np.diag(np.dot(ideal_chi,chirot)))
    IF = np.real(1-np.sum(ideal_chi*chirot.transpose(),axis=1))
    if return_chi:
        return chirot
    else:
        return IF

def rotated_fidelity(rot_vector, ideal_chi, my_chi):
    "calculate the rotated fidelity"
    fid_vec = chirotinfidel(rot_vector, ideal_chi, my_chi)
    fidelity = np.sum(1 - fid_vec)
    return fidelity

def rotated_chi(rot_vector, my_chi):
    my_chi = chirotinfidel(rot_vector, my_chi, my_chi, return_chi=True)
    return my_chi


def optimize_chi_dict(chi_dict, id_chi, mean_key='all', return_chi=False):
    """Takes a dictionary of Chi matrices
    It optimizes a given entry of this dictionary and
    calculates the fidelity of all items of this dictionary with the same operation
    """
    mean_chi = chi_dict[mean_key]
    key_list = chi_dict.keys()
    key_list.sort()
    rot_vec, mean_fid = findoptimizedrot(id_chi, mean_chi)
    fid_dict = {}
    if return_chi:
        chi_dict = {}
    for my_key in key_list:
        my_chi = chi_dict[my_key]
        my_fid = rotated_fidelity(rot_vec, id_chi, my_chi)
        if return_chi:
            chi_dict[my_key] = rotated_chi(rot_vec, my_chi)
        fid_dict[my_key] = my_fid
        print str(my_key) + " fidelity: " +str(my_fid)
    if return_chi:
        return fid_dict, chi_dict
    else:
        return fid_dict

if  __name__ == "__main__":
    # chi = np.array([[0.9,0,0,0],[0, 0.1,0,0],[0,0,0,0],[0,0,0,0]])
    # chi_id = np.array([[0,0,0,0],[0,0,0,0],[0,0,1,0],[0,0,0,0]])
    # findoptimizedrot(chi_id, chi, return_chi = False)

    chi2 = np.zeros((16,16))
    chi2[1,1] = 0.7
    chi2[0,0] = 0.3
    chi2_id = np.zeros((16,16))
    chi2_id[0,0] = 1
    findoptimizedrot(chi2_id, chi2, return_chi = False)
