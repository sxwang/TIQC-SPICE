#!/usr/bin/env python
# -*- mode: Python; coding: latin-1 -*-
# Time-stamp: "2010-09-19 16:50:12 tomi"

#  file       quantum_tools.py
#  copyright  (c) Thomas Monz 2010

import numpy as np

def PartialTrace(psirho, dim_list, trace_elements):
    """
    given input state (state or density matrix) psirho, defined by
    the dimensions dim_list, trace out the elements of that list
    provided in trace_elements
    e.g. PartialTrace(rho, [2,4,2], [1,2]) will trace out dimensions
    [4,2]
    ---
    counting is done from the left !!! have to fix code later on for that
    """
    psirho = np.array(psirho)
    dim_list = np.array(dim_list)
    # because humans start at one
    trace_elements = np.array(trace_elements)-1
    # check arguments
    if np.any( trace_elements > len(dim_list)-1) or np.any(trace_elements < 0):
        print 'Invalid trace element'

    n = len(dim_list)
    rdim = dim_list[::-1]
    keep = np.delete(np.arange(n),trace_elements)
    dimtrace = np.prod(dim_list[trace_elements])
    dimkeep = len(psirho)/dimtrace

    if len(psirho.shape) == 1:
        # state vector
        perm = n-1-np.array(keep[::-1].tolist()+trace_elements.tolist())
        print perm
        print rdim
        x = np.reshape(np.transpose(np.reshape(psirho, rdim,'FORTRAN'),perm), [dimkeep,dimtrace], 'FORTRAN')
        rho_traced = np.dot(x,x.conjugate().transpose())

    else:
        perm = n-1-np.array(keep[::-1].tolist()+(keep[::-1]-n).tolist()+trace_elements.tolist()+(trace_elements-n).tolist())
        x = np.reshape(np.transpose(np.reshape(psirho, rdim.tolist()+rdim.tolist(),'FORTRAN'), perm), [dimkeep,dimkeep,dimtrace**2], 'FORTRAN')
        x = np.sum(x[:,:,np.arange(1,dimtrace**2+1,dimtrace+1)-1],axis = 2)
        rho_traced = x

    # norm
    rho_traced /= np.trace(rho_traced)

    return rho_traced

def sqrtm_dm(x):
    """
    Compute the square root y of x such that y*y = x, where x is a general
    positive semi-definite matrix and * is standard matrix multiplication.
    Computation is performed by eigenvalue decomposition.

    Will not work for defective inputs (and no error is generated)!
    """

    from scipy.linalg import eigh

    d,e = eigh(x)
    # sometimes we still do get veeery small negative
    # eigenvalues. so i cut off at 10**-12 and set them
    # to zero. that way the sqrt on real number doesn't
    # produce nans
    for k in xrange(len(d)):
        if d[k] < 10**-12:
            d[k] = 0

    d = (d**0.5)
    # replace inv with trans.conj for faster evaluation
    return np.dot(e,np.dot(np.diag(d),e.transpose().conjugate()))

# rho = np.zeros((8,8))
# rho[0,0] = 1
# rho[7,7] = 1
# rho[0,7] = 1
# rho[7,0] = 1
# #rho = np.array([0,0,0,0, 1,0,0,1])

# PartialTrace(rho,[2,2,2],[2,1])


# quantum_tools.py ends here
