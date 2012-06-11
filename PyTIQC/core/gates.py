''' define a set of ideal gates for constructing algorithms. This is a unitary operator parallel to simulateevolution. Also include functions to perform state evolution and fidelity evaluation. '''

import numpy as np
from scipy import linalg
import copy
import matplotlib.pylab as pl
from mpl_toolkits.mplot3d import Axes3D
import PyTIQC.tools.quantum_tools as qtls

########################################
# basic matrices
#
# note these are Pauli not Spin matrices (no factors of 1/2)
########################################

X = np.array([[0,1],[1,0]])
Y = np.array([[0,-1j],[1j,0]])
Z = np.array([[1,0],[0,-1]])
H = 1/np.sqrt(2)*np.array([[1,1],[1,-1]])
I = np.array([[1,0],[0,1]])
S = np.array([[1,0],[0,1j]])
T = np.array([[1,0],[0,np.exp(1j*np.pi/4)]])

########################################
# composite gates
#
########################################

def dotN(Us):
    ''' multi-gate dot product '''
    U = Us[0];
    for i in range(1, len(Us)):
        U = np.dot(Us[i], U)
    return U

def kronN(Us):
    ''' multi-gate kron '''
    n = len(Us)
    if n == 1:
        return Us[0]
    elif n == 2:
        return np.kron(Us[0], Us[1])
    else:
        return np.kron(Us[0], kronN(Us[1:]))

def R(gate, n, i):
    ''' perform gate on ith qubit in n-qubit space '''
    # count qubits from right
    i = (n-1)-i

    Is = np.zeros([n,2,2], np.complex64)
    Is[:] = np.array([[1,0],[0,1]])
    Is[i] = gate
    return kronN(Is)

def G(gate, n):
    ''' global gate on n qubits, no expm'''
    N = 2**n
    U = np.zeros([N,N])
    for i in range(n):
        U = U + R(gate, n, i)
    return U

def control(n, U, i, j):
    ''' n = # qubits, U = controlled operator,
        i = control index, j = target index '''

    # convert to counting qubits from the right
    i = (n-1) - i
    j = (n-1) - j

    Is = np.zeros([n,2,2], np.complex64)
    Is[:] = np.array([[1,0],[0,1]])
    I0 = np.copy(Is)
    I0[i] = np.array([[1,0],[0,0]])
    I1 = np.copy(Is)
    I1[i] = np.array([[0,0],[0,1]])
    I1[j] = U
    return kronN(I0) + kronN(I1)

def controlMN(n, U, i, ni, j, nj):
    ''' generalized control gate
        n qubits, U = controlled operator,
        control: ni qubits starting from i,
        target: nj qubits starting from j '''

    # convert to counting qubits from the right
    #i = (n-1) - i
    #j = (n-1) - j

    # control qubit sub-matrix (consecutive)
    C0n = np.ones(2**ni)
    C0n[-1]=0
    C0 = np.diag(C0n)
    # target qubit sub-matrix
    C1n = np.zeros(2**ni)
    C1n[-1]=1
    C1 = np.diag(C1n)

    I0 = []
    k = 0
    while(k<n):
        if k != i: I0.append(I);  k=k+1;
        else:      I0.append(C0); k=k+ni;

    I1 = []
    k = 0
    while(k<n):
        if k == i:   I1.append(C1); k=k+ni;
        elif k == j: I1.append(U);  k=k+nj;
        else:        I1.append(I);  k=k+1;

    # convert to counting qubits in reverse
    I0 = I0[::-1]
    I1 = I1[::-1]

    K1 = kronN(I0)
    K2 = kronN(I1)
    return K1+K2

def swap(n, i, j):
    ''' swap gate between qubits i and j '''
    N = 2**n
    U = np.zeros([N,N])
    # convert to counting qubits from the right
    i = (n-1) - i
    j = (n-1) - j
    for k in range(N):
        bin = dec2bin(k, n)
        bin[i], bin[j] = bin[j], bin[i]
        dec = bin2dec(bin)
        vi = np.zeros([N, 1])
        vi[k,0] = 1
        vo = np.zeros([1, N])
        vo[0,dec] = 1
        U = U + np.kron(vo, vi)
    return U

def Toffoli(n, i, j, k):
    ''' Toffoli gate: CNOT on qubit k controlled by qubits i, j '''
    if n < 3:
        raise ValueError("need at least 3 qubits")
    uni = []

    # in controlMN(..) I assume that qubits i,j are next to each other
    # if they're not, have to fix it. 3 possibilities:
    if abs(i-j) >= 1:
        # just swap i and i+1 if it doesn't collide with j or run off the end
        if i+1 != k and i+1 != n:
            uni.append(swap(n, j, i+1))
            uni.append(controlMN(n, X, i, 2, k, 1))
            uni.append(swap(n, j, i+1))
        # otherwise swap i and j first
        else:
            tmp = i
            i = j
            j = tmp
            uni.append(swap(n, j, i+1))
            uni.append(controlMN(n, X, i, 2, k, 1))
            uni.append(swap(n, j, i+1))
    else:
        # here i, j are already next to each other
        uni.append(controlMN(n, X, i, 2, k, 1))        

    uni = dotN(uni)
    return uni

def Fredkin(n, i, j, k):
    ''' Fredkin gate: SWAP on qubit j,k controlled by qubit i '''
    if n < 3:
        raise ValueError("need at least 3 qubits")
    uni = []

    # need to move j,k to be next to each other like in controlMN
    if abs(j-k) >= 1:
        if j+1 != i and j+1 != n:
            uni.append(swap(n, j+1, k))
            uni.append(controlMN(n, swap(2, 0, 1), i, 1, j, 2))
            uni.append(swap(n, j+1, k))
        else:
            tmp = j
            j = k
            k = tmp
            uni.append(swap(n, j+1, k))
            uni.append(controlMN(n, swap(2, 0, 1), i, 1, j, 2))
            uni.append(swap(n, j+1, k))
    else:
        uni.append(controlMN(n, swap(2, 0, 1), i, 1, j, 2))        

    uni = dotN(uni)
    return uni

def CNOT(N=2, i=0, j=1):
    ''' controlled-not gate '''
    return control(N, X, i, j)

# helper functions

def dec2bin(k, n):
    ''' k = integer, n = number of bits.
    return an array of bits of width n '''
    s = np.binary_repr(k)
    s = s.zfill(n)
    return np.uint8(list(s))

def bin2dec(s):
    ''' s = list of integer digits.
    return decimal representation of s '''
    nstr = ''
    for i in range(len(s)):
        nstr = nstr + str(s[i])
    return int(nstr, 2)

############################################
# gates with arbitrary angles
#
############################################

def Rg(th, ph, n, hiddenions=[]):
    ''' global rotation with arbitrary phase '''
    OP = np.cos(ph) * X + np.sin(ph) * Y
    return Ug(OP, th, n, hiddenions) 

def Rx(th):
    ''' global X rotation '''
    return linalg.expm(-0.5j * th * X)

def Ry(th):
    ''' global Y rotation '''
    return linalg.expm(-0.5j * th * Y)

def Ug(gate, th, n, hiddenions=[]):
    ''' global gate on n qubits'''
    N = 2**n
    U = np.zeros([N,N])
    for i in range(n):
        if i not in hiddenions:
            U = U + R(gate, n, i)
    return linalg.expm(-0.5j * th * U)

def Ri(th, ph, n, i, hiddenions=[]):
    ''' individual rotation on qubit i '''
    if i not in hiddenions:
        OP = np.cos(ph) * X + np.sin(ph) * Y
        return linalg.expm(-0.5j * th * R(OP, n, i) )
    else:
        return G(I, n)/n

def Ux(th, n):
    ''' global X gate '''
    return Ug(X, th, n)

def Uy(th, n):
    ''' global Y gate '''
    return Ug(Y, th, n)

def Uz(th, n, i, hiddenions=[]):
    ''' individual Z gate on qubit i '''
    if i not in hiddenions:
        return linalg.expm(-0.5j * th * R(Z, n, i))
    else:
        return G(I, n)/n

def MS(th, ph, n, hiddenions=[]):
    ''' MS gate with arbitrary phase '''
    # by Daniel
    N = 2**n
    Jx = np.zeros([N,N], np.complex64)
    for i in range(n): 
        if i not in hiddenions:
            Jx += R(X, n, i)
    Jy = np.zeros([N,N], np.complex64)
    for i in range(n): 
        if i not in hiddenions:
            Jy += R(Y, n, i)

    term = np.cos(ph)*Jy + np.sin(ph)*Jx
    
    return linalg.expm( - 0.5j * th/2 * np.dot(term, term) )

def MSx(th, n):
    ''' global MS-x gate '''
    return MS(th, np.pi/4, n)

def MSy(th, n):
    ''' global MS-y gate '''
    return MS(th, -np.pi/4, n)



#########################################
# calculate evolution and get overall unitary
#########################################

def totalUnitary(seq):
    uni = 1
    for op in seq:
        uni = np.dot(op, uni)

    return uni

def calculateevolution(seq, nqubits, y0=None):
    ''' the equivalent of qc.simulateevolution, if desired, construct and simulate a sequence completely with terminology contained within this module. '''
 
    N = nqubits
    numpulses = np.shape(seq)[0]

    Useq = np.around(dotN(seq),10)
    runU = Useq

    # calculate the states
    psi = np.zeros(2**N, np.complex64)
    psi[2**N-1] = 1
    psi = psi.T

    # time evolution
    Y = np.zeros([numpulses+1, 2**N], np.complex128)
    if y0 == None:
        Y[0,:] = psi
    else: 
        Y[0,:] = y0
    i=1
    while i<= numpulses:
      Y[i,:] = np.dot(seq[i-1], Y[i-1,:])
      i = i+1

    # final state
    psiout = np.dot(runU, psi)

    # populations and densmats
    YR = abs(Y)**2
    rho = np.zeros((numpulses+1, 2**N, 2**N), np.complex128)
    for i in range(numpulses+1):
      rho[i,:,:] = np.outer(np.conj(Y[i,:]), Y[i,:])

    return Y, YR, rho

##################
# functions for evaluation
##################

def fidelity(rho1, rho2):
    ''' fidelity =  np.real(np.trace(np.dot(rho1, rho2))) '''
    # this is faster
    return np.real(np.sum(rho1*rho2.transpose()))

def jozsafid(rho1, rho2):
    ''' jozsa fidelity '''
    tmp = qtls.sqrtm_dm(rho2)
    return np.real(np.trace(qtls.sqrtm_dm(np.dot(np.dot(tmp,rho1),tmp)))**2)

def sso(p1, p2):
    ''' squared statistical overlap '''
    # take diagonal if the input looks like a densmat
    if len(np.shape(p1))==2 and len(np.shape(p2))==2:
        p1 = np.diag(p1)
        p2 = np.diag(p2)
    return np.real(np.sum(np.sqrt(p1) * np.sqrt(p2))**2)

def tracedist(rho1, rho2):
    ''' trace distance. Can be calculated for either populations or densmats - decide based on dimensions of input.'''
    if len(np.shape(rho1))==2 and len(np.shape(rho2))==2:
        tmp = rho1-rho2
        return 1-np.real(np.trace(qtls.sqrtm_dm(np.dot(tmp.conjugate().transpose(),tmp))))/2
    else:
        return 1-0.5*np.sum(np.abs(rho1-rho2))

fidelities_dict = {
    'jozsafid': lambda rho1, rho2: jozsafid(rho1, rho2),
    'tracedist-rho': lambda rho1, rho2: tracedist(rho1, rho2),
    'tracedist-pop': lambda rho1, rho2: tracedist(np.diag(rho1), np.diag(rho2)),
    'sso': lambda rho1, rho2: sso(np.diag(rho1), np.diag(rho2))
}
''' map the appropriate fidelity measure to the keyword. '''


##########################
# functions for plotting
##########################

def dispmtx(M, digits=5):
    ''' print the matrix rounded to (default) 5 decimal points. '''
    print np.around(M, digits)

def displaystates(psi, N=5, pop=False):
  ''' print the population of each state '''
  for i in range(len(psi)):
    if np.around(abs(psi[i]), 5) > 0:
        if pop:
            print np.binary_repr(i).zfill(N), ": ", np.around(psi[i],3)
        else:
            print np.binary_repr(i).zfill(N), ": ", np.around(abs(psi[i]**2),3), np.around(psi[i], 5)

def displaytracedstates(psi, N=5, tracemask="11000", pop=False):
  ''' print the population of states, trace out certain qubits '''
  # count number of bits not traced over in tracemask
  # note: assumes N = len(tracemask)
  y = np.zeros(2**tracemask.count('0'))
  # check each index
  for i in range(len(psi)):
      ind = []  # pack the important bits into a char list
      for [j, digit] in enumerate(np.binary_repr(i).zfill(N)):
          if tracemask[j] == '0':
              ind.append(digit) # only pack if digit j is not traced out
      if pop:
          y[ int(''.join(ind), 2) ] += psi[i]
      else:
          y[ int(''.join(ind), 2) ] += abs(psi[i]**2)

  for i in range(len(y)):
    if np.around(abs(y[i]), 5) > 0:
      print i, ': ', np.binary_repr(i).zfill(tracemask.count('0')), ": ", np.around(y[i],3)

  return y

def plotBar3D(rho, displ=1):
  ''' make 3D bar plots of density matrices or unitaries '''

  # this doesn't seem to work get rid of the warnings about dividing by zero
  j = np.seterr(divide='ignore')
  # so do this as a workaround
  cheat =  0.000001
  realrho = np.real(rho)
  realrho[np.nonzero(realrho==0)] = cheat
  imagrho = np.imag(rho)
  imagrho[np.nonzero(imagrho==0)] = cheat
  rhoplot = realrho + imagrho*1j

  dim = len(rho)

  fig = pl.figure(displ,figsize=(12,5))

  xpos, ypos = np.meshgrid(np.arange(dim), np.arange(dim))
  xpos = xpos.flatten()
  ypos = ypos.flatten()
  zpos = np.zeros(len(xpos.flatten()))
  dx = 0.7*np.ones_like(xpos)
  dy = dx.copy()

  cl='y'
  
  rect1 = fig.add_subplot(1,3,1).get_position() 
  ax1 = Axes3D(fig, rect1)
  ax1.bar3d(xpos, ypos, zpos, dx, dy, abs(rhoplot).flatten(), color=cl)
  ax1.set_zlim3d([0,1])
#  ax1.set_zlim3d([-np.pi,np.pi])

  rect2 = fig.add_subplot(1,3,2).get_position() 
  ax2 = Axes3D(fig, rect2)
  ax2.bar3d(xpos, ypos, zpos, dx, dy, realrho.flatten(), color=cl)
  ax2.set_zlim3d([-0.5,0.5])
#  ax2.set_zlim3d([-np.pi, np.pi])

  rect3 = fig.add_subplot(1,3,3).get_position() 
  ax3 = Axes3D(fig, rect3)
  ax3.bar3d(xpos, ypos, zpos, dx, dy, imagrho.flatten(), color=cl)
  ax3.set_zlim3d([-0.5,0.5])

  pl.show()
