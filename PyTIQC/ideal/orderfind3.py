''' calculate the unitary operator and state evolution for 
 the 5-qubit order finding algorithm '''

import PyTIQC.core.gates as U
import numpy as np
import matplotlib.pyplot as pl

pi = np.pi

# some permutations. Select one of these
#Perm = U.swap(2,0,1)  # Fredkin
#Perm = U.CNOT(2,1,0)   # order1
#Perm = np.kron(U.X, U.I)   # order2
Perm = np.array([[1,0,0,0],[0,0,1,0],[0,0,0,1],[0,1,0,0]],np.complex64) # order3
#Perm = np.array([[0,1,0,0],[0,0,1,0],[0,0,0,1],[1,0,0,0]],np.complex64) # order4
#Perm = U.Ux(pi/2, 2)

Perm2 = np.dot(Perm, Perm)
Perm4 = np.dot(Perm2, Perm2)

# define the gate sequence
# currently counting qubits from the RIGHT
# i.e. "01000" == |y> = |01>, H and QFT on |000>, input state = int(01000,2) = 8
Useq = [ \
  U.R(U.H, 5, 0), \
  U.R(U.H, 5, 1), \
  U.R(U.H, 5, 2), \
  # target qubits 3,4, control on 0,1,2
  U.controlMN(5, Perm, 2, 1, 3, 2), \
  U.controlMN(5, Perm2, 1, 1, 3, 2), \
  U.controlMN(5, Perm4, 0, 1, 3, 2), \
  # Uqft part
  U.R(U.H, 5, 0), \
  U.control(5, U.S, 1, 0), \
  U.R(U.H, 5, 1), \
  U.control(5, U.T, 2, 0), \
  U.control(5, U.S, 2, 1), \
  U.R(U.H, 5, 2)]#, \
  #U.swap(5, 0, 2)  ]

Uqft = [ \
  U.R(U.H, 5, 0), \
  U.control(5, U.S, 1, 0), \
  U.R(U.H, 5, 1), \
  U.control(5, U.T, 2, 0), \
  U.control(5, U.S, 2, 1), \
  U.R(U.H, 5, 2), \
  U.swap(5, 0, 2)  ]


# calculate overall unitary
Uof = U.dotN(Useq)
Uof = np.around(Uof, 10)


### batch process a set of 4 states
psilist = [0, 8, 16, 24] # binary for 00000,01000,10000,11000
for ind in psilist:
    # calculate the states
    psi = np.zeros(32, np.complex64)
    psi[ind] = 1
    psi = psi.T

    # time evolution
    [Y,YR,rho] = U.calculateevolution(Useq, 5, y0=psi.T)

    # final state
    psiout = np.dot(Uof, psi)

    # print the results
    print
    print "input state", np.binary_repr(ind).zfill(5)
    #U.displaystates(psiout)
    #print "traced"
    U.displaytracedstates(psiout, 5, tracemask='11000')

