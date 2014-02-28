''' CNOT on 2 out of 3 ions based on Volckmarizer sequence '''

import PyTIQC.core.gates as U
import numpy as np

import matplotlib.pyplot as pl

pi = np.pi

N = 3 #number of qubits

# the Innsbruck pulse sequence
Useq = [ \
    U.Ux(pi/2, N),
    U.Uz(pi/2, N, 0),
    U.MS(pi/4, 0, N),
    U.Ux(pi/4, N),
    U.Uz(pi, N, 2),
    U.Ux(pi/4, N),
    U.MS(pi/4, 0, N),
    U.Uz(pi/2, N, 0),
    U.Ux(pi/2, N),
    U.Uz(pi, N, 2) ]

Ucnot = U.dotN(Useq)
Ucnot = np.around(Ucnot, 10)

# calculate the states
psi = np.zeros(2**N, np.complex64)
psi[-1] = 1
psi = psi.T

# time evolution
Y = np.zeros([2**N, len(Useq)+1], np.complex128)
Y[:,0] = psi.T
i=1
while i<=len(Useq):
  Y[:,i] = np.dot(Useq[i-1], Y[:,i-1])
  i = i+1

# final state
psiout = np.dot(Ucnot, psi)

# print the results
U.displaystates(psiout, N)
pl.plot(abs(Y.T)**2)
pl.show()
