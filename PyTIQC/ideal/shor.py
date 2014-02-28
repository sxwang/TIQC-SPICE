''' generate Shor 15 unitaries from ideal gates '''

import PyTIQC.core.gates as U
import numpy as np
import matplotlib.pyplot as pl

pi = np.pi

N = 7

# easy Shor case
Useq = [ \
    U.R(U.H, N, 0), \
    U.R(U.H, N, 1), \
    U.R(U.H, N, 2), \
    # target qubits 3,4,5,6, control on 0,1,2
    U.control(N, U.X, 2, 4), \
    U.control(N, U.X, 2, 5), \
    # Uqft part
    U.R(U.H, N, 0), \
    U.control(N, U.S, 1, 0), \
    U.R(U.H, N, 1), \
    U.control(N, U.T, 2, 0), \
    U.control(N, U.S, 2, 1), \
    U.R(U.H, N, 2) ]
Ushor11 = np.around(U.dotN(Useq), 10)

# hard Shor case
Useq = [ \
    U.R(U.H, N, 0), \
    U.R(U.H, N, 1), \
    U.R(U.H, N, 2), \
    # target qubits 3,4,5,6, control on 0,1,2
    U.control(N, U.X, 2, 4), \
    U.control(N, U.X, 2, 5), \
    # what we want here: Fredkin(N, 0, swap(3,5)), Fredkin(N, 0, swap(4,6))
    U.CNOT(N, 3, 5), \
    U.Toffoli(N, 1, 5, 3), \
    U.CNOT(N, 3, 5), \
    U.CNOT(N, 6, 4), \
    U.Toffoli(N, 1, 4, 6), \
    U.CNOT(N, 6, 4), \
    # Uqft part
    U.R(U.H, N, 0), \
    U.control(N, U.S, 1, 0), \
    U.R(U.H, N, 1), \
    U.control(N, U.T, 2, 0), \
    U.control(N, U.S, 2, 1), \
    U.R(U.H, N, 2) ]
Ushor7 = np.around(U.dotN(Useq), 10)

# calculate the states
psi = np.zeros(2**N, np.complex64)
psi[U.bin2dec('1000000')] = 1
psi = psi.T

psiout = np.dot(Ushor7, psi)

U.displaytracedstates(psiout, N, tracemask='1111000')

# for Ushor11, should get 0,4
# for Ushor7, should get 0,2,4,6
# see Vandersypen's thesis p.207


# now get unitaries for the Kitaev version
N = 5

a11 = [
    U.control(N, U.X, 0, 1), \
    U.control(N, U.X, 0, 3) ]
Ua11 = U.dotN(a11)

a7a = [
    U.control(N, U.X, 0, 2), \
    U.control(N, U.X, 0, 3) ]
Ua7a = U.dotN(a7a)

a7b = [
    U.CNOT(N, 1, 3), \
    U.Toffoli(N, 0, 3, 1), \
    U.CNOT(N, 1, 3), \
    U.CNOT(N, 4, 2), \
    U.Toffoli(N, 0, 2, 4), \
    U.CNOT(N, 4, 2) ]
Ua7b = U.dotN(a7b)

#np.savetxt('Ushor_a11.txt', real(Ua11), fmt='%1d')
#np.savetxt('Ushor_a7a.txt', real(Ua7a), fmt='%1d')
#np.savetxt('Ushor_a7b.txt', real(Ua7b), fmt='%1d')
