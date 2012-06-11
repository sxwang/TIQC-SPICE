#!/usr/bin/env python
# -*- mode: Python; coding: latin-1 -*-
# Time-stamp: "2011-05-27 14:46:12 c704252"

#  file       SchroedingerEqSolvers.py
#  author     Thomas Monz

# starting in singer et al.
# quolloquium: trapped ions as qubits: essential numerical tools

# for correct chebychev
# Numerical approaches to time evolution of complex quantum systems
# arXiv:0907.3022

# there is something wrong/off with the factors 1j in killians paper
# with the numerical approach paper definitions, everything works fine

import numpy as np
import numpy.matlib as npml
import scipy.integrate as scint
import pylab as pl
import scipy.linalg as splnlg
import scipy.special as spsp
import time

# ------------------
# time evolution for time depending hamiltonians
# 1: standard ode
# 2: chebyshev

class ODE_timeevo:
    """ standard ODE solver included in Scipy: ZVODE, a complex-valued variable-coefficient ODE sovler. It uses the BDF method with adaptive timesteps. For more details see http://www.netlib.org/ode/zvode.f. """
    def __init__(self, H, tstart, tend, delta_t, psi_in, pbar, pulsestarttime, pulseendtime):
        self.H = H
        self.tstart = tstart
        self.tend = tend
        self.delta_t = delta_t
        self.psi_in = psi_in
        self.pbar = pbar
        self.pstarttime = pulsestarttime
        self.pendtime = pulseendtime

        self.psidot = lambda t, psi: -1j*np.dot(H(t), psi)
        self.jacobi = lambda t, psi: -1j*H(t)

        self.precalc()
        self.do_timesteps()

    def precalc(self):
        self.time = np.arange(self.tstart,self.tend, self.delta_t)
        # arange gives values in range [start,stop) with delta_t ... so the last one is missing
        # add by hand 
        self.time = np.append(self.time, self.tend)
        self.steps = len(self.time)

        self.dim = len(self.psi_in)

        # reserve some mem for states
        self.y = np.zeros((self.dim, self.steps), dtype = 'complex128')
        self.pop = np.zeros((self.dim, self.steps))

        self.y[:,0] = self.psi_in
        self.pop[:,0] = np.abs(self.y[:,0])**2


    def do_timesteps(self):
        eq = scint.ode(self.psidot,jac = self.jacobi).set_integrator('zvode')
        eq.set_initial_value(self.psi_in, self.tstart)

        k=1

        while eq.successful() and eq.t < self.tend:
            eq.integrate(eq.t+self.delta_t)
            if not eq.successful():
                print "ODE solver failed, exiting"
                break
            if k >= np.shape(self.y)[1]:
                #print "warning: index exceeded, exiting"
                break
            self.y[:,k] = eq.y
            self.pop[:,k] = np.abs(eq.y)**2

            if self.pbar:
                self.pbar.update(int(1.*(self.time[k]-self.pstarttime)*100/(self.pendtime-self.pstarttime)))
            k += 1 # increase running index

        # something bad happens if the pulse endtime happens to come very close to a data.Y point: the last point doesn't get updated. To get around this:
        if np.sum(np.abs(self.y[:,-1])**2) == 0 and np.shape(self.y)[1] > 1:
            self.y[:,-1] = self.y[:,-2]

class Chebyshev_timeevo:
    """ Chebyshev ODE solver. for details see "Numerical approaches to time evolution of complex quantum systems",  arXiv:0907.3022. """
    def __init__(self, H, tstart, tend, delta_t, psi_in, chebyshev_order=7, alpha = 0.01, progbar=True):
        self.H = H
        self.tstart = tstart
        self.tend = tend
        self.delta_t = delta_t
        self.psi_in = psi_in
        self.chebyshev_order = chebyshev_order
        self.alpha = alpha
        self.printprogbar = progbar

        self.dim = len(self.psi_in)

        self.time = np.arange(self.tstart,self.tend, self.delta_t)
        # arange gives numbers in [tstart,tend)
        # so we add tend by hand
        self.time = np.append(self.time, self.tend)
        self.steps = len(self.time)

        # reserve some mem for states
        self.y = np.zeros((self.dim, self.steps), dtype = 'complex128')
        self.pop = np.zeros((self.dim, self.steps))

        self.y[:,0] = self.psi_in
        self.pop[:,0] = np.abs(self.psi_in)**2

        self.precalc()
        self.do_definitions()

        self.do_timeevo()


    def precalc(self):
        t1 = time.time()
        # get max/min eigenvalues of all H(t)
        self.Emin = 0
        self.Emax = 0
        for k in self.time:
            eig_energ = splnlg.eigvalsh(self.H(k))

            if self.Emin > np.min(eig_energ):
                self.Emin = np.min(eig_energ)
            if self.Emax < np.max(eig_energ):
                self.Emax = np.max(eig_energ)

        print 'precalc took: ', time.time()-t1, ' s'
        # arg for part of time evo

    def do_definitions(self):
        # redo hamilton
        self.b = 1./2 * (self.Emax+self.Emin)
        self.epsilon = self.alpha * (self.Emax-self.Emin)
        self.a = 1./2 * (self.Emax-self.Emin + self.epsilon)

        self.coef_arg = self.a*self.delta_t

#        self.Hcheb = lambda t: 2* (self.H(t)-self.Emin*np.diag(np.ones(self.dim)))/(self.Emax-self.Emin) - np.diag(np.ones(self.dim))
        self.Hcheb = lambda t: (self.H(t)-self.b*np.diag(np.ones(self.dim)) )/self.a

    def do_timestep(self, psiin, t):

        # reserve mem
        psi_cheby = -1*np.zeros((self.chebyshev_order,self.dim), dtype='complex128')
        a_cheby = -1*np.zeros(self.chebyshev_order, dtype='complex128')

        Hnow = self.Hcheb(t)

        psi_cheby[0,:] = psiin
        psi_cheby[1,:] = np.dot(Hnow, psi_cheby[0,:])

        a_cheby[0] = spsp.jn(0, self.coef_arg)
        tmp1 = np.arange(1,self.chebyshev_order,1)
        a_cheby[1:] = 2*(-1j)**(tmp1) * spsp.jn(tmp1,self.coef_arg)

        for k in xrange(self.chebyshev_order-2):
            psi_cheby[k+2,:] = 2*np.dot(Hnow, psi_cheby[k+1,:]) - psi_cheby[k,:]


        # sum of coeffiecents times the state
        psiout = np.exp(-1j*(self.Emax+self.Emin)*self.delta_t/2) * np.tensordot(a_cheby, psi_cheby,(0,0))

        return psiout

    def do_timeevo(self):
        for k in xrange(self.steps-1):
            self.y[:,k+1] = self.do_timestep(self.y[:,k], self.time[k])
            self.pop[:,k+1] = np.abs(self.y[:,k+1])**2



if __name__ == "__main__":

    # ------------------
    # settings

    doplot = True

    delta_t = 0.01
    T0 = 0
    Tmax = 1

    T = np.arange(delta_t,Tmax,delta_t)
    Y = -1* np.ones((len(T),4), dtype='complex128')
    Pop = -1* np.ones((len(T),4))

    # ------------------
    # interaction/hamiltonian

    rabi = 2*np.pi*2
    beatfreq = 2*np.pi*4
    id = np.diag(np.array([1,1]))
    sigmax = 1./2*np.array([[0,1],[1,0]])
    Hint = rabi*npml.kron(id,sigmax)
    H0 = np.diag(np.array([0,0,1,1])) * 100
    H = lambda t: (1+np.cos(beatfreq * t))/2 * Hint + H0

    # input state

    psiin = np.array([0,1,0,0])
    #psiin = 1./np.sqrt(2)* np.array([1,1])
    tcur = T0

    b = ODE_timeevo(H, T0, Tmax, delta_t, psiin)

    t2 = time.time()
    a = Chebyshev_timeevo(H, T0, Tmax, delta_t, psiin)
    print 'chebyshev evo: ',time.time()-t2


    # ------------------
    # plot result

    if doplot:
        fig1 = pl.figure(1)
        fig1.clf()
        ax1 = fig1.add_subplot(111)
        pl.hold(True)
        ax1.plot(b.time,b.pop[0,:], label = 'pop 0, ode', color = 'red')
        ax1.plot(b.time,b.pop[1,:], label = 'pop 1, ode', color = 'green')

        ax1.plot(a.time,a.pop[0,:], label = 'pop 0, cheby', color='blue', linestyle = 'dashed')
        ax1.plot(a.time,a.pop[1,:], label = 'pop 1, cheby', color='black', linestyle = 'dashed')


        pl.hold(False)
        pl.legend(loc = 'upper right')
        pl.axis([0, Tmax,0,1])
        fig1.show()



# SchroedingerEqSolvers.py ends here
