#!/usr/bin/env python
# -*- mode: Python; coding: latin-1 -*-
# Time-stamp: "2011-06-22 09:36:04 c704252"

#  file       densitymatrix.py
#  copyright  (c) Thomas Monz 2009
#  url        http://wiki.havens.de

from numpy import *
#from scipy.linalg import *
import numpy as np
import PyTIQC.tools.quantum_tools as qtls
import scipy.linalg as splnlg

class DensityMatrixObject:
  def __init__(self,exp_rho):
    # check whether it's just a state
    if (1 in exp_rho.shape) or (len(exp_rho.shape)==1):
      print('State provided, fixed into density matrix')
      exp_rho=np.outer(exp_rho,conjugate(exp_rho))

    # absolute measures

    self.exprho=exp_rho
    self.normalized=trace(self.exprho)
    self.purity = np.real(np.trace(np.dot(self.exprho,self.exprho)))
    self.linentropy = 1-float(real(trace(dot(self.exprho,self.exprho))))
    # easier to calculate over eigenvalues
#    self.vNeumannEntropy = -trace(dot(self.exprho,logm(self.exprho)))

    self.idealrho= None
    self.fidelity= None
    self.jozsafidelity= None
    self.tracedistance= None
    self.relativeentropy = None

  # relative measures

  def fid(self,idealrho=None):
    self.check_idealrho(idealrho)
    self.fidelity=trace(dot(self.idealrho,self.exprho))
    return self.fidelity

  def jozsafid(self,idealrho=None):
    self.check_idealrho(idealrho)
    tmp=qtls.sqrtm_dm(self.exprho)
    self.jozsafidelity=real(trace(qtls.sqrtm_dm(dot(dot(tmp,self.idealrho),tmp)))**2).tolist()
    return self.jozsafidelity
#    print('JozsaFidelity: %.5f' % self.jozsafidelity)

  def trdistance(self,idealrho=None):
    self.check_idealrho(idealrho)
    tmp=self.idealrho-self.exprho
    self.tracedistance=np.real(np.trace(qtls.sqrtm_dm(np.dot(tmp.conjugate().transpose(),tmp))))/2
    return self.tracedistance
#    print('TraceDistance: %.5f' % self.tracedistance)

  def sso(self,idealrho=None):
    self.check_idealrho(idealrho)
    p1 = np.diag(self.idealrho)
    p2 = np.diag(self.exprho)
    return np.real(np.sum(np.sqrt(p1) * np.sqrt(p2))**2)

  def trdistancePop(self, idealrho=None):
    self.check_idealrho(idealrho)
    p1 = np.diag(self.idealrho)
    p2 = np.diag(self.exprho)
    return np.sum(np.abs(p1-p2))

  def relentropy(self,idealrho=None):
    # attention - log not working on pure matrices with 0-eigenvalues
    self.check_idealrho(idealrho)
    a = np.trace(np.dot(self.idealrho,splnlg.logm(self.idealrho)))
    b = np.trace(np.dot(self.idealrho,splnlg.logm(self.exprho)))
    self.relativeentropy=np.real(a-b)
    return self.relativeentropy
#    print('RelativeEntropy: %.5f' % self.relativeentropy)

  def fixsinglequbitphases(self):
    print('Trying to optimise single qubit phases')

  def tracequbits(self,list_ions):
    print('Tracing over the specified qubits')

  def check_idealrho(self,newrho=None):
    if self.idealrho==None and newrho==None:
      print('No ideal density matrix provided, stopping here')
    else:
      if (not self.idealrho==None) and (not newrho==None):
#        print('Replacing old ideal matrix with new one')
        self.idealrho=newrho
      if self.idealrho==None:
#        print('Using provided matrix as ideal matrix')
        self.idealrho=newrho

  def check_GHZ_entanglement(self):
    rho_size = self.exprho.shape[0]
    off_diag = abs(self.exprho[rho_size-1,0])
    sum_array = np.zeros(rho_size/2.0-1)
    for index in xrange(rho_size/2.0-1):
      sum_array[index] = self.exprho[index+1,index+1] + \
          self.exprho[rho_size-index-1][rho_size-index-1]
    print "max unwanted: " +str(np.max(sum_array))
    print "off_diag: " + str(off_diag)



if __name__ == "__main__":

  a=1/sqrt(2)*array([1,0,0,1])
  idrho=outer(a,conjugate(a))
  tmp=idrho.copy()
  tmp[0,3]*=0.95
  tmp[1,1]=0.01
  tmp[2,2]=0.02
  tmp[3,0]*=0.95
  tmp[3,3]=0.47
#  print(tmp)
  exrho=tmp


  tst=DensityMatrixObject(exrho)
  tst.fid(idrho)
  print tst.jozsafid()

  tst=DensityMatrixObject(idrho)
  tst.fid(exrho)
  print tst.jozsafid()


  tst.trdistance()
#  tst.relentropy()

# densitymatrix.py ends here

