"""IterML
module for iterative maximum likelihood reconstruction

Example for a iterative state tomo

>>> import qpython.densitmatrixreconstruction as dmr
>>> data = rd.ReadData('1607')
>>> rho=dmr.IterML.iterfun(data, 100)

or

>>> import qpython.densitmatrixreconstruction as dmr
>>> rho=dmr.IterML.iterfun('1607', 100)

"""
# psyco magic?
# import psyco
# psyco.full()

import numpy as np
import numpy.matlib as npml
import time
try:
    import PyTIQC.evaluation.readdata as rd
except:
    pass
import densitymatrix
import scipy.io

if id(np.dot) == id(np.core.multiarray.dot):
    print " Not using blas/lapack!"
    print " switching to blas/lapack!"
    np.alterdot()

# ------------------------------------

def __iter_readdata(filename):
    data=np.loadtxt(filename,skiprows=1)
    return data

# ------------------------------------

def iterfun(data,NumberOfIterations, path=None):
    """iterative maximum likelihood state tomography
    data can be either a matrix,data_object or a time string
    """
    if type(data)== str:
        data = rd.ReadData(data, path=path)
    try:
        data = data.data_dict['cprb']
    except AttributeError:
        pass
#    tic=time.time()
    NumberRows=len(data[:,1])
    NumberCols=len(data[1,:])
    # check if number of rows and columns fit
    # numberofcols= 2^NumberOfIons+1
    # numberrows=3^numberofions
    if not 3**np.log2(NumberCols-1)==NumberRows:
        print("dataset not complete")
    NumberOfIons=np.int(np.log2(NumberCols-1))
    probs=data[:,1:]
    pulsecode=data[:,0]
    NoExp=len(pulsecode)
    # pulses = double(dec2base(pulsecode,3,NoI))-double('0');
    pulses=np.zeros((3**NumberOfIons,NumberOfIons))
#  pulses.shape=
    for i in xrange(3**NumberOfIons):
        pulses[i,:]=np.array(__dec2base(i,3,NumberOfIons))
    #print(pulses)
    # pulsetable = reshape(kron(ones(1,2^NoI),pulses)',NoI,NoExp*2^NoI)';
    # first part kron(ones(1,2^NoI),pulses)'  =
    #  = mtlb.transpose(mtlb.kron(ones(2**NumberOfIons),pulses))
    # reshape = reshape(<whatever>,(NumberOfIons,3**NumberOfIons*2**NumberOfIons),'FORTRAN')
    a=npml.kron(np.ones(2**NumberOfIons),pulses).transpose()
    pulsetable = np.reshape(a,(NumberOfIons,3**NumberOfIons*2**NumberOfIons),'FORTRAN').transpose()
    #print(pulsetable)
    # Now the experimental data are stored in the same way:
    # probs = reshape(probs',1,NoExp*2^NoI)';
    probs=np.reshape(probs.transpose(),(1,3**NumberOfIons*2**NumberOfIons),'FORTRAN').transpose()
    probs=probs[:,0]
    # For each experimental data point, a measurement of the D state has to be
    # labeled with +1, a measurement of the S-state with 0;
    # meas = (double(dec2bin(2^NoI-1:-1:0)')-double('0'))';
    # meastable = kron(ones(NoExp,1),meas);
    meas=np.zeros((2**NumberOfIons,NumberOfIons))
    k=0
    for i in xrange(2**NumberOfIons-1,-1,-1):
        meas[k,:]=np.array(__dec2base(i,2,NumberOfIons),dtype=int)
        k+=1
    a=np.ones((3**NumberOfIons,1))
    meastable=npml.kron(a,meas)
    #meastable = meastable + 2 * pulsetable + 1;
    #Ntable=length(meastable);
    meastable+= 2*pulsetable
    Ntable=len(meastable)

    #Here are the corresponding projectors:
    #P(:,:,1) = [0;1]*[0;1]';               % -
    #P(:,:,2) = [1;0]*[1;0]';               % +
    #P(:,:,4) = [-1;1]*[-1;1]'/2;           % -
    #P(:,:,3) = [1;1]*[1;1]'/2;             % +
    #P(:,:,6) = [i;1]*[i;1]'/2;             % -
    #P(:,:,5) = [-i;1]*[-i;1]'/2;           % +
    P=np.zeros((6,2,2),dtype=complex)
    P[0,:,:]=np.outer(np.array([0,1]),np.array([0,1]))
    P[1,:,:]=np.outer(np.array([1,0]),np.array([1,0]))
    P[2,:,:]=np.outer(0.5*np.array([-1,1]),np.array([-1,1]))
    P[3,:,:]=np.outer(0.5*np.array([1,1]),np.array([1,1]))
    P[4,:,:]=np.outer(0.5*np.array([1j,1]),np.conjugate(np.array([1j,1])))
    P[5,:,:]=np.outer(0.5*np.array([-1j,1]),np.conjugate(np.array([-1j,1])))
    # about to start iterations ...
    rho = np.identity(2**NumberOfIons)/2**NumberOfIons
    # AllOp=zeros(2^NoI,2^NoI,Ntable);
#  AllOp=zeros((Ntable,2**NumberOfIons,2**NumberOfIons),dtype=complex)
    AllOp=[]
#    toc2=time.time()
#    print toc2-tic,'seconds for initialisation have elapsed'
    AllOp2 = np.zeros((Ntable,2**NumberOfIons,2**NumberOfIons),dtype=complex)
    AllOpTransposed = np.zeros((Ntable,2**NumberOfIons,2**NumberOfIons),dtype=complex)
    for k in xrange(Ntable):
        ind=meastable[k,:].copy()
        Op=P[ind[0],:,:]
        for m in xrange(1,NumberOfIons):
            Op=npml.kron(Op,P[ind[m],:,:])
        AllOp.append(Op)
        AllOp2[k,:,:] = Op
        AllOpTransposed[k,:,:] = Op.transpose()
    # really starting with the iterations now
#    toc3=time.time()
#    print toc3-toc2,'seconds for operator initalisation'
    ROp_start=np.zeros((2**NumberOfIons,2**NumberOfIons),dtype=complex)
    list_probOp = np.zeros(Ntable)
    for i in xrange(NumberOfIterations):
        ROp=ROp_start.copy()
        list_probOp2 = np.sum(np.sum(npml.multiply(rho, AllOpTransposed), axis = 1), axis = 1)
        # for k in xrange(Ntable):
        #     Op=AllOp[k]
        #     # the following is tons faster because it relies on element wise multiplication only
        #     probOp=(rho*Op.transpose()).sum()
        #     list_probOp[k] = probOp
            # okay. if probOp would be zero, it would be a problem. but i
            # never got a zero here, so i'll just skip the if
            # if probOp > 0:
            #     ROp += probs[k]/probOp * Op
        # tensordot results in a factor of 2 faster evaluation of
        # the w4 data compared to the old loop approach
        # print (list_probOp2[0] -list_probOp[0])

        ROp2 = np.tensordot(probs/list_probOp2,AllOp2,(0,0))
        rho=np.dot(np.dot(ROp2,rho),ROp2)
#        rho=np.dot(np.dot(ROp,rho),ROp)
        rho/=rho.trace()
#    toc=time.time()
#    print time.time()-toc3,' seconds for iterations'
    return rho


# ------------------------------------

def iterfun_obj(otherself, data, NumberOfIterations):
    rho = iterfun(data, NumberOfIterations)
    a1 = densitymatrix.DensityMatrixObject(rho)
    return a1

# ------------------------------------


def __dec2base_precise(n, b):
    (n1,n2) = divmod(n, b)
    if n1 == 0:
        return [n2]
    else:
        return __dec2base_precise(n1, b) + [n2]

def __dec2base(n,b,l):
    list=__dec2base_precise(n,b)
    for i in xrange(l-len(list)):
        list.insert(0,0)
    return list

# ------------------------------------

if __name__ == "__main__":

    w2data=__iter_readdata('w2_global.dat')
    # rho_w2=iterfun(w2data,100)
    # rho_w2_ideal=array([[0,0,0,0],[0,0.5,0.5,0],[0,0.5,0.5,0],[0,0,0,0]])
    # print 'overlap with w2: ', real(trace(dot(rho_w2_ideal,absolute(rho_w2))))

    w4data=__iter_readdata('w4_global.dat')
    rho_w4=iterfun(w4data,100)
    ideal_psi_w4=1/np.sqrt(4)*np.array([0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]);
    rho_w4_ideal=np.outer(ideal_psi_w4,ideal_psi_w4)
    overlap = np.real(np.trace(np.dot(rho_w4_ideal,np.absolute(rho_w4))))
    if overlap > 0.75:
        print 'passed'

    #print 'overlap with w4: ', overlap

    # w6data=__iter_readdata('w6_addressed.dat')
    # rho_w6=iterfun(w6data,100)

# IterML.py ends here
