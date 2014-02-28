''' simulation for the 5-qubit (Kitaev) Shor experiment '''

import numpy as np
import PyTIQC.core.simtools as sim
import PyTIQC.core.qctools as qc
import PyTIQC.core.gates as U
import matplotlib.pyplot as pl

import time, sys, pickle, datetime
import shelve, copy

import Kitaev
reload(Kitaev)

pi = np.pi

### run params ##################################
pulseseqfileShor = 'experiments/Shorseq.py'

doRun = True
doIdeal = True
printpulse = False
select_a = 7 # 7 or 11
saveKitaev = False
doPPKitaev = True

NumberOfIons = 5
NumberOfPhonons = 0 if doIdeal else 7
hspace = sim.hspace(NumberOfIons,2,NumberOfPhonons,0)
hspace.initial_state("quantum", qstate='DSSSS')
params = sim.parameters(hspace)

params.use_servers( ['local'] )
params.shortestMS = 16
params.calcPulseParams()
params.progbar = False
params.saveallpoints = False
params.coherenceTime = 20000
params.hidingerr = 1

params.printpulse = printpulse
params.progbar = printpulse
params.pplog = False

dec = sim.decoherence(params)
dec.doRandNtimes = 1
dec.doPP = True

if doPPKitaev: params.ppcpus = 2

Kit = Kitaev.Kitaev()

##########################################
# load the pulse sequences
# change ion order and define permutations
##########################################
execfile(pulseseqfileShor)

shor7b1 = copy.deepcopy(Fredkin)
shor7b2 = copy.deepcopy(Fredkin)
shor7b1.changeions(params, (0,1,3))
shor7b2.changeions(params, (0,2,4))
shor7b = sim.PulseSequence([ \
    sim.Hide(params, 2, True),
    sim.Hide(params, 4, True),
    shor7b1,
    sim.Hide(params, 2, False),
    sim.Hide(params, 4, False),
    sim.Hide(params, 1, True),
    sim.Hide(params, 3, True),
    shor7b2,
    sim.Hide(params, 1, False),
    sim.Hide(params, 3, False)])

    
### general pulse sequence
def GeneratePulseSeq(shor7a, shor7b, shor11, a):

    NOP = sim.Delay(params, 0.1)
    if a == 7:
        a2modN = shor7b
        amodN = shor7a
    elif a == 11:
        a2modN = sim.Delay(params, 0.1)
        amodN = shor11

    pulseseq_group = Kit.GeneratePulseSeq(params, NOP, a2modN, amodN)
    return pulseseq_group

#######################################################

pulseseq = GeneratePulseSeq(shor7a, shor7b, shor11, select_a)

if doIdeal:
    for ps in pulseseq:
        for pulse in ps:
            pulse.use_ideal = True

### run it
if doRun:
    tic = time.time()
    result = Kit.simulateevolution(pulseseq, params, dec, doPP=doPPKitaev)
    if not doIdeal: print "runtime: ", time.time()-tic, "sec"
    print np.around(result,3)

if saveKitaev:
    timestr = datetime.datetime.now().strftime('%Y%m%d-%H%M%S-%f')
    Kit.getQFTall()
    Kit.datagroupsave('Shor'+str(select_a)+'-data-'+fnameappend+timestr+'.pk', 
                      RhoOnly=True)
    Kit.resultSave('Shor'+str(select_a)+'-res-'+fnameappend+timestr+'.npy')
    f = open('Shor'+str(select_a)+'-params-'+fnameappend+timestr+'.pk', 'wb')
    pickle.dump(dec.dict_params_static, f)
    f.close()
