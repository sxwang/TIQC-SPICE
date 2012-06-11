''' Order-finding experiment. 
3 or 5 ions, using 4 pulse sequences for 4 permutation functions '''

import numpy as np
import PyTIQC.core.simtools as sim
import PyTIQC.core.qctools as qc
import PyTIQC.core.gates as U
import matplotlib.pyplot as pl

import time, sys, copy

import Kitaev 
reload(Kitaev)

pi = np.pi

order3ideal = np.array([0.3438, 0.0145, 0.0625 , 0.2355 ,0.0312, 0.2355, 0.0625, 0.0145])

### run params ##################################
pulseseqfileOF = 'experiments/Orderfindseq.py'
pulseseqfileQFT = 'experiments/QFTseq.py'

doRun = True
doIdeal = True
doAll = False
printpulse = False
select_order = 3
doKitaev = True
saveKitaev = False

NumberOfIons = 5 if not doKitaev else 3
NumberOfPhonons = 0 if doIdeal else 7
hspace = sim.hspace(NumberOfIons,2,NumberOfPhonons,0)
if select_order == 1 or select_order == 2:
    hspace.initial_state("quantum", qstate='DDS' if doKitaev else 'DDSSS')
else:
    hspace.initial_state("quantum", qstate='SSS' if doKitaev else 'SSSSS')
params = sim.parameters(hspace)

params.use_servers( ['local'] )
params.shortestMS = 16
params.calcPulseParams()
params.progbar = False
params.saveallpoints = False
params.coherenceTime = 25000

params.printpulse = printpulse
params.progbar = printpulse

dec = sim.decoherence(params)
dec.doRandNtimes = 1
dec.dict['all'] = True
#dec.doPP = True

##########################################
# load the pulse sequences
# change ion order and define permutations
##########################################
execfile(pulseseqfileOF)
execfile(pulseseqfileQFT)
if not doKitaev:
    order1.changeions(params, (4,3,2))
    order2.changeions(params, (4,3,2))
    order31 = copy.deepcopy(order3)
    order34 = copy.deepcopy(order3)
    order31.changeions(params, (4,3,2))
    order32.changeions(params, (4,3,1))
    order34.changeions(params, (4,3,0))
    order4.changeions(params, (4,3,2))
    order42.changeions(params, (4,3,1))
    QFTseq.changeions(params, (2,1,0))
else:
    order1.changeions(params, (2,1,0))
    order2.changeions(params, (2,1,0))
    order31 = copy.deepcopy(order3)
    order34 = copy.deepcopy(order3)
    order31.changeions(params, (2,1,0))
    order32.changeions(params, (2,1,0))
    order34.changeions(params, (2,1,0))
    order4.changeions(params, (2,1,0))
    order42.changeions(params, (2,1,0))
    
Kit = Kitaev.Kitaev()

def SelectPermutation(select_order, doKitaev):
    if select_order == 0:
        Perm = sim.Delay(params, 0.1)
        Perm2 = sim.Delay(params, 0.1)
        Perm4 = sim.Delay(params, 0.1)
    elif select_order == 1:
        Perm = order1
        Perm2 = sim.Delay(params, 0.1)
        Perm4 = sim.Delay(params, 0.1)
    elif select_order == 2:
        Perm = order2
        Perm2 = sim.Delay(params, 0.1)
        Perm4 = sim.Delay(params, 0.1)
    elif select_order == 3:
        Perm = order31
        if not doKitaev:
            Perm2 = sim.PulseSequence([ \
                sim.Hide(params, 1, False),
                sim.Hide(params, 2, True),
                order32 ])
            Perm4 = sim.PulseSequence([ \
                sim.Hide(params, 2, True),
                sim.Hide(params, 1, True),
                sim.Hide(params, 0, False),
                order34, ])
        else:
            Perm2 = order32
            Perm4 = order34
    elif select_order == 4:
        Perm = order4
        if not doKitaev:
            Perm2 = sim.PulseSequence([ \
                sim.Hide(params, 1, False),
                sim.Hide(params, 2, True),
                order42 ])
        else:
            Perm2 = order42
        Perm4 = sim.Delay(params, 0.1)

    return Perm, Perm2, Perm4

### general pulse sequence
def GeneratePulseSeq(Perm, Perm2, Perm4, QFTseq, doKitaev):
    if not doKitaev:
        pulseseq = sim.PulseSequence([ \
            sim.Hide(params, 3, True),
            sim.Hide(params, 4, True),
            sim.Rcar(params, pi/2, -pi/2),
            sim.Hide(params, 3, False),
            sim.Hide(params, 4, False),
            sim.Hide(params, 0, True),
            sim.Hide(params, 1, True,),
            Perm, 
            Perm2,
            Perm4,
            sim.Hide(params, 0, False),
            sim.Hide(params, 1, False),
            sim.Hide(params, 2, False),
            sim.Hide(params, 3, True),
            sim.Hide(params, 4, True),
            QFTseq,
        ])
        return pulseseq
    else:
        pulseseq_group = Kit.GeneratePulseSeq(params, Perm4, Perm2, Perm)

        return pulseseq_group

#######################################################

[Perm, Perm2, Perm4] = SelectPermutation(select_order, doKitaev)
pulseseq = GeneratePulseSeq(Perm, Perm2, Perm4, QFTseq, doKitaev)

if not doRun: sys.exit()

if doIdeal:
    if not doKitaev:
        for pulse in pulseseq:
            pulse.use_ideal = True
    else:
        for ps in pulseseq:
            for pulse in ps:
                pulse.use_ideal = True


### run it
if doAll:
    qstates_list = ('DDSSS', 'SDSSS', 'DSSSS', 'SSSSS') if not doKitaev \
        else ('DDS', 'SDS', 'DSS', 'SSS')
    orders_list = (1, 2, 3, 4)
#    for y0, select_order in zip(qstates_list, orders_list):
    for select_order in orders_list:
      for y0 in qstates_list:
        [Perm, Perm2, Perm4] = SelectPermutation(select_order, doKitaev)
        pulseseq = GeneratePulseSeq(Perm, Perm2, Perm4, QFTseq, doKitaev)
        if doIdeal:
            if not doKitaev:
                for pulse in pulseseq:
                    pulse.use_ideal = True
            else:
                for ps in pulseseq:
                    for pulse in ps:
                        pulse.use_ideal = True

        print "initial state", y0, "| expected order", select_order
        params.initial_state("quantum", qstate = y0)

        tic = time.time()
        if doKitaev:
            result = Kit.simulateevolution(pulseseq, params, dec)
        else:
            data = qc.simulateevolution(pulseseq, params, dec)
            result = U.displaytracedstates(data.YRPN[-1], pop=True)
        if not doIdeal: print "runtime: ", time.time()-tic, "sec"
        

else:
    tic = time.time()
    if doKitaev:
        result = Kit.simulateevolution(pulseseq, params, dec)
        print np.around(result, 3)

    else:
        data = qc.simulateevolution(pulseseq, params, dec)
        result = U.displaytracedstates(data.YRPN[-1], pop=True)
        print np.around(result, 3)
    if not doIdeal: print "runtime: ", time.time()-tic, "sec"

if saveKitaev:
    qc.saveRun(pulseseq, params, dec, Kit.data_group, 'OF'+str(select_order)+'-werr.shlv')
    Kit.getQFTall()
    Kit.resultSave('OF'+str(select_order)+'-result.npy')
