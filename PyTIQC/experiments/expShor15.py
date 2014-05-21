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

def main():
    ### run params ##################################
    pulseseqfileShor = 'experiments/Shorseq.py'
   
    firstMultMap = True # Do first multiplication with CNOT 
                        # mapping instead of Multiplier 
                        # (does nothing in easy cases)
    select_a = 7 # Factoring base: 2, 4, 7, 8, 13, 11, or 14

    doRun = True
    doIdeal = True
    printpulse = False
    saveKitaev = False
    doPPKitaev = True

    print 'a =',select_a
    
    NumberOfIons = 5
    NumberOfPhonons = 0 if doIdeal else 7 
    hspace = sim.hspace(NumberOfIons,2,NumberOfPhonons,0)
    hspace.initial_state("quantum", qstate='DSSSS')
    params = sim.parameters(hspace)
    
    params.use_servers( ['all'] )
    params.ppcpus = 16
    params.shortestMS = 16
    params.calcPulseParams()
    params.progbar = True
    params.saveallpoints = False
    params.coherenceTime = 15000
    params.hidingerr = 1
    
    params.printpulse = printpulse
    params.progbar = printpulse
    params.pplog = False
    
    dec = sim.decoherence(params)
    dec.doRandNtimes = 16
    dec.doPP = True
    dec.dict['all'] = True
#    if doPPKitaev: params.ppcpus = 2
    
    Kit = Kitaev.Kitaev()
    
    ##########################################
    # load the pulse sequences
    # change ion order and define permutations
    ##########################################
    execfile(pulseseqfileShor,locals(),globals())

    Fred12 = copy.deepcopy(Fredkin)
    Fred23 = copy.deepcopy(Fredkin)
    Fred34 = copy.deepcopy(Fredkin)
    Fred24 = copy.deepcopy(Fredkin)
    Fred13 = copy.deepcopy(Fredkin)
    Fred12.changeions(params, (1,2,0))
    Fred23.changeions(params, (2,3,0))
    Fred34.changeions(params, (3,4,0))
    Fred24.changeions(params, (2,4,0))
    Fred13.changeions(params, (1,3,0))
   
    cnot13 = copy.deepcopy(cnot12)
    cnot24 = copy.deepcopy(cnot12)
    cnot13.changeions(params,(0,1,3))
    cnot24.changeions(params,(0,2,4))

    times2 = sim.PulseSequence([ \
        sim.Hide(params, 3, True),
        sim.Hide(params, 4, True), #vis: 1,2
        Fred12,
        sim.Hide(params, 3, False),
        sim.Hide(params, 1, True), #vis: 2,3
        Fred23,
        sim.Hide(params, 4, False),
        sim.Hide(params, 2, True), #vis: 3,4
        Fred34,
        sim.Hide(params, 1, False),
        sim.Hide(params, 2, False) #vis: 1,2,3,4
        ])

    times4 = sim.PulseSequence([ \
        sim.Hide(params, 1, True),
        sim.Hide(params, 3, True), #vis: 2,4
        Fred24,
        sim.Hide(params, 1, False),
        sim.Hide(params, 3, False),
        sim.Hide(params, 2, True),
        sim.Hide(params, 4, True), #vis: 1,3
        Fred13,
        sim.Hide(params, 2, False),
        sim.Hide(params, 4, False)
        ])

    times8 = sim.PulseSequence([ \
        sim.Hide(params, 1, True),
        sim.Hide(params, 2, True), #vis: 3,4
        Fred34,
        sim.Hide(params, 2, False),
        sim.Hide(params, 4, True), #vis: 2,3
        Fred23,
        sim.Hide(params, 1, False),
        sim.Hide(params, 3, True), #vis: 1,2
        Fred12,
        sim.Hide(params, 3, False),
        sim.Hide(params, 4, False)
        ])

    times13 = copy.deepcopy(times8)
    times13.append(cnot1234)

    times7 = copy.deepcopy(times2)
    times7.append(cnot1234)


    ### general pulse sequence
    def GeneratePulseSeq(a):
    
        NOP = sim.PulseSequence([sim.Delay(params, 0.1)])

        if a in (2,7,8,13): # hard cases
            if firstMultMap:
                a2modN = cnot24
            else:
                a2modN = times4
        else:               # easy cases
            a2modN = NOP
        
        if   a == 2:
            amodN = times2
        elif a == 7:
            amodN = times7
        elif a == 8:
            amodN = times8
        elif a == 13:
            amodN = times13
        elif a == 4:
            amodN = cnot24
        elif a == 11:
            amodN = cnot13
        elif a == 14:
            amodN = cnot1234

        pulseseq_group = Kit.GeneratePulseSeq(params, [NOP, a2modN, amodN])
        return pulseseq_group
    
    #######################################################
    
    pulseseq = GeneratePulseSeq(select_a)
    
    if doIdeal:
        for ps in pulseseq:
            for pulse in ps:
                pulse.use_ideal = True
    
    ### run it
    if doRun:
        tic = time.time()
        result = Kit.simulateevolution(pulseseq, params, dec, doPP=doPPKitaev)
        if not doIdeal: print "runtime: ", time.time()-tic, "sec"
        print np.around(result,6)
    
    if saveKitaev:
        timestr = datetime.datetime.now().strftime('%Y%m%d-%H%M%S-%f')
        fnameappend = ''
        Kit.getQFTall()
        Kit.datagroupsave('Shor'+str(select_a)+'-data-'+fnameappend+timestr+'.pk', 
                          RhoOnly=True)
        Kit.resultSave('Shor'+str(select_a)+'-res-'+fnameappend+timestr+'.npy')
        f = open('Shor'+str(select_a)+'-params-'+fnameappend+timestr+'.pk', 'wb')
        pickle.dump(dec.dict_params_static, f)
        f.close()

    return True#dataobj

if __name__ == '__main__':
    main()
