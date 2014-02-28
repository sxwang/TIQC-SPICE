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

def go():
    ### run params ##################################
    pulseseqfileShor = 'experiments/Shorseq2.py'
   
    OLD = False # Do pulse sequences from experiment (does nothing if a=11)
    HARD = True # Do first multiplication with CNOTs or Fredkins (HARD=Fredkins) (does nothing if a = 11)
    select_a = 7 # 2, 7, 8, 13, or 11

    doRun = True
    doIdeal = False
    printpulse = False
    saveKitaev = False
    doPPKitaev = True

    print 'a =',select_a,'OLD =',OLD,'HARD =',HARD
    
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
#    execfile(pulseseqfileShor)

    Fredkin = sim.PulseSequence([ \
    sim.Rcar(params, 0.5*pi, 0.5*pi),
    sim.Rac(params, 0.5*pi, 0, 0),
    sim.RMS(params, 0.5*pi, pi/2),
    sim.Rac(params, 0.5*pi, 0, 1),
    sim.Rac(params, -0.5*pi, 0, 0),
    sim.Rcar(params, 0.75*pi, 1*pi),
    sim.RMS(params, -0.25*pi, pi/2),
    sim.Rac(params, 0.5*pi, 0, 1),
    sim.RMS(params, 0.5*pi, pi/2),
    sim.Rcar(params, 0.5*pi, 0*pi),
    sim.Rac(params, -0.25*pi, 0, 1),
    sim.Rac(params, 0.5*pi, 0, 0),
    sim.RMS(params, 0.5*pi, pi/2),
    sim.Rac(params, 0.5*pi, 0, 1),
    sim.Rac(params, 0.5*pi, 0, 2),
    sim.Rcar(params, 0.5*pi, 0*pi),
    sim.Rac(params, -0.5*pi, 0, 2),
    sim.Rac(params, -0.5*pi, 0, 1)  ])
    
    cnot1234 = sim.PulseSequence([ \
    sim.Rcar(params, 1.5*pi, 0*pi),
    sim.Rac(params, 1.5*pi, 0, 0),
    sim.RMS(params, 0.25*pi, 0.5*pi),
    sim.Rcar(params, 0.75*pi, 0*pi),
    sim.Rac(params, 1*pi, 0, 0),
    sim.RMS(params, 1.75*pi, 0.5*pi),
    sim.Rcar(params, 0.75*pi, 0*pi),
    sim.Rac(params, 1.5*pi, 0, 0),
    sim.Rcar(params, 0.5*pi, 0*pi) ])

    cnot1234_old = sim.PulseSequence([ \
    sim.Rcar(params, 0.6959*pi, 1.25*pi),
    sim.Rac(params, 0.6667*pi, 0, 0),
    sim.Rcar(params, 0.6959*pi, 0.25*pi),
    sim.RMS(params, 0.75*pi, 0.5*pi),
    sim.Rac(params, 1.0*pi, 0, 0),
    sim.RMS(params, 1.25*pi, 0.5*pi),
    sim.Rac(params, 1.5*pi, 0, 0),
    sim.Rcar(params, 0.75*pi, 0.0*pi),
    sim.Rcar(params, 0.75*pi, 0.0*pi)
    ])
    
    
    shor11 = sim.PulseSequence([ \
    sim.Rcar(params, pi/2, 0),
    sim.Rac(params, -pi/2, 0, 0),
    sim.Rcar(params, pi/4, 0),
    sim.RMS(params, pi/8, pi/2),
    sim.Rac(params, -pi, 0, 1),
    sim.Rac(params, -pi, 0, 3),
    sim.Rcar(params, pi/4, 0),
    sim.RMS(params, -pi/8, pi/2),
    sim.Rac(params, -pi, 0, 0),
    sim.RMS(params, pi/8, pi/2),
    sim.Rac(params, -pi, 0, 1),
    sim.Rac(params, -pi, 0, 3),
    sim.RMS(params, -pi/8, pi/2),
    sim.Rcar(params, pi/2, pi),
    sim.Rac(params, -pi/2, 0, 0),
    sim.Rcar(params, pi/2, pi) ])
    
    
    cnot13 = sim.PulseSequence([ \
    sim.Rcar(params, pi/2, 0),
    sim.Rac(params, -pi/2, 0, 0),
    sim.Rcar(params, pi/4, 0),
    sim.RMS(params, pi/8, pi/2),
    sim.Rac(params, -pi, 0, 1),
    sim.Rac(params, -pi, 0, 3),
    sim.Rcar(params, pi/4, 0),
    sim.RMS(params, -pi/8, pi/2),
    sim.Rac(params, -pi, 0, 0),
    sim.RMS(params, pi/8, pi/2),
    sim.Rac(params, -pi, 0, 1),
    sim.Rac(params, -pi, 0, 3),
    sim.RMS(params, -pi/8, pi/2),
    sim.Rcar(params, pi/2, pi),
    sim.Rac(params, -pi/2, 0, 0),
    sim.Rcar(params, pi/2, pi) ])
    
    cnot24 = sim.PulseSequence([ \
    sim.Rcar(params, pi/2, 0),
    sim.Rac(params, -pi/2, 0, 0),
    sim.Rcar(params, pi/4, 0),
    sim.RMS(params, pi/8, pi/2),
    sim.Rac(params, -pi, 0, 2),
    sim.Rac(params, -pi, 0, 4),
    sim.Rcar(params, pi/4, 0),
    sim.RMS(params, -pi/8, pi/2),
    sim.Rac(params, -pi, 0, 0),
    sim.RMS(params, pi/8, pi/2),
    sim.Rac(params, -pi, 0, 2),
    sim.Rac(params, -pi, 0, 4),
    sim.RMS(params, -pi/8, pi/2),
    sim.Rcar(params, pi/2, pi),
    sim.Rac(params, -pi/2, 0, 0),
    sim.Rcar(params, pi/2, pi) ])
    
    Fred12 = copy.deepcopy(Fredkin)
    Fred23 = copy.deepcopy(Fredkin)
    Fred34 = copy.deepcopy(Fredkin)
    Fred12.changeions(params, (1,2,0))
    Fred23.changeions(params, (2,3,0))
    Fred34.changeions(params, (3,4,0))
    
    shor7bcomp = sim.PulseSequence([ \
        cnot13,
        sim.Hide(params, 3, True),
        sim.Hide(params, 4, True),
        Fred12,
        sim.Hide(params, 3, False),
        sim.Hide(params, 4, False),
        sim.Hide(params, 1, True),
        sim.Hide(params, 2, True),
        Fred34,
        sim.Hide(params, 1, False),
        sim.Hide(params, 2, False)])
    
    shor13b = sim.PulseSequence([ \
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
        sim.Hide(params, 4, False),
        cnot1234
        ])
    
    shor8b = sim.PulseSequence([ \
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
    
    shor7b = sim.PulseSequence([ \
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
        sim.Hide(params, 2, False), #vis: 1,2,3,4
        cnot1234 ])
    
    shor2b = sim.PulseSequence([ \
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
        sim.Hide(params, 2, False)]) #vis: 1,2,3,4
   
    shor7a = sim.PulseSequence([ \
        sim.Rcar(params, pi/2, 0),
        sim.Rac(params, -pi/2, 0, 0),
        sim.Rcar(params, pi/4, 0),
        sim.RMS(params, pi/8, pi/2),
        sim.Rac(params, -pi, 0, 2),
        sim.Rac(params, -pi, 0, 3),
        sim.Rcar(params, pi/4, 0),
        sim.RMS(params, -pi/8, pi/2),
        sim.Rac(params, -pi, 0, 0),
        sim.RMS(params, pi/8, pi/2),
        sim.Rac(params, -pi, 0, 2),
        sim.Rac(params, -pi, 0, 3),
        sim.RMS(params, -pi/8, pi/2),
        sim.Rcar(params, pi/2, pi),
        sim.Rac(params, -pi/2, 0, 0),
        sim.Rcar(params, pi/2, pi) ]) 

    shor7b1 = copy.deepcopy(Fredkin)
    shor7b2 = copy.deepcopy(Fredkin)
    shor7b1.changeions(params, (0,1,3))
    shor7b2.changeions(params, (0,2,4))
    shor7bold = sim.PulseSequence([ \
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
    def GeneratePulseSeq(a):
        ideal = None
    
        NOP = sim.Delay(params, 0.1)
        if a == 7:
            if(OLD):
                a2modN = shor7bold
                amodN = shor7a
            elif(HARD):
                a2modN = shor7bold
                amodN = shor7b
            else:
                a2modN = cnot24
                amodN = shor7b

            ideal = [0, 0,   0, .25, 0, 0, 0, 0, \
                     0, .25, 0, 0,   0, 0, 0, .25,
                     0, 0,   0, 0,   0, 0, 0, 0,
                     0, 0,   0, .25, 0, 0, 0, 0]
        elif a == 2:
            a2modN = cnot24
            amodN = shor2b
            ideal  = [0, 0,   0, 0,   0, 0,   0, 0, \
                      0, 0,   0, 0,   0, 0,   0, .25,
                      0, 0,   0, 0,   0, 0,   0, .25,
                      0, 0,   0, .25, 0, .25, 0, 0]
        elif a == 8:
            a2modN = cnot24
            amodN = shor8b
            ideal  = [0, 0,   0, 0,   0, 0,   0, 0, \
                      0, 0,   0, 0,   0, 0,   0, .25,
                      0, 0,   0, 0,   0, 0,   0, .25,
                      0, 0,   0, .25, 0, .25, 0, 0]
        elif a == 13:
            a2modN = cnot24
            amodN = shor13b
            ideal = [0, 0,   0, .25, 0, 0, 0, 0, \
                     0, .25, 0, 0,   0, 0, 0, .25,
                     0, 0,   0, 0,   0, 0, 0, 0,
                     0, 0,   0, .25, 0, 0, 0, 0]
        elif a == 11:
            a2modN = sim.Delay(params, 0.1)
            amodN = shor11
            ideal = [0, 0, 0, 0, 0, .5, 0, 0, \
                     0, 0, 0, 0, 0, 0,  0, .5,
                     0, 0, 0, 0, 0, 0,  0, 0,
                     0, 0, 0, 0, 0, 0,  0, 0]
    
        pulseseq_group = Kit.GeneratePulseSeq(params, NOP, a2modN, amodN)
        return pulseseq_group,ideal
    
    #######################################################
    
    pulseseq,ideal = GeneratePulseSeq(select_a)
    
    if doIdeal:
        for ps in pulseseq:
            for pulse in ps:
                pulse.use_ideal = True
    
    ### run it
    if doRun:
        tic = time.time()
        result = Kit.simulateevolution(pulseseq, params, dec, doPP=doPPKitaev)#, ideal=ideal
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
