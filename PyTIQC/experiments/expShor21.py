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
    pulseseqfileShor = 'experiments/Shorseq_6qubit.py'
   
    firstMultMap = True # Do first multiplication with CNOT 
                        # mapping instead of Multiplier 
                        # (does nothing in easy cases)

    select_a = 4 # Factoring base: 2,4,5,8,10,11,13,16,17,19

    doRun = True
    doIdeal = True
    printpulse = False
    saveKitaev = False
    doPPKitaev = True

    print 'N = 21, a =',select_a
    
    NumberOfIons = 6
    NumberOfPhonons = 0 if doIdeal else 7 
    hspace = sim.hspace(NumberOfIons,2,NumberOfPhonons,0)
    hspace.initial_state("quantum", qstate='DSSSSS')
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
    dec.doRandNtimes = 1 #16
    dec.doPP = True
    dec.dict['all'] = False
#    if doPPKitaev: params.ppcpus = 2
    
    Kit = Kitaev.Kitaev()
    
    ##########################################
    # load the pulse sequences
    # change ion order and define permutations
    ##########################################
    execfile(pulseseqfileShor,locals(),globals())

    cnot12_6 = sim.PulseSequence([ cnot12 ])

    Fred35 = copy.deepcopy(Fredkin)
    Fred13 = copy.deepcopy(Fredkin)
    Fred35.changeions(params, (3,5,0))
    Fred13.changeions(params, (1,3,0))
   
    cnot15 = copy.deepcopy(cnot12_6)
    cnot23 = copy.deepcopy(cnot12_6)
    cnot24 = copy.deepcopy(cnot12_6)
    cnot25 = copy.deepcopy(cnot12_6)
    cnot35 = copy.deepcopy(cnot12_6)
    cnot15.changeions(params,(0,1,5))
    cnot23.changeions(params,(0,2,3))
    cnot24.changeions(params,(0,2,4))
    cnot25.changeions(params,(0,2,5))
    cnot35.changeions(params,(0,3,5))


    times4 = sim.PulseSequence([ \
        sim.Hide(params, 1, True),
        sim.Hide(params, 2, True),
        sim.Hide(params, 4, True), #vis: 3,5
        Fred35,
        sim.Hide(params, 1, False),
        sim.Hide(params, 5, True), #vis: 1,3
        Fred13,
        sim.Hide(params, 2, False),
        sim.Hide(params, 4, False),
        sim.Hide(params, 5, False)
        ])

    times16 = sim.PulseSequence([ \
        sim.Hide(params, 2, True),
        sim.Hide(params, 4, True),
        sim.Hide(params, 5, True), #vis: 1,3
        Fred35,
        sim.Hide(params, 5, False),
        sim.Hide(params, 1, True), #vis: 3,5
        Fred13,
        sim.Hide(params, 1, False),
        sim.Hide(params, 2, False),
        sim.Hide(params, 4, False)
        ])

    ### general pulse sequence
    def GeneratePulseSeq(a):
    
        NOP = sim.PulseSequence([sim.Delay(params, 0.1)])

        if a in (2,5,16,19):    # a^2=4 cases
            if firstMultMap:
                a16modN = cnot15
            else:
                a16modN = times16
            a8modN = times4
            a4modN = times16
            a2modN = times4
            if a == 16:
                a1modN = times16
            else:
                a1modN = cnot24 # only cheating method implimented...

        elif a in (4,10,11,17): # a^2=16 cases
            if firstMultMap:
                a16modN = cnot13
            else:
                a16modN = times4
            a8modN = times16
            a4modN = times4
            a2modN = times16
            if a == 4:
                a1modN = times4
            else:
                a1modN = cnot24 # only cheating method implimented...

        elif a in (8,13):
            a16modN = NOP
            a8modN  = NOP
            a4modN  = NOP
            a2modN  = NOP
            if   a == 8:
                a1modN = cnot25
            elif a == 13:
                a1modN = cnot23

        oplist = [a16modN, a8modN, a4modN, a2modN, a1modN]

        pulseseq_group = Kit.GeneratePulseSeq(params, oplist)
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
