''' run simulation of QFT, load ideal densmat from file
and display resulting fidelity vs pulse '''

import numpy as np
import PyTIQC.core.simtools as sim
import PyTIQC.core.qctools as qc
import PyTIQC.core.gates as U
import matplotlib.pyplot as pl

import PyTIQC.evaluation.EvaluateData as evd

import time
import shelve, copy, inspect

pi = np.pi

### run params ##################################
pulseseqfile = 'experiments/QFTseq.py'
initstate = 1  # 1=sss, 2=p2, 3=p3
plottitle = "QFT"
displ = 1

doRun = True
doIdeal = False
calcfinal = True
doPP = True
savedata = False
calc_spontdecay_during_measurement = False

#################################################
files_dict = {}
idealrho = shelve.open('experiments/qft-ideal-rho.shlv')
if initstate == 1:     # SSS
    idealdata = idealrho['sss']
elif initstate == 2:   # P2
    idealdata = idealrho['p2']
elif initstate == 3:   # P3
    idealdata = idealrho['p3']
else:
    idealdata=None


NumberOfIons = 3
NumberOfPhonons = 7
hspace = sim.hspace(NumberOfIons,2,NumberOfPhonons,0)
params = sim.parameters(hspace)

params.use_servers( ['local'] )
params.shortestMS = 16
params.calcPulseParams()
params.progbar = False
params.savedata = savedata
#params.printpulse = True
params.addswitchtime = True

dec = sim.decoherence(params)
dec.doRandNtimes = 8
dec.dict['spontdecay'] = True
dec.doPP = doPP
dec.doPPprintstats=False

execfile(pulseseqfile)
pulseseq = QFTseq
if initstate == 2:
    QFTseqP2.extend(pulseseq)
    pulseseq = QFTseqP2
elif initstate == 3:
    QFTseqP3.extend(pulseseq)
    pulseseq = QFTseqP3

# turn off noise for some pulses
if doIdeal:
    for pulse in pulseseq:
        pulse.use_ideal = True

if doRun:
    tic = time.time()
    data = qc.simulateevolution(pulseseq, params, dec)
    data.RhoPNAll = np.array(data.RhoPNAll)
    toc = time.time()
    print "runtime: ", toc-tic, "seconds"
  
if params.savedata: 
    print "saving run data to file"
    qc.saveRun(pulseseq, params, dec, data, params.savedataname)

### evaluate and plot fidelities
if calcfinal:
    evalobj = evd.EvaluateData()
    if doRun:
        evalobj.loadidealdata(idealdata)
        evalobj.loadsimdata(data)
        evalobj.calculatePlotFidelities(displ,plottitle)#, grouperrors=100)
    else:
        print "no data included, doRun must be set to True"

### include spontaneous decay in measurement
if calc_spontdecay_during_measurement:
    pop_decayed = sim.DecayedPopulations_CCD(data.YRPN[-1], params)
    print pop_dcayed
