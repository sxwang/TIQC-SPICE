''' scan of carrier transition varying detuning, showing spectrum '''

import numpy as np
import PyTIQC.core.simtools as sim
import PyTIQC.core.qctools as qc
import PyTIQC.tools.progressbar as progbar
import matplotlib.pyplot as pl

pi = np.pi

NumberOfIons = 1
NumberofPhonons = 1
hspace = sim.hspace(NumberOfIons,2,NumberofPhonons,0)
hspace.initial_state("thermal", nBar=0.5)
params = sim.parameters(hspace)
dec = sim.decoherence(params)

params.stepsize = 600

detuning = np.arange(-1.2*params.omz, 1.2*params.omz, params.omz/300)
#dets1 = np.arange(-1.1*params.omz, -0.9*params.omz, 2*pi*0.001)
#detc = np.arange(-2*pi*0.1, 2*pi*0.1, 2*pi*0.001)
#dets2 = np.arange(0.9*params.omz, 1.1*params.omz, 2*pi*0.001)
#detuning = np.hstack([ dets1, detc, dets2 ])
YR = np.zeros([len(detuning),hspace.dimensions], np.complex128)

widgets = [progbar.Percentage(), ' ', progbar.Bar(),' ', progbar.ETA()]
pbar = progbar.ProgressBar(widgets=widgets).start()

for i in range(len(detuning)):
    params.detuning = detuning[i]
    pulseseq = sim.PulseSequence( [ \
        sim.Rcar(params, pi/params.eta[0], 0), \
        ] )

    data = qc.simulateevolution(pulseseq, params, dec)
    YR[i,:] = data.YR[-1,:]

    pbar.update(int(1.*(i+1)*100/(len(detuning))))

data1 = sim.database(detuning, YR, hspace, pulseseq=None, statesvalid = False)
#data1.displaypopulation(1)
data1.tracedpopulation(0)
pl.plot(detuning/2/pi, 1-data1.YtrI)
pl.xlabel('relative frequency (2$\pi$ MHz)')
pl.ylabel('D state population')
pl.show()
