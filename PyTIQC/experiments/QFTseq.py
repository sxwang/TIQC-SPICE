''' pulse sequences for (fully-quantum) QFT '''

import PyTIQC.core.simtools as sim

QFTseqP1 = sim.PulseSequence( [ sim.Rcar(params, pi/2, -pi/2) ] )
QFTseqP2 = sim.PulseSequence( [ sim.Rcar(params, pi/2, -pi/2), 
                                sim.Rac(params, pi/2, 0, 0), 
                                sim.Rcar(params, pi/2, 0) ] )
QFTseqP3 = sim.PulseSequence( [ \
    sim.Rcar(params, pi/4, pi/2),  #U.Uy(pi/4, N),
    sim.Rac(params, pi, 0, 1),    #U.Uz(pi, N, 1),
    sim.Rcar(params, pi/4, -pi/2), #  U.Uy(-pi/4,  N),
    sim.Rcar(params, pi/2, pi/2),   #U.Uy(pi/2,N),
    sim.Rac(params, pi, 0, 2),     #U.Uz(pi, N, 2),
    sim.Rcar(params, pi/2, -pi/2),    #U.Uy(-pi/2,N),
    sim.RMS(params, pi/4, 0),    #U.MS(pi/4,0,N),
    sim.Rac(params, pi, 0, 1),   #U.Uz(pi, N, 1),
    sim.RMS(params, pi/4, 0),    #U.MS(pi/4,0,N)
]   )
QFTseqP4 = sim.PulseSequence( [ sim.Rcar(params, pi/4, -pi/2), 
                                sim.Rac(params, pi, 0, 2), 
                                sim.Rcar(params, pi/4, pi/2) ] )


QFTseq = sim.PulseSequence( [ \
    sim.Rcar(params, pi/2, 0),       #1
    sim.Rac(params, pi, 0, 1),       #2
    sim.Rac(params, pi/2, 0, 2),     #3
    sim.RMS(params, pi/8, 0),        #4
    sim.Rac(params, pi, 0, 2),       #5
    sim.RMS(params, pi/16, 0),       #6
    sim.Rcar(params, pi/16, -pi/2),  #7
    sim.Rac(params, pi, 0, 1),       #8
    sim.RMS(params, 3*pi/16, 0),     #9
    sim.Rcar(params, 3*pi/16, pi/2), #10
    sim.Rac(params, 3*pi/2, 0, 1),   #11
    sim.Rcar(params, pi/4, pi/2),    #12
    sim.RMS(params, pi/8, 0),        #13
    sim.Rac(params, pi, 0, 2),       #14
    sim.RMS(params, pi/8, 0),        #15
    sim.Rac(params, pi/2, 0, 0),     #16
    sim.Rac(params, pi, 0, 1),       #17
    sim.Rcar(params, pi/2, 0)        #18
    ] )
