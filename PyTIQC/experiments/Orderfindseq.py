''' pulse sequences for order-finding '''

import PyTIQC.core.simtools as sim

# order 1 - Toffoli(1,3,2) - from Volkmarizer
order1 = sim.PulseSequence( [ \
sim.Rcar(params, pi/2, pi/2),  #0
sim.Rac(params, pi/4, 0, 1),  #1
sim.RMS(params, pi/2, pi/2),   #2
sim.Rcar(params, pi/2,pi), #3
sim.Rac(params, 3*pi/2, 0, 1), #4
sim.Rcar(params, pi/4,1*pi), #5
sim.RMS(params, pi/4, pi/2), #6
sim.Rac(params, pi/2, 0, 1), #7
sim.RMS(params, pi/2, pi/2), #8
sim.Rcar(params, pi/2,0*pi), #9
sim.Rcar(params, pi/2,-pi/2) ])

# order 2 - CNOT(1,3) - from paper (IV.A first one)
order2 = sim.PulseSequence( [ \
sim.Rcar(params, pi/2, 0),
sim.Rac(params, 7*pi/2, 0, 2), # 3pi/2 in experiment
sim.RMS(params, pi/4, pi/2),
sim.Rcar(params, pi/4, 0),
sim.Rac(params, 3*pi, 0, 1),
sim.Rcar(params, pi/4, 0),
sim.RMS(params, pi/4, pi/2),
sim.Rac(params, 7*pi/2, 0, 2),
sim.Rcar(params, pi/2, 0),
sim.Rac(params, 3*pi, 0, 1) ])

# order 3
order3 = sim.PulseSequence( [ \
sim.Rac(params, 7*pi/2, 0, 0),
sim.Rcar(params, 1.5*pi, 0.5*pi),
sim.Rac(params, 7*pi/2, 0, 0),
sim.RMS(params, pi/4, pi/2),
sim.Rcar(params, 1.25*pi, 0*pi),
sim.Rac(params, pi/2, 0, 2),
sim.Rcar(params, 0.5*pi, 1*pi),
sim.Rcar(params, 0.25*pi, 0.5*pi),
sim.Rac(params, 3*pi, 0, 1),
sim.Rcar(params, 0.25*pi, 0.5*pi),
sim.RMS(params, -pi/2, pi/2),
sim.Rac(params, 7*pi/2, 0, 1),
sim.Rac(params, 13*pi/4, 0, 0),
sim.Rac(params, pi, 0, 2),
sim.RMS(params, -pi/4, pi/2),
sim.Rcar(params, 0.19591327406*pi, -0.5*pi),
sim.Rac(params, 10*pi/3, 0, 1),
sim.Rcar(params, 0.19591327406*pi, 0.5*pi),
sim.Rcar(params, 0.25*pi, 1*pi),
sim.RMS(params, pi/2, pi/2),
sim.Rac(params, 0.25*pi, 0, 1),
sim.Rcar(params, 0.5*pi, 0.5*pi),
sim.Rac(params, 7*pi/2, 0, 2) ])

# order 3 ^2
order32 = sim.PulseSequence( [ \
sim.Rcar(params, 0.5*pi, 0.5*pi),
sim.Rac(params, 15*pi/4, 0, 1),
sim.Rcar(params, 0.5*pi, 0*pi),
sim.RMS(params, pi/2, pi/2),
sim.Rac(params, 0.5*pi, 0, 1),
sim.RMS(params, -pi/4, pi/2),
sim.Rcar(params, 0.25*pi, 1*pi),
sim.Rac(params, 15*pi/4, 0, 0),
sim.Rcar(params, 0.5*pi, 1*pi),
sim.Rac(params, 0.5*pi, 0, 1),
sim.RMS(params, pi/2, pi/2),
sim.Rac(params, pi/2, 0, 0),
sim.Rcar(params, 0.25*pi, 1*pi),
sim.RMS(params, pi/4, pi/2),
sim.Rac(params, pi/2, 0, 2),
sim.Rac(params, 0.5*pi, 0, 1),
sim.Rcar(params, 0.5*pi, 0*pi) ])

# order 4
order4 = sim.PulseSequence( [ \
sim.Rcar(params, 0.5*pi, 0*pi),
sim.Rcar(params, 1.0*pi, 0.5*pi),
sim.Rac(params, pi/2, 0, 2),
sim.RMS(params, -pi/8, pi/2),
sim.Rac(params, 3*pi, 0, 0),
sim.RMS(params, pi/8, pi/2),
sim.Rcar(params, 0.5*pi, 0.5*pi),
sim.Rcar(params, 1.5*pi, 0*pi),
sim.Rac(params, 7*pi/2, 0, 1),
sim.Rac(params, pi/2, 0, 2),
sim.RMS(params, -pi/4, pi/2), #10
sim.Rac(params, (4-2*13523./20456)*pi, 0, 0),
sim.Rcar(params, 0.19591327406*pi, -0.5*pi),
sim.Rac(params, 2*378./7999*pi, 0, 0),
sim.Rac(params, 0.666666666667*pi, 0, 2),
sim.Rcar(params, 0.19591327406*pi, 0.5*pi),
sim.Rcar(params, 0.25*pi, 0*pi),
sim.Rac(params, (4-2*63./15667)*pi, 0, 0), #17
sim.Rcar(params, 0.5*pi, -0.5*pi),
sim.RMS(params, -pi/2, pi/2),
sim.Rac(params, (4-1./3)*pi, 0, 1),  #20
sim.RMS(params, pi/2, pi/2),
sim.Rac(params, (2*1038./6827)*pi, 0, 0),
sim.Rcar(params, 0.5*pi, 0.5*pi),
sim.Rac(params, 0.25*pi, 0, 0) ])

# order 42 - CNOT(1,2) - from paper (IV.A first one)
order42 = sim.PulseSequence( [ \
sim.Rcar(params, pi/2, 0),
sim.Rac(params, pi/2, 0, 2),
sim.RMS(params, pi/4, pi/2),
sim.Rcar(params, pi/4, 0),
sim.Rac(params, pi, 0, 0),
sim.Rcar(params, pi/4, 0),
sim.RMS(params, pi/4, pi/2),
sim.Rac(params, pi/2, 0, 2),
sim.Rcar(params, pi/2, 0),
sim.Rac(params, pi, 0, 0) ])
