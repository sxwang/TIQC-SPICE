''' pulse sequences for Shor '''

import PyTIQC.core.simtools as sim

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
