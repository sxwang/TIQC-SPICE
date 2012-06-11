import numpy as np
import shelve

pi = np.pi

root = '/home/sxwang/TIQC-SPICE/MathCode/dephasefid'

dephase_files_dict = {
'3'+str(('M', pi/8)): root+'3q-Mpi8.csv',
#'3'+str(('M', pi/16)): root+'3q-Mpi16.csv',
'3'+str(('M', 3*pi/16)): root+'3q-Mpi316.csv',
'3'+str(('R', pi/2)): root+'3q-pi2.csv',
'3'+str(('R', pi/4)): root+'3q-pi4.csv',
'3'+str(('R', pi/16)): root+'3q-pi16.csv',
'2'+str(('M', pi/2)): root+'2q-Mpi2.csv',
'2'+str(('R', pi/2)): root+'2q-pi2.csv',
'2'+str(('R', pi/4)): root+'2q-pi4.csv',
'2'+str(('R', pi)): root+'2q-pi.csv',
'1'+str(('R', pi)): root+'1q-pi.csv',
'1'+str(('R', pi/2)): root+'1q-pi2.csv',
'1'+str(('R', pi/4)): root+'1q-pi4.csv'
}

dephase_dict = {}

for key, fname in dephase_files_dict.items():
    dephase_dict[key] = np.loadtxt(fname, delimiter=',')

d = shelve.open("evaluation/dephasefidconv.shlv")
for key, csv in dephase_dict.items():
    print key, np.shape(csv)
    d[key] = csv
d.close()
