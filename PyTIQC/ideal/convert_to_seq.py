# this is the original Volckmarizer-to-PyTIQC translator written by Philipp.

"""convert_to_seq - Converts output from the Volckmarizer to a sequence in QCSL

It needs the number of ions and the smallestMS fraction that is set up.

Example usage:    nr_of_ions = 3
>>    smallest_ms = 1/4.
>>    hamil_str = ' [  1 2 4 3 6 3 2 1 5   4 3 6 2 6 4 2 1 6 2 3 5 3 6 2 6]'
>>    hamil_list = get_list_from_matlab_string(hamil_str)
>>    ampl_str = '[ 1/4 1/2 1/4 -1/32 -1/2 1/32 1/4 3/4 -1/4 1/4 -1/16 -13523/20456 -628/6411 378/7999 1/3 628/6411 1/8 -63/15667  -1/4 -1/8 -1/6 1/8 1038/6827 1/4 1/8]'
>>    ampl_list = get_list_from_matlab_string(ampl_str, operation=float)
>>    my_string = generate_sequence(hamil_list, ampl_list, nr_of_ions, smallest_ms)

Following conventions are used:

X(theta) - Rcar3(theta,[0,1 9,8)
Y(theta) - Rcar3(theta,[-.5,.5],8)
 %% theta_z = abs(theta) for theta<0 theta = 2-theta for theta>0
z(theta, ion) - Rzred(theta_z,0,ion)
MS(+theta) - RBichro1(theta/smallestMS*#2,#5,8,3)
MS(-theta) - RBichro2(theta/smallestMS*#3,#6,8,3)
"""

from __future__ import division
import copy

class QOperation:
    "Base class for quantum operations"
    phase_list = None
    is_ms = False
    ion = 8
    def __init__(self, theta, smallest_ms=1/4.):
        self.theta = theta
        self.smallest_ms = smallest_ms

    def get_string(self):
        if self.theta > 0:
            my_phase = self.phase_list[0]
        else:
            my_phase = self.phase_list[1]
        my_str = "Rcar3("+str(abs(self.theta))+","+my_phase+","+str(int(self.ion))+")"
        return my_str               

    def get_PyTIQC_string(self):
        if self.theta > 0:
            my_phase = self.phase_list[0]
        else:
            my_phase = self.phase_list[1]
        if int(self.ion) == 8:
            py_str = "sim.Rcar(params, "+str(abs(self.theta*2))+"*pi, "+my_phase+"*pi"+"),"
        else:
            py_str = "sim.Rcar(params, "+str(abs(self.theta*2))+"*pi, "+my_phase+"*pi,"+str(int(self.ion))+"),"
        return py_str 
        

class XOp(QOperation):
    phase_list = ['0','1']
    def __str__(self):
        return ('X ' + self.theta)

class YOp(QOperation):
    phase_list = ['0.5','-0.5']

    def get_string(self):
        my_str = "Rcar("+str(self.theta)+",-0.5,"+str(int(self.ion))+")"
        return my_str

    def __str__(self):
        return ('Y ' + self.theta)

class ZOp(QOperation):
    phase_list = ['0','0']
    def __init__(self, ion):
        self.ion = ion

    def __call__(self, theta):
        self.theta = theta
        return self

    def get_string(self):
        if self.theta < 0:
            theta = abs(self.theta)
        else:    
            theta = 2 - self.theta 
        my_str = "Rzred("+str(theta)+",-0.5,"+str(int(self.ion))+")"
        return my_str

    def get_PyTIQC_string(self):
        theta2 = self.theta*2
        if self.theta < 0:
            theta = 2 - abs(theta2)
        else:    
            theta = theta2
        py_str = "sim.Rac(params, "+str(theta)+"*pi, 0, "+str(2-(int(self.ion)-1))+"),"
        return py_str        

    def __str__(self):
        return ('z ' + self.theta + ' ' +self.ion)


class MSOp(QOperation):
    is_ms = True
    var_name_ms1_amp = '#2'
    var_name_ms1_phase = '#5'
    var_name_ms2_amp = '#3'
    var_name_ms2_phase = '#6'

    phase_list = ['pi/2', '0']

    def __init__(self, smallest_ms=1/4.):
        self.smallest_ms = smallest_ms
    
    def __call__(self, theta):
        self.theta = theta
        return self

    def get_string(self):
        ms_multiplication = abs(self.theta/self.smallest_ms)
        if self.theta > 0:
            my_str = "Rbichro1("+str(ms_multiplication)+"*" +str(self.var_name_ms1_amp) \
                +","+str(self.var_name_ms1_phase)+","+str(int(self.ion))+")"
        else:
            self.theta = abs(self.theta)
            my_str = "Rbichro2("+str(ms_multiplication)+"*" +str(self.var_name_ms2_amp) \
                +","+str(self.var_name_ms2_phase)+","+str(int(self.ion))+")"            
        return my_str               

    def get_PyTIQC_string(self):
        if self.theta > 0:
            my_phase = self.phase_list[0]
        else:
            my_phase = self.phase_list[1]
        py_str = "sim.RMS(params, "+str(self.theta*4)+"*pi, "+my_phase+"),"
        return py_str

    def __str__(self):
        return ('MS ' + self.theta)


def init_operator_dict(nr_of_ions, smallest_ms=1/4.):
    op_dict = {}
    op_dict[1] = XOp
    op_dict[2] = YOp
    op_dict[3] = MSOp(smallest_ms=smallest_ms)
    for index in xrange(nr_of_ions):
        op_dict[4+index] = ZOp(index+1)
    return op_dict

def generate_sequence(hamil_list, ampl_list, nr_of_ions, smallest_ms):
    op_dict = init_operator_dict(nr_of_ions, smallest_ms)
    pulse_list = []
    for index in xrange(len(hamil_list)):
        my_hamil = hamil_list[index]
        my_ampl = ampl_list[index]
        my_pulse = op_dict[my_hamil](my_ampl)
        pulse_list.append(my_pulse)
    sequence_string = '\n %%%%% Automatically generated sequence - recheck !!!! %%%%%%%%% \n\n'
    for item in pulse_list:
        sequence_string += item.get_string() + '\n'
    for item in pulse_list:
        sequence_string += item.get_PyTIQC_string() + '\n'
    return sequence_string

def get_list_from_matlab_string(mat_str, operation=int):
    mat_str = mat_str.lstrip().rstrip()
    str_list = mat_str.lstrip('[').rstrip(']').split()
    number_list =[]
    for item in str_list:
        number_list.append(operation(eval(item)))
    return number_list

    
def Fredkin():
    nr_of_ions = 3
    smallest_ms = 1/16
    hamil_str = '[ 2 6 3 5 6 1 3 5 3 1 5 6 3 5 4 1 4 5]'
    ampl_str = '[ 1/4 1/4 1/8 1/4 -1/4 -3/8 -1/16 1/4 1/8 1/4 -1/8 1/4 1/8 1/4 1/4 1/4 -1/4 -1/4]'
    hamil_list = get_list_from_matlab_string(hamil_str)
    ampl_list = get_list_from_matlab_string(ampl_str, operation=float)
    my_string = generate_sequence(hamil_list, ampl_list, nr_of_ions, smallest_ms)
    print "\n%%%%%%%% Fredkin %%%%%%%%%%%%%%%"
    print my_string

def Toffoli():
    nr_of_ions = 3
    smallest_ms = 1/16
    hamil_str = '[ 2 5 3 1 5 1 3 5 3 1 2 ]'
    ampl_str = '[1/4 1/8 1/8 -1/4 -1/4 -1/8 1/16 1/4 1/8 1/4 -1/4]'
    hamil_list = get_list_from_matlab_string(hamil_str)
    ampl_list = get_list_from_matlab_string(ampl_str, operation=float)
    my_string = generate_sequence(hamil_list, ampl_list, nr_of_ions, smallest_ms)
    print "\n%%%%%%%% Toffoli %%%%%%%%%%%%%%%"
    print my_string

def Cnot13():
    nr_of_ions = 3
    smallest_ms = 1/32
    hamil_str = '[ 1 4 3 6 3 1 4 3 6 3 1 4 1]'
    ampl_str = '[1/4 -1/4 -1/32 1/2 1/32 1/8 1/2 -1/32 1/2 1/32 -1/8 -1/4 -1/4]'
    hamil_list = get_list_from_matlab_string(hamil_str)
    ampl_list = get_list_from_matlab_string(ampl_str, operation=float)
    my_string = generate_sequence(hamil_list, ampl_list, nr_of_ions, smallest_ms)
    print "\n%%%%%%%% CNOT(1,3) %%%%%%%%%%%%%%%"
    print my_string


def Order3():
    nr_of_ions = 3
    smallest_ms = 1/16
    hamil_str = '[ 6 2 6 3 1 4 1 2 5 2 3 5 6 4 3 2 5 2 1 3 5 2 4] '
    ampl_str = '[ -1/4 3/4 -1/4 1/16 5/8 1/4 -1/4 1/8 -1/2 1/8 -1/8 -1/4 -3/8 1/2 -1/16 -628/6411 -1/3 628/6411 -1/8 1/8 1/8 1/4 -1/4]'
    hamil_list = get_list_from_matlab_string(hamil_str)
    ampl_list = get_list_from_matlab_string(ampl_str, operation=float)
    my_string = generate_sequence(hamil_list, ampl_list, nr_of_ions, smallest_ms)
    print "\n%%%%%%%% Order3 %%%%%%%%%%%%%%%"
    print my_string


def Order4():
    nr_of_ions = 3
    smallest_ms = 1/16
    hamil_str = '[1 2 4 3 6 3 2 1 5 4 3 6 2 6 4 2 1 6 2 3 5 3 6 2 6]'
    ampl_str = '[ 1/4 1/2 1/4 -1/32 -1/2 1/32 1/4 3/4 -1/4 1/4 -1/16 -13523/20456 -628/6411 378/7999 1/3 628/6411 1/8 -63/15667 -1/4 -1/8 -1/6 1/8 1038/6827 1/4 1/8 ]'
    hamil_list = get_list_from_matlab_string(hamil_str)
    ampl_list = get_list_from_matlab_string(ampl_str, operation=float)
    my_string = generate_sequence(hamil_list, ampl_list, nr_of_ions, smallest_ms)
    print "\n%%%%%%%% Order4 %%%%%%%%%%%%%%%"
    print my_string

def Order32():
    nr_of_ions = 3
    smallest_ms = 1/16
    hamil_str = '[ 2 5 1 3 5 3 1 6 1 5 3 6 1 3 4 5 1]'
    ampl_str = '[ 1/4 -1/8 1/4 1/8 1/4 -1/16 -1/8 -1/8 -1/4 1/4 1/8 1/4 -1/8 1/16 1/4 1/4 1/4] '
    hamil_list = get_list_from_matlab_string(hamil_str)
    ampl_list = get_list_from_matlab_string(ampl_str, operation=float)
    my_string = generate_sequence(hamil_list, ampl_list, nr_of_ions, smallest_ms)
    print "\n%%%%%%%% Order3^2 %%%%%%%%%%%%%%%"
    print my_string

def Shora7a():
    nr_of_ions = 5
    smallest_ms= 1/8
    hamil_str = '[2 1 4 7 1 6 5 3 5 3 2 6 7 3 6 1 3 4 3 5 3 5 7 3 1 6 3 5 6 2 7 2 6 5]'
    ampl_str = '[ -1/4 1/8 1/2 1053/7385 1/8 1/4 -1/4 1/32 1/2 1/32 -1/2 -1/2 1/2 1/32 1/2 -1/8 1/32 1/2 1/32 1/2 1/32 -1/2 1/2 1/32 -1/8 -1/2 1/32 1/4 -1/4 -1123/7245 -22108/30707 -1517/15969 1/4 1/4]'
    hamil_list = get_list_from_matlab_string(hamil_str)
    ampl_list = get_list_from_matlab_string(ampl_str, operation=float)
    my_string = generate_sequence(hamil_list, ampl_list, nr_of_ions, smallest_ms)
    print "\n%%%%%%%% Shor a=7(a) %%%%%%%%%%%%%%%"
    print my_string

def Shora11():
    nr_of_ions = 5
    smallest_ms = 1/8
    hamil_str = '[1 8 3 7 3 1 8 3 7 3 1 8 1 1 8 3 5 3 1 8 3 5 3 1 8 1]'
    ampl_str = '[1/4     -1/4    -1/32      1/2     1/32      1/8      1/2    -1/32      1/2     1/32     -1/8     -1/4     -1/4    1/4     -1/4    -1/32      1/2     1/32      1/8      1/2    -1/32      1/2     1/32     -1/8     -1/4     -1/4]'
    hamil_list = get_list_from_matlab_string(hamil_str)
    ampl_list = get_list_from_matlab_string(ampl_str, operation=float)
    my_string = generate_sequence(hamil_list, ampl_list, nr_of_ions, smallest_ms)
    print "\n%%%%%%%% Shor a=11 %%%%%%%%%%%%%%%"
    print my_string
 
if __name__ == '__main__':
#    Fredkin()
    Toffoli()
#    Cnot13()
#    Order3()
#    Order4()
#    Order32()
#    Shora7a()
#    Shora11()
