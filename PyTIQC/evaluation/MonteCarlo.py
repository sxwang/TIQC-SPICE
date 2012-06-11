"""MonteCarlo
A module for Monte Carlo analysis of State and Process tomographies

MonteCarloState: Monte Carlo analysis of a state tomography

Example:
>>> mc_obj = MonteCarloState('1610', rho_id)

MonteCarloProcess: Monte Carlo analysis of a process tomography:

Example:
>>> mc_obj = MonteCarloProcess('1610', chi_id)
"""

import types
import numpy as np
import readdata as rd
import densitymatrixreconstruction.IterML as IterML
import processtomography.proctom as proctom
import processtomography.quantumprocess as quantumprocess


class MonteCarloState:
    """Monte Carlo Analysis of a state tomography

    parameters: 
    
    raw_data - Either data object or time string
    ideal_object - density matrix to compare the result to
    
    nr_of_cycles=None - 
    """
    tomography_fun = IterML.iterfun_obj
    def __init__(self, raw_data, ideal_object, measure='fid',
                 nr_of_samples=100, path=None, nr_of_cycles=None, verbose=True):
        self.get_data(raw_data, path, nr_of_cycles)
        self.ideal_object= ideal_object
        self.measure = measure
        self.finished_samples = 0
        self.nr_of_samples = nr_of_samples
        self.verbose = verbose
        self.run_montecarlo()
        self.analyse_run()

    def add_samples(self, nr_of_samples=100):
        self.nr_of_samples = nr_of_samples
        self.run_montecarlo()
        self.analyse_run()

    def run_montecarlo(self, save_matrices=True):
        if save_matrices:
            self.matrix_list = []
        self.distance_list = []
        try:
            for index in xrange(self.nr_of_samples):
                if self.verbose:
                    print("Processing sample " + str(index) + " of " +str(self.nr_of_samples))
                data = self.generate_data()
                #print data
                matrix = self.tomography_fun(data,100)
                #print matrix.exp_chi
                distance = self.get_distance(matrix)
                #print distance
                self.distance_list.append(distance)
                if save_matrices:
                   self.matrix_list.append(matrix)
                self.finished_samples += 1
        except KeyboardInterrupt:
            print "Fetched Keyboard interupt. Exiting"

    def generate_data(self):
        multi_data = np.zeros(self.data_array.shape)
        nr_rows =self.data_array.shape[0]
        for row in xrange(nr_rows):
            this_data = self.data_array[row, 1:]
            multi_data[row, 1:] = np.random.multinomial(self.nr_of_cycles,this_data,1) / float(self.nr_of_cycles)
        multi_data[:,0] = self.data_array[:,0]
        return multi_data

    def get_distance(self, matrix):
        if self.measure == 'jozsafid':
            return abs(matrix.fid(self.ideal_object))
        elif self.measure == 'tracedist-rho':
            return abs(matrix.trdistance(self.ideal_object))
        elif self.measure == 'tracedist-pop':
            return abs(matrix.trdistancePop(self.ideal_object))
        elif self.measure == 'sso':
            return abs(matrix.sso(self.ideal_object))
        else:
            print "ERROR: measure", measure, "not defined"
            return 0

    def get_data(self, raw_data, path, nr_of_cycles=None):
        self.nr_of_cycles = nr_of_cycles
        if isinstance(raw_data, np.ndarray):
            self.data_array = raw_data
            self.nr_of_cycles = nr_of_cycles#
            if not nr_of_cycles:
                raise RuntimeError('Error: you mus specify the number of cycles')
            # Need to get nr of ions
            return
        if type(raw_data)== str:
            self.data_obj = rd.ReadData(raw_data, path=path)
        else:
            self.data_obj = raw_data
        self.data_array = self.data_obj.data_dict['cprb']

    def analyse_run(self):
        self.mean_distance = np.mean(self.distance_list)
        self.std_distance = np.std(self.distance_list)
        if self.verbose:
            print "Mean distance: " + str(self.mean_distance)
            print "Distance std deviation: " + str(self.std_distance)

class MonteCarloProcess(MonteCarloState):
    """Monte Carlo Analysis for a process tomography"""
    tomography_fun = quantumprocess.proctomo_obj
