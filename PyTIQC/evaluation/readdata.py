#!/usr/bin/env python
# -*- mode: Python; coding: latin-1 -*-
# Time-stamp: "04-May-2010 14:31:55 c704215"

#  file       readdata.py
#  copyright  (c) Philipp Schindler 2008

"""readdata.py - Tools for file reading and path management
PathObject - Reads the Params.dat file and gets the information on the individual scans
ReadData - Reads a single datafile. It reads qc, csingle and cprb files if available
ReadDataMultiple - Reads multiple files at once
"""
import re
import os.path
import numpy
import copy
import numpy as np

class PathObject:
    """
    Opens a day of experiment data and loads the parameter file
    The path name is set globally if the aus() command is invoked
    
    PathObject.list_files - list files of the day with additional scan information
    """
    def __init__(self, path_name=None):
        assert type(path_name)==str, "Path name must be a string"
        self.parameters_file = "Params.txt"
        self.parameters = None
        self.current_file = None
        self.path_name = path_name
        self.readdata_instance = ReadData
        self.get_parameters()
        self.set_default()


    def list_files(self):
        "Prints the files of this day"
        self.reload()
        print(self.parameters)

    def get_range(self, start_time, stop_time=None):
        "Returns a list of filenames for a given range"
        qc_list = self.parameters.keys()
        qc_list.sort()
        time_list = []
        if type(start_time) == list:
            stop_time = start_time[1]
            start_time = start_time[0]
        start_time = int(start_time)
        stop_time = int(stop_time)
        for filename in  qc_list:
            time = int(filename[2:6])
            if time >= start_time and time <= stop_time:
                time_list.append(str(filename[2:-4]))
        return time_list

    def reload(self):
        "Reloads from parameter file"
        self.get_parameters()

    def get_parameters(self):
        "Sets the parameters dictionary"
        self.parameters = ReloadDict(self)
        file_name = self.path_name + self.parameters_file
        file_obj = open(file_name)
        file_string = file_obj.read()
        file_obj.close()
        params_array = file_string.split("Params for ")
        params_array.pop(0)
        for param_string in params_array:
            param_list = param_string.split("\n")
            scan_name = param_list.pop(0).strip()
            param_dict = {}
            info_string = param_list.pop(0)
            param_dict["scan information"] = info_string
            even_index = False
            sequence_name = param_list.pop(2)
            param_dict["sequence_name"] = sequence_name
            sequence_var_names = param_list[2].split("\t")
            param_dict["sequence_variables"] = sequence_var_names
            for param_line in param_list:
                param_line = param_line.strip()
                if even_index == False:
                    even_index = True
                    header_list = param_line.split("\t")
                else:
                    even_index = False
                    data_list = param_line.split("\t")

                    for index in range(len(data_list)):
                        try:
                            param_name = header_list[index]
                            if data_list[index] == "":
                                break
                            param_val = float(data_list[index])
                            param_dict[param_name] = param_val
                        except IndexError:
                            pass

            self.parameters[scan_name] = param_dict

    def get_next_file(self):
        "Outdated"
        name_array = self.parameters.keys()
        name_array.sort()
        if self.current_file == None:
            self.current_file = name_array[0]
            file_name = self.current_file
        else:
            try:
                index = name_array.index(self.current_file)
                file_name = name_array[index +1]
            except SyntaxError:
                file_name = self.current_file
        data_obj = self.readdata_instance(file_name)
        self.current_file = file_name
        return data_obj

    def get_last_file(self):
        "Returns the latest file of the day Outdated?"
        self.get_parameters()
        name_array = self.parameters.keys()
        name_array.sort()
        file_name = name_array.pop()
        data_obj = self.readdata_instance(file_name)
        return data_obj

    def get_files(self):
        "Returns a sorted list of the filenames of this day"
        try:
            a = self.parameters.keys()
            a.sort
            return a
        except:
            raise RuntimeError("Error: parameters file not defined")

    def set_default(self):
        "Set the default path for all ReadData objects"
        ReadData.default_path = self

class ReloadDict(dict):
    "A dictionary class which tries to reload if a key is not found"
    def __init__(self, parent=None):
        self.parent = parent

    def __getitem__(self, key):
        if len(key) < 10:
            key = "qc" + key + ".dat"
        try:
            ret_val = dict.__getitem__(self, key)
        except KeyError:
            self.parent.reload()
            ret_val = dict.__getitem__(self.parent.parameters, key)
        return ret_val

    def __str__(self):
        my_list = self.keys()
        my_list.sort()
        my_str = ""
        for item in my_list:
            this_str = item + " | "
            seq_info = copy.copy(self.__getitem__(item)['scan information'])
            seq_info = seq_info.split(":")[0].split(" ")[-1].strip('"') + " | "
            seq_name = copy.copy(self.__getitem__(item)['sequence_name'])
            seq_name = seq_name.split("\\")[-1]

            my_str += this_str + seq_info + seq_name + "\n"
        my_str.rstrip("\n")
        return my_str

class StringList(list):
    "Prints the elements of the list with the __str__method "
    def __str__(self):
        my_str = ""
        for item in self:
            my_str += str(item) + "\n"
        return my_str


def ReadDataMultiple(file_list, exclude_list=[], path=None,
                     is_range=False, sequence_name = None):
    """Loads multiple files at once
    file_list : List of the times to load
    exclude_list : List of the times to exclude
    is_range : file_list contains only begin and end time
    sequence_name : Add only files with this sequence file name
    path : Püath of the files

    Returns a list of Readdata objects
    """
    data_list = StringList()
    if type(file_list) != list:
        file_list = [file_list]
    if type(exclude_list) != list:
        exclude_list = [exclude_list]
    map(str, exclude_list)
    if path == None:
        path = ReadData.default_path
    if is_range:
        file_list = path.get_range(file_list)
    for time in file_list:
        if str(time) not in exclude_list:
            my_data = ReadData(time, path)
            seq_name = my_data.parameters["sequence_name"].split('\\')[-1]
            if seq_name[0:-1] == sequence_name or sequence_name == None:
                data_list.append(my_data)
    return data_list

class ReadData:
    """Base class for ReadData functions

    objects:

    xaxis_name : String for the x axis
    yaxis_name : String for the y axis
    data_dict  : Dictionary which contains the data
               Depending on the datafiles available following
               keys are possible:
                   qc      - Single ion PMT data
                   cprb    - Quantum state files
                   csingle - Single ion excitatiotns
    """

    default_path = None
    def __init__(self, file_name=None, path=None):
        "Do the whole magic"
        # Check the parameters
        self.header_dict = {}
#        assert type(file_name) == str, "Filename must be a string"
        # Init class variable defaults
#        self.type_string = None
        self.type_list = []
        self.file_array = None
        self.time_string = None
        self.xaxis_name = ""
        self.yaxis_name = ""
        self.data_dict = {}
        # The basic data types:
        self.data_type_list = ["qc", "csingle", "cprb"]
        if path == None:
            self.get_default_path()
        else:
            self.path = path
        if self.path == None:
            raise RuntimeError("No path object given")
        if file_name == None or type(file_name) == int:
            self.path.reload()
            file_list = self.path.parameters.keys()
            file_list.sort()
            if file_name == None:
                file_name = 1
            file_name = file_list[-file_name]
        
        self.file_name = file_name
        self.analyze_filename()
        self.get_data()
        self.plot_data()
        try:
            self.parameters = self.path.parameters['qc'+self.time_string+'.dat']
        except KeyError:
            raise RuntimeError("Error: file not found in params.txt: " +\
                                   'qc'+self.time_string+'.dat')

    def search_param(self, regex_string):
        """Searches the parameter dictionary for a certain string"""
        keys = self.parameters.keys()
        matched_keys = []
        p=re.compile(".*"+regex_string,re.I)
        for item in keys:
            if p.match(item) != None:
                matched_keys.append(item)
        for item in matched_keys:
            print item + ": " + str(self.parameters[item])


    def plot_data(self):
        "Dummy function"
        return None

    def get_default_path(self):
        "Some better path handling is needed"
        self.path = self.default_path

    def analyze_filename(self):
        "Check what type of file we've got"
        if self.file_name[-3:] != "dat" :
            self.time_string = str(self.file_name)
            self.file_name = "qc" + str(self.file_name) + ".dat"
            mytype_string = "qc"
        else:
            re_exp1 = re.compile("[0-9]")
            name_array = re_exp1.split(self.file_name)
            mytype_string = name_array[0]
            #now lets get the time argument. split by type_string and then by .dat
            re_exp2 = re.compile(mytype_string)
            re_exp3 = re.compile("\.dat")
            #take stuff after the type
            name_array1 = re_exp2.split(self.file_name)
            #and now take everything before the .dat
            name_array2 = re_exp3.split(name_array1[1])
            self.time_string = name_array2[0]
        for type_string in self.data_type_list:
            filename = self.path.path_name + str(type_string) + \
                       str(self.time_string) + ".dat"
            if os.path.isfile(filename):
                self.type_list.append(type_string)

    def load_file(self, type_string):
        "Load a file name and put it into a list"
        assert type(self.path.path_name)==str, "Path name must be a string: "\
                                + str(self.path.path_name)
        file_name = type_string + self.time_string + ".dat"
        full_name = self.path.path_name + file_name
        try:
            file_obj = open(full_name)
        except SyntaxError:
            raise RuntimeError("file not found: "+full_name)
        file_string = file_obj.read()
        file_obj.close()
        file_string = file_string.strip()
        #to fix the faulty header line of the cprb files
        if type_string == "cprb":
            file_string = file_string.replace("\r\n", "\n")
            file_string = file_string.replace("\r", "\n")
        file_list1 = file_string.split("\n")
        file_array = []
        for line in file_list1:
            file_array.append(line.split("\t"))

        header_list = file_array.pop(0)
        self.xaxis_name = header_list[0]
        self.yaxis_name = header_list[1]
        self.header_dict[type_string] = header_list
        array_lines = len(file_array)
        array_rows = len(file_array[0])
        my_data = numpy.zeros((array_lines, array_rows))
        for line in xrange(array_lines):
            for elem in xrange(array_rows):
                my_data[line, elem] = float(file_array[line][elem])
        self.data_dict[type_string] = my_data


    def get_data(self):
        """Reads all available datafiles"""
        if self.type_list == []:
            raise RuntimeError("File not found: " + 'qc'+ \
                                   self.time_string+'.dat')
        else:
            for type_string in self.type_list:
                self.load_file(type_string)

    def append(self, other):
        """Frontend for append single - also handles lists of data"""
        if type(other) != list:
            other = [other]
        for another in other:
            if another is not self:
                self.append_single(another)

    def append_single(self, other):
        """Append two datasets"""
        if self.parameters['sequence_name'] != other.parameters['sequence_name']:
            raise RuntimeError("Cannot append files of different sequences")
        for key in self.data_dict.keys():
            my_data = self.data_dict[key]
            other_data = other.data_dict[key]
            size1 = numpy.shape(my_data)
            size2 = numpy.shape(other_data)
            all_data = numpy.zeros((size1[0]+size2[0], size1[1]))
            all_data[0:size1[0], :] = my_data
            all_data[size1[0]:, :] = other_data
            self.data_dict[key] = all_data
        self.time_string += " app " + other.time_string

    def add_multiple(self, other_list):
        """Adds a list of data objects
        This is faster than the __add__ function"""
        otherme = copy.deepcopy(self)
        if len(other_list) == 0:
            other_list = [other_list]
        for other_data in other_list:
            if other_data != otherme: 
                otherme.__add__(other_data, use_copy=False)
            else:
                print "cannot add myself"
        return otherme

    def __add__(self, otherme, use_copy=True):
        #other = copy.deepcopy(otherme)
        other = otherme
        if use_copy:
            myself = copy.deepcopy(self)
        else:
            myself = self
        if (other == 0) or other == None:
            return myself
        if myself.parameters["scan information"] != \
               other.parameters["scan information"]:
            raise RuntimeError("Cannot add these two files")

        my_cycles = myself.parameters['cycles']
        his_cycles = copy.copy(other.parameters['cycles'])
        total_cycles = my_cycles + his_cycles
        for key in myself.data_dict.keys():
            my_data = myself.data_dict[key] * float(my_cycles)
            his_data = copy.copy(other.data_dict[key]) * float(his_cycles)
            myself.data_dict[key] = (my_data + his_data) / total_cycles
        myself.time_string += " + " + copy.copy(other.time_string)
        myself.parameters['cycles'] = total_cycles
        return myself

    def __str__(self):
        return "data object: qc" + self.time_string + ".dat"


def readcqst(filename, path =None, nr_of_splits=2, nr_of_ions=4):
    """20100810
    32 - SSSSS -1442
    31 - SSSSD -1449
    """
    if path == None:
        myqc = ReadData(filename, path=None)
        path = myqc.path.path_name
    filename = "cqst"+filename+'.dat'
    data = np.loadtxt(path+filename, skiprows=1)
    seq_par_list = data[:,0]
    nr_of_cycles = data.shape[1]
    print nr_of_cycles
    cprb_split_list = []
    for split_item in xrange(nr_of_splits):
        frame_cprb = np.zeros((len(seq_par_list), 2**nr_of_ions+1))
        frame_cprb[:,0] = data[:,0]
        cprb_split_list.append(frame_cprb)
        
    for seq_index in xrange(len(seq_par_list)):
        my_row = data[seq_index, 1:]
        mydat_array = [[],[]]
        for mydat in my_row:
            my_cqst_array =  seperate_data(mydat)
            for myindex in xrange(len(my_cqst_array)):
                mydat_array[myindex].append(my_cqst_array[myindex])

        for mysplit_index in xrange(len(mydat_array)):
            mydat = mydat_array[mysplit_index]
            my_count = 0
            for my_single in mydat:
                if my_single > 0:
                    my_count += 1
                    my_pos = convert_to_cprb_pos(my_single)
#                    print str(my_single) + " - " + str(my_pos)
                    cprb_split_list[mysplit_index][seq_index][my_pos+1] += 1
#                    cprb_split_list[0][seq_index][my_pos+1] += 1
#            print my_count
            if my_count == 0 :                   
                raise RuntimeError("Zero counts on one possibility detected - Need more data - duh")
    return cprb_split_list

def convert_to_cprb_pos(cqst_nr, nr_of_ions=4):
    """converts cqst number to cprb position"""
 #   cprb_pos = 2**(nr_of_ions) - cqst_nr 
    cprb_pos = cqst_nr - 1
    return cprb_pos

def seperate_data(mydat):        
    if mydat > 2**4:
        return [mydat - 2**4,-1]
    else:
        return [-1, mydat]

if __name__ == "__main__":
#    data = '/mnt/labeval/home/c704215/LinDaten/20100810/cqst1442.dat'
    data = '/mnt/labeval/home/c704215/LinDaten/20100804/cqst2325.dat'
    atarray = readcqst(data)
            
##
## readdata.py
## Login : <viellieb@ohm>
## Started on  Fri Aug  1 21:14:37 2008 Philipp Schindler
## $Id$
##
## Copyright (C) 2008 Philipp Schindler
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
##


# readdata.py ends here
