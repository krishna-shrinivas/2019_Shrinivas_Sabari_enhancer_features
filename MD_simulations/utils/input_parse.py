#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 13 19:05:28 2018

@author: krishna
"""

import os

def read_particle_type_information(pt_file):
    
    """
        read_particle_type_information takes in the path of the particle
        information file (pt_file)
        
        A dictionary of particle information is read line by line for each monomer,
        using the following syntax for any line:
            
            type,mass,diameter,charge
            
        TODO: An important check needs to be introduced for default definitions of
        monomer types
    """
    
    
    pt_info ={};
    pt_info['type'] = [];
    pt_info['mass'] = {};
    pt_info['diameter'] = {};
    pt_info['charge'] = {};
    pt_info['is_active_DNA'] = {};

    with open(pt_file, 'r') as f:
        for line in f:
            if line[0] is not '%':
                line=line.strip();
                var = line.split(',');
                pt_info['type'].append(var[0]);
                pt_info['mass'][str(var[0])]=(float(var[1]));
                pt_info['diameter'][str(var[0])]=(float(var[2]));
                pt_info['charge'][str(var[0])]=(float(var[3]));               
                pt_info['is_active_DNA'][str(var[0])]=(int(var[4]));               

                
    return(pt_info)
    
    
def write_input_params(file_output,input_params):
    """
        write_input_params writes all the input information from
        input_params into file_output in the following syntax
            key,value\n
    """    

    with open(file_output,'a') as f:
        for key in input_params.keys():
            f.write( ''.join(key)+','+str(input_params[key])+'\n');
    f.close()

def input_parse(filename):

    

    input_parameters  ={};
    with open(filename, 'r') as f:
        for line in f:
            line=line.strip();
            var_name,var_value = line.split(',');
            input_parameters[var_name] = (var_value);
#            print("Value of {} is {}".format(var_name,input_parameters[var_name]))
    return input_parameters;

def params_parse(filename):
    """
        Reads input parameters from a file of form    param_rep,param_value\n
        into a list and returns the name as well as values.
    """

    param_list  =[];
    with open(filename, 'r') as f:
        param_name = f.readline().strip();
        for line in f:
            if param_name.find('seq') > -1:
                param_list.append(str(line.strip()))
            else:
                param_list.append(float(line.strip()))
    return (param_name,param_list);

def energetics_parse(filename):
    """
        Reads pairwise LJ interactions from the file passed
        and returns a dictionary with the value of each interaction
        stored under a key tuple of form (typei,typej)
    """
    interaction_list = {};
    with open(filename,'r') as f:
        for line in f:
            if line[0] is not '%':
                line = line.strip();
                var  =  line.split(',');
                interaction_list[var[0],var[1]] = float(var[2]);
                interaction_list[var[1],var[0]] = float(var[2]);
                
    return interaction_list;