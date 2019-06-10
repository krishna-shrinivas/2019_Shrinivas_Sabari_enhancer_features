#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Created on Tue Jul  3 19:19:40 2018
#@author: krishna-shrinivas
"""
Module of functions parse and pre-process the output data
structures from MD simulations. This step is required before
running analysis on data using process_MD_data
"""
import os
import gsd
import gsd.fl
import gsd.hoomd
import numpy

def convert_snap_shot_to_monomer_chain_map(snap,N):
    """
    Takes in a HOOMD snapshot and total number of monomers
    to construct a dictionary that maps each unique monomer
    to a parent chain (polymer). Further, each unique chain
    is mapped to it's type (A==DNA, B==TF, C ==coactivator)

    monomer_to_chain_map is a dictionary whose keys are are the
    monomer ID's (fixed during simulation), and corresponding values are
    unique chain IDs.

    chain_type maps each each chain ID to a chain type. The simulation
    chain types broadly stand for

    -   A == DNA chain
    -   B == TF molecule
    -   C == coactivator molecule

    """

    b_list = snap.bonds;
    p_list = snap.particles;
    monomer_to_chain_map = {};
    chain_type = {};
    N_bonds = b_list.N;

    if N>0:
        mol_count = 0;
        ptype = p_list.types[p_list.typeid[0]];
        monomer_to_chain_map[0] = (mol_count);
        if ptype is not 'B' and ptype is not 'C':
            chain_type[0] ='A' ;
        elif ptype is 'B':
            chain_type[0] = 'B';
        elif ptype is 'C':
            chain_type[0] = 'C';
    i=1;

    while ((i<=N_bonds-2)):
        flag =0;
        while((i<=N_bonds-1) and (flag==0)):
            if sum(b_list.group[i,:]-b_list.group[i-1,:]-1) ==0:
                monomer_to_chain_map[b_list.group[i,0]] =(mol_count);
                monomer_to_chain_map[b_list.group[i,1]] =(mol_count);

                i =i +1;

            else:
                flag =1;
                mol_count = mol_count +1;
                ptype = p_list.types[p_list.typeid[b_list.group[i,0]]];
                monomer_to_chain_map[b_list.group[i,0]] =(mol_count);
                monomer_to_chain_map[b_list.group[i,1]] =(mol_count);
                chain_type[mol_count] =ptype;
                i= i+1;
    for i in numpy.arange(N):
        if i not in monomer_to_chain_map.keys():
            mol_count=mol_count +1;
            monomer_to_chain_map[i] = mol_count;
            chain_type[mol_count] = p_list.types[p_list.typeid[int(i)]]

    return(monomer_to_chain_map,chain_type)

def return_data(folder):


    """
    Takes in the path to the parent output folder (folder)
    and automatically obtains the sorted list of
    parameters, paths to files, and input_parameters.

    **Input**

    -   folder = path to outermost listing of output folders

        Typically, the folder/* list various tested conditions
        and folder/condition/Data/* contains Observables, Figures, and Trajectory
        sub-directory with those specific types of files

    **Output**

    -   param_name =   The string name of the variable that is swept over. Default
        values is 'default' when only one condition exists.

    -   sorted_param_range = sorted list of the different values of parameters

    -   sorted_unique_paths =  list of paths to each Observable folder folder/condition/Data/Observables
        or folder/condition/Data/Observables/sub-condition/

    -   sorted_input_parameters =   list of input_params for each condition

    """

    sub_folder_format =['N_A','N_B','N_C','L','seed','seq_A','seq_B','seq_C'];
    a = sorted(os.listdir(folder));
    param_name = [];
    param_range = [];
    common_path = [];
    if len(a) > 1:
        first_folder = a[0].split('_')
        second_folder = a[1].split('_')
        for count in numpy.arange(len(first_folder)):
            if first_folder[count] != second_folder[count]:
                p =count;

        param_name = sub_folder_format[p];
        param_range = [];
        for sub_folder in a:
            param_range.append(sub_folder.split('_')[p])
            common_path.append(folder + sub_folder + '/Data/Observables/')
    else:
        sub_folder = a[0]
        common_path.append(folder + sub_folder + '/Data/Observables/')

    unique_paths =[];
    for path in common_path:
        a = os.listdir(path);
        if len(a) > 1:
            first_folder = a[0].split('_')
            second_folder = a[1].split('_')
            p =(list(set(first_folder) - set(second_folder))[0])
            index =  first_folder.index(p);
            param_name = ''.join(first_folder[0:index]);
            param_range = [];
            for sub_folder in a:
                param_range.append(sub_folder.split('_')[index])
                unique_paths.append(path + sub_folder + '/')
        else:
            unique_paths.append(path+a[0] +'/')


    """   If only single param_range , set default values """
    if (len(unique_paths)==1):
        param_name = 'default';
        param_range.append(1.0);

    input_parameters = [];
    for path in unique_paths:
        a = os.listdir(path);
        count =0;
        while (count<len(a)):
            if a[count].find('input') > -1:
                filename = path + a[count]

                count = len(a) + 1;
            count = count+1;
        input_param_list = {};
        with open(filename, 'r') as f:
            prefix ='';
            data = f.readlines();
        for line in data:
            line=line.strip();
            var_name,var_value = line.split(',');
            if var_name is '%':
                prefix = str(var_value);
            elif var_name.find('seq') > -1:
                input_param_list[prefix+var_name] = str(var_value);
            else:
                input_param_list[prefix+var_name] = float(var_value);
        input_parameters.append(input_param_list);

    if param_name.find('seq') == -1:
        param_range = [float(x) for x in param_range];

    count = 0;
    for path in unique_paths:
        traj_path = path.replace('Observables','Trajectory')
        p =os.listdir(traj_path);
        if p[0] == 'Movies':
            del(p[0]);
        if p[-1] == 'Movies':
            del(p[-1]);
        gsd_file =traj_path + p[-1];
        f = gsd.fl.GSDFile(gsd_file, 'rb')
        t = gsd.hoomd.HOOMDTrajectory(f)
        init_snap = t[0];
        N_monomers = input_parameters[count]['N_A']*len(input_parameters[count]['seq_A']) + input_parameters[count]['N_B']*len(input_parameters[count]['seq_B']) + input_parameters[count]['N_C']*len(input_parameters[count]['seq_C'])
        input_parameters[count]['N_monomers'] = int(N_monomers);
        (monomer_to_chain_map,chain_type) = convert_snap_shot_to_monomer_chain_map(init_snap,N_monomers);
        input_parameters[count]['MC_map'] = monomer_to_chain_map;
        input_parameters[count]['CT_map'] = chain_type;
        count= count+1;


    p = [];
    if param_name.find('seq') > -1:
        a = [len(x) for x in param_range];
        p = numpy.argsort(a)
    else:
        p = numpy.argsort(param_range)

    sorted_param_range = []
    sorted_input_parameters = [];
    sorted_unique_paths = [];
    for count in numpy.arange(len(param_range)):
        sorted_param_range.append(param_range[p[count]]);
        sorted_input_parameters.append(input_parameters[p[count]]);
        sorted_unique_paths.append(unique_paths[p[count]]);

    return (param_name, sorted_param_range,sorted_unique_paths,sorted_input_parameters)
