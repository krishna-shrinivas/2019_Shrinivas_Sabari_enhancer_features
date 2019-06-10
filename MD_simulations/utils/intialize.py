#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 13 19:03:14 2018

@author: krishna
"""
import numpy
import itertools

def random_initialize(position,L,N_min,N_end,min_gap):
    mol =N_min;
    while(mol<N_end):
        p = min_gap*numpy.floor(numpy.random.randint(low=  numpy.ceil(-L/2) ,high= numpy.floor(L/2), size=(1,3))/min_gap);
        if not any(numpy.equal(position,p).all(1)):
            position[mol,:]=p ;
            mol = mol+1;

    return position;

def generate_cubic_position(r):
    """
    This function generates the cubic positions along (+-r,0,0), (+-r,+-r,0),
    and (+-r,+-r,+-r) as a numpy of zeros (26,3) returns a;
    
    """


    a = numpy.zeros((26,3));
    a[0:3,:]=-1*((numpy.array(list(set(itertools.permutations([0,0,r],3))))));
    a[3:6,:]= numpy.array((list(set(itertools.permutations([0,0,r],3)))));
    a[6:12,:]=numpy.array((list(set(itertools.permutations([0,r,-r],3)))));
    a[12:15,:]=numpy.array((list(set(itertools.permutations([0,-r,-r],3)))));
    a[15:18,:] = numpy.array((list(set(itertools.permutations([0,r,r],3)))));
    a[18:21,:] = numpy.array((list(set(itertools.permutations([-r,r,r],3)))));
    a[21:24,:] = numpy.array((list(set(itertools.permutations([-r,-r,r],3)))));
    a[24,:] = numpy.array((list(set(itertools.permutations([-r,-r,-r],3)))));
    a[25,:] = numpy.array((list(set(itertools.permutations([r,r,r],3)))));
    return a;

def generate_cubic_face_position(l,d):

    """
    This function returns the cubic patches along the poles of a sphere 
    of radius l  of patch size d     
    """    
    
    a = numpy.zeros((6,3));
    a[[0,2,4],:]=-1* numpy.array((list(set(itertools.permutations([0,0,0.5*(l-d)],3)))));
    a[[1,3,5],:]= numpy.array((list(set(itertools.permutations([0,0,0.5*(l-d)],3)))));
    return a;


def initialize_system(N_chains,seq_chains,pt_info,position,particle_info,pos_params,min_gap_multi=1):

    
    """
    This function initializes the box with N_chains of sequence seq_chains with
    monomer information in pt_info whose positions are stored in the position
    matrix, particle_info, and pos_params gives the size of the box passed. 
    """
    
    particle_types = pt_info['type'];
    diameter_list = pt_info['diameter'];
    mass_list = pt_info['mass'];
    charge_list = pt_info['charge'];
    is_active_DNA_list = pt_info['is_active_DNA']
    pos = generate_cubic_position(diameter_list['A']/3);

    mol_start = particle_info['Nchain_min'];
    bond_start = particle_info['bond_start'];
    monomer = int(particle_info['N_mono_min']);
    L = pos_params['L'];

    mol = mol_start;
    
    # Iterate over all the chains
    while (mol < mol_start + N_chains):
        
        # Starts at the beginning of each chain and obtains the relevant
        #   details of the first monomer of the chain
        chain_pos =0;
        mon_type=  particle_types.index(seq_chains[chain_pos]);
        particle_info['typeid'][monomer] = mon_type;
        particle_info['diameter'][monomer] = diameter_list[seq_chains[chain_pos]];
        particle_info['mass'][monomer] = mass_list[seq_chains[chain_pos]];
        particle_info['charge'][monomer] = charge_list[seq_chains[chain_pos]];
        particle_info['is_active_DNA'][monomer] = is_active_DNA_list[seq_chains[chain_pos]];

        flag = 0;
        while(flag==0):

            p = L*numpy.random.random_sample((1,3)) - L/2;
            mon = 0;

            while((mon<monomer) and (flag==0)):
                                            
                if seq_chains[chain_pos] == 'A':                    
                    if particle_types[int(particle_info['typeid'][mon])] == 'A':

                        scaffold = 0;
                        scaffold_mon =0
                        while ((scaffold < pos.shape[0]) and (flag==0)):
                            while((scaffold_mon < pos.shape[0]) and (flag==0)):
                               
                                pos_monomer = pos[scaffold,:] + p;
                                pos_mon = pos[scaffold_mon,:] + position[mon,:] ;

                                min_gap = 0.5* (particle_info['diameter'][monomer]/3+particle_info['diameter'][mon]/3);
                                dist = numpy.power(numpy.sum(numpy.power(pos_monomer-pos_mon,2)),0.5);
                                if dist< min_gap:
                                    flag =1;
                                
                                scaffold_mon = scaffold_mon+1;
                            
                            scaffold = scaffold+1
                    else:
                        scaffold = 0;
                        while ((scaffold < pos.shape[0]) and (flag==0)):
                            pos_monomer = pos[scaffold,:] + p;
                            min_gap = 0.5* (particle_info['diameter'][monomer]/3+particle_info['diameter'][mon]);
                            dist = numpy.power(numpy.sum(numpy.power(position[mon,:] - pos_monomer,2)),0.5);
                            if dist< min_gap:
                                flag =1;
                            scaffold =scaffold+1;
                else:
                    if particle_types[int(particle_info['typeid'][mon])] == 'A':
                        scaffold_mon =0;
                        while((scaffold_mon < pos.shape[0]) and (flag==0)):
                            pos_mon = pos[scaffold_mon,:] + position[mon,:];
                            min_gap = 0.5* (particle_info['diameter'][monomer]+particle_info['diameter'][mon]/3);
                            dist = numpy.power(numpy.sum(numpy.power(pos_mon - p,2)),0.5);
                            if dist< min_gap:
                                flag =1;
                            
                            scaffold_mon = scaffold_mon+1;

                dist = numpy.power(numpy.sum(numpy.power(position[mon,:] - p,2)),0.5);
                min_gap = 0.5* (particle_info['diameter'][monomer]+particle_info['diameter'][mon]);
                if dist< min_gap:
                    flag =1;    


                mon = mon+1;
            
            flag = flag-1;

            if flag<0:
                position[monomer,:] = p;                    
                monomer = monomer+1;
                chain_pos = chain_pos+1;

        
        while(chain_pos < len(seq_chains)):
        
            mon_type= particle_types.index(seq_chains[chain_pos]);
            particle_info['typeid'][monomer] = mon_type;
            particle_info['diameter'][monomer] = diameter_list[seq_chains[chain_pos]];
            particle_info['mass'][monomer] = mass_list[seq_chains[chain_pos]];
            particle_info['charge'][monomer] = charge_list[seq_chains[chain_pos]];
            particle_info['is_active_DNA'][monomer] = is_active_DNA_list[seq_chains[chain_pos]];

            flag = 0;
            while(flag==0):
                bond_gap = 1.2 * (particle_info['diameter'][monomer]+particle_info['diameter'][monomer-1]);
                p = position[monomer-1,:]+bond_gap*numpy.random.random_sample((1,3)) - 0.5*bond_gap;
                
                if ((p.min()< -L/2) or (p.max() > L/2)):
#                     print("This is not valid {}".format(p));
                    flag=1;
                    
                mon = monomer-1;
                while((mon>=0) and (flag==0)):
                                                
                    if seq_chains[chain_pos] == 'A':                    
                        if particle_types[int(particle_info['typeid'][mon])] == 'A':
 
                            # C.O.M distance between A centers
#                            min_gap = 0.5* (particle_info['diameter'][monomer]+particle_info['diameter'][mon]);
#                            dist = numpy.power(numpy.sum(numpy.power(p-position[mon,:],2)),0.5);
#                            if dist< min_gap:
#                                flag =1;
                           
                           #Now repeat this calculation for the scaffolding particles 
                            
                            
                            scaffold = 0;
                            scaffold_mon =0
                            while ((scaffold < pos.shape[0]) and (flag==0)):
                                while((scaffold_mon < pos.shape[0]) and (flag==0)):
                                   
                                    pos_monomer = pos[scaffold,:] + p;
                                    pos_mon = pos[scaffold_mon,:] + position[mon,:] ;
    
                                    min_gap = 0.5* (particle_info['diameter'][monomer]/3+particle_info['diameter'][mon]/3);
                                    dist = numpy.power(numpy.sum(numpy.power(pos_monomer-pos_mon,2)),0.5);
                                    if dist< min_gap:
                                        flag =1;
                                    
                                    scaffold_mon = scaffold_mon+1;
                                
                                scaffold = scaffold+1
                        else:
                            scaffold = 0;
                            while ((scaffold < pos.shape[0]) and (flag==0)):
                                pos_monomer = pos[scaffold,:] + p;
                                min_gap = 0.5* (particle_info['diameter'][monomer]/3+particle_info['diameter'][mon]);
                                dist = numpy.power(numpy.sum(numpy.power(position[mon,:] - pos_monomer,2)),0.5);
                                if dist< min_gap:
                                    flag =1;
                                scaffold =scaffold+1;
                    else:
                        if particle_types[int(particle_info['typeid'][mon])] == 'A':
                            scaffold_mon =0;
                            while((scaffold_mon < pos.shape[0]) and (flag==0)):
                                pos_mon = pos[scaffold_mon,:] + position[mon,:];
                                min_gap = 0.5* (particle_info['diameter'][monomer]+particle_info['diameter'][mon]/3);
                                dist = numpy.power(numpy.sum(numpy.power(pos_mon - p,2)),0.5);
                                if dist< min_gap:
                                    flag =1;
                                
                                scaffold_mon = scaffold_mon+1;
    
                
                    dist = numpy.power(numpy.sum(numpy.power(position[mon,:] - p,2)),0.5);
                    min_gap = 0.5* (particle_info['diameter'][monomer]+particle_info['diameter'][mon]);
                    if dist< min_gap:
                        flag =1;    
    
    
                    mon = mon-1;
                
                flag = flag-1;

                
                if flag<0:
#                     print("Dist between {} and {} is {}".format(monomer-1,monomer,numpy.power(numpy.sum(numpy.power(position[monomer-1,:] - p,2)),0.5)))
                  
                    position[monomer,:] = p;
                    particle_info['bonds'][bond_start] = numpy.array([monomer-1,monomer]);
                    bond_type = seq_chains[chain_pos]+'_'+seq_chains[chain_pos-1];
                    if bond_type not in particle_info['bond_types']:
                        particle_info['bond_types'].append(bond_type) ;
                    
                    bond_id = particle_info['bond_types'].index(bond_type);
                    particle_info['bond_id'][bond_start] = bond_id;
                    
                    bond_start = bond_start +1;
                    monomer = monomer+1;
                    chain_pos = chain_pos+1;

        
        mol = mol+1;
        
        
    
    
    particle_info['Nchain_min'] = mol;
    particle_info['N_mono_min'] = monomer;  
    particle_info['bond_start'] =bond_start;
    
    return (particle_info,position)




def initialize_system_complete(N_chains,seq_chains,pt_info,position,particle_info,pos_params,min_gap_multi=1):

    
    """
    This function initializes the box with N_chains of sequence seq_chains with
    monomer information in pt_info whose positions are stored in the position
    matrix, particle_info, and pos_params gives the size of the box passed. 
    """
    
    """
        Store all the information on particle type, diameter, mass, charge, and whether
        the monomer is an active piece of DNA or not.

        mol_start stores the starting point of chains to intialize, so do bond_start
        and monomer respectively.
    
    """    
    particle_types = pt_info['type'];
    diameter_list = pt_info['diameter'];
    mass_list = pt_info['mass'];
    charge_list = pt_info['charge'];
    is_active_DNA_list = pt_info['is_active_DNA']
    pos = generate_cubic_position(diameter_list['A']/3);

    mol_start = particle_info['Nchain_min'];
    bond_start = particle_info['bond_start'];
    monomer = int(particle_info['N_mono_min']);
    L = pos_params['L'];

    mol = mol_start;
    
    # Iterate over all the chains
    while (mol < mol_start + N_chains):
        
        #   Starts at the beginning of each chain and obtains the relevant
        #   details of the first monomer of the chain
        chain_pos =0;
        mon_type=  particle_types.index(seq_chains[chain_pos]);
        particle_info['typeid'][monomer] = mon_type;
        particle_info['diameter'][monomer] = diameter_list[seq_chains[chain_pos]];
        particle_info['mass'][monomer] = mass_list[seq_chains[chain_pos]];
        particle_info['charge'][monomer] = charge_list[seq_chains[chain_pos]];
        particle_info['is_active_DNA'][monomer] = is_active_DNA_list[seq_chains[chain_pos]];

        flag = 0;
        while(flag==0):

            p = L*numpy.random.random_sample((1,3)) - L/2;
            mon = 0;

            """
                Loops over all the previous monomers effectively to ensure that placement
                of current monomer does not have any excluded volume with already intialized
                monomers.
            """
            while((mon<monomer) and (flag==0)):
                

                """
                    Special care needs to be specified for "active" DNA particles
                    as they will contain excluded volume. Instead of checking for 
                    'A' type particles, we can instead check for is_active_DNA type
                    particles
                """
                            
                if particle_info['is_active_DNA'][monomer] is 1:                    
                    if particle_info['is_active_DNA'][mon] is 1:

                        scaffold = 0;
                        scaffold_mon =0
                        while ((scaffold < pos.shape[0]) and (flag==0)):
                            while((scaffold_mon < pos.shape[0]) and (flag==0)):
                               
                                pos_monomer = pos[scaffold,:] + p;
                                pos_mon = pos[scaffold_mon,:] + position[mon,:] ;

                                min_gap = 0.5* (particle_info['diameter'][monomer]/3+particle_info['diameter'][mon]/3);
                                dist = numpy.power(numpy.sum(numpy.power(pos_monomer-pos_mon,2)),0.5);
                                if dist< min_gap:
                                    flag =1;
                                
                                scaffold_mon = scaffold_mon+1;
                            
                            scaffold = scaffold+1
                    else:
                        scaffold = 0;
                        while ((scaffold < pos.shape[0]) and (flag==0)):
                            pos_monomer = pos[scaffold,:] + p;
                            min_gap = 0.5* (particle_info['diameter'][monomer]/3+particle_info['diameter'][mon]);
                            dist = numpy.power(numpy.sum(numpy.power(position[mon,:] - pos_monomer,2)),0.5);
                            if dist< min_gap:
                                flag =1;
                            scaffold =scaffold+1;
                else:
                    if particle_info['is_active_DNA'][mon] is 1:
                        scaffold_mon =0;
                        while((scaffold_mon < pos.shape[0]) and (flag==0)):
                            pos_mon = pos[scaffold_mon,:] + position[mon,:];
                            min_gap = 0.5* (particle_info['diameter'][monomer]+particle_info['diameter'][mon]/3);
                            dist = numpy.power(numpy.sum(numpy.power(pos_mon - p,2)),0.5);
                            if dist< min_gap:
                                flag =1;
                            
                            scaffold_mon = scaffold_mon+1;

                dist = numpy.power(numpy.sum(numpy.power(position[mon,:] - p,2)),0.5);
                min_gap = 0.5* (particle_info['diameter'][monomer]+particle_info['diameter'][mon]);
                if dist< min_gap:
                    flag =1;    


                mon = mon+1;
            
            flag = flag-1;

            if flag<0:
                position[monomer,:] = p;                    
                monomer = monomer+1;
                chain_pos = chain_pos+1;

        
        while(chain_pos < len(seq_chains)):
        
            mon_type= particle_types.index(seq_chains[chain_pos]);
            particle_info['typeid'][monomer] = mon_type;
            particle_info['diameter'][monomer] = diameter_list[seq_chains[chain_pos]];
            particle_info['mass'][monomer] = mass_list[seq_chains[chain_pos]];
            particle_info['charge'][monomer] = charge_list[seq_chains[chain_pos]];
            particle_info['is_active_DNA'][monomer] = is_active_DNA_list[seq_chains[chain_pos]];

            flag = 0;
            while(flag==0):
                bond_gap = 1.2 * (particle_info['diameter'][monomer]+particle_info['diameter'][monomer-1]);
                p = position[monomer-1,:]+bond_gap*numpy.random.random_sample((1,3)) - 0.5*bond_gap;
                
                if ((p.min()< -L/2) or (p.max() > L/2)):
#                     print("This is not valid {}".format(p));
                    flag=1;
                    
                mon = monomer-1;
                while((mon>=0) and (flag==0)):
                                                
                    if particle_info['is_active_DNA'][monomer] is 1:                    
                        if particle_info['is_active_DNA'][mon] is 1:
 
                            # C.O.M distance between A centers
#                            min_gap = 0.5* (particle_info['diameter'][monomer]+particle_info['diameter'][mon]);
#                            dist = numpy.power(numpy.sum(numpy.power(p-position[mon,:],2)),0.5);
#                            if dist< min_gap:
#                                flag =1;
                           
                           #Now repeat this calculation for the scaffolding particles 
                            
                            
                            scaffold = 0;
                            scaffold_mon =0
                            while ((scaffold < pos.shape[0]) and (flag==0)):
                                while((scaffold_mon < pos.shape[0]) and (flag==0)):
                                   
                                    pos_monomer = pos[scaffold,:] + p;
                                    pos_mon = pos[scaffold_mon,:] + position[mon,:] ;
    
                                    min_gap = 0.5* (particle_info['diameter'][monomer]/3+particle_info['diameter'][mon]/3);
                                    dist = numpy.power(numpy.sum(numpy.power(pos_monomer-pos_mon,2)),0.5);
                                    if dist< min_gap:
                                        flag =1;
                                    
                                    scaffold_mon = scaffold_mon+1;
                                
                                scaffold = scaffold+1
                        else:
                            scaffold = 0;
                            while ((scaffold < pos.shape[0]) and (flag==0)):
                                pos_monomer = pos[scaffold,:] + p;
                                min_gap = 0.5* (particle_info['diameter'][monomer]/3+particle_info['diameter'][mon]);
                                dist = numpy.power(numpy.sum(numpy.power(position[mon,:] - pos_monomer,2)),0.5);
                                if dist< min_gap:
                                    flag =1;
                                scaffold =scaffold+1;
                    else:
                        if particle_info['is_active_DNA'][mon] is 1:
                            scaffold_mon =0;
                            while((scaffold_mon < pos.shape[0]) and (flag==0)):
                                pos_mon = pos[scaffold_mon,:] + position[mon,:];
                                min_gap = 0.5* (particle_info['diameter'][monomer]+particle_info['diameter'][mon]/3);
                                dist = numpy.power(numpy.sum(numpy.power(pos_mon - p,2)),0.5);
                                if dist< min_gap:
                                    flag =1;
                                
                                scaffold_mon = scaffold_mon+1;
    
                
                    dist = numpy.power(numpy.sum(numpy.power(position[mon,:] - p,2)),0.5);
                    min_gap = 0.5* (particle_info['diameter'][monomer]+particle_info['diameter'][mon]);
                    if dist< min_gap:
                        flag =1;    
    
    
                    mon = mon-1;
                
                flag = flag-1;

                
                if flag<0:
#                     print("Dist between {} and {} is {}".format(monomer-1,monomer,numpy.power(numpy.sum(numpy.power(position[monomer-1,:] - p,2)),0.5)))
                  
                    position[monomer,:] = p;
                    particle_info['bonds'][bond_start] = numpy.array([monomer-1,monomer]);
                    bond_type = seq_chains[chain_pos]+'_'+seq_chains[chain_pos-1];
                    if bond_type not in particle_info['bond_types']:
                        particle_info['bond_types'].append(bond_type) ;
                    
                    bond_id = particle_info['bond_types'].index(bond_type);
                    particle_info['bond_id'][bond_start] = bond_id;
                    
                    bond_start = bond_start +1;
                    monomer = monomer+1;
                    chain_pos = chain_pos+1;

        
        mol = mol+1;
        
        
    
    
    particle_info['Nchain_min'] = mol;
    particle_info['N_mono_min'] = monomer;  
    particle_info['bond_start'] =bond_start;
    
    return (particle_info,position)

