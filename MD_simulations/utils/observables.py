#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 13 19:08:26 2018

@author: krishna
"""

import freud
import numpy

class cluster_size_analyze:
    """
        Observable call-back register function to estimate size
        of the largest cluster
    """
    def __init__(self,system,cluster_size=None,cluster_rg=None):
        self.system = system;


    def __call__(self,step):
        snap = self.system.take_snapshot();
        N = self.system.N_total;
        pos =snap.particles.position[0:N,:];
        box = freud.box.Box(snap.box.Lx, snap.box.Ly,snap.box.Lz);
        cluster_g = freud.cluster.Cluster(box,rcut=1.4);
        cluster_g.computeClusters(pos)
        cluster_idx = cluster_g.getClusterIdx();
        cluster_prop = freud.cluster.ClusterProperties(box);
        cluster_prop.computeProperties(pos, cluster_idx)
        a = cluster_prop.getClusterSizes();
#        a_index = numpy.where(a==max(a))
##        rg = cluster_prop.getClusterG();
#        MI_rel = rg[a_index,:,:];
#        eig, eig_val = numpy.linalg.eig(MI_rel);

        return (max(a))

class cluster_rg_analyze:
    """
        Observable call-back register function to estimate radius
        of the largest cluster
    """
    
    
    def __init__(self,system,cluster_size=None,cluster_rg=None):
        self.system = system;

    def __call__(self, step):
        snap = self.system.take_snapshot();
        r_c = (snap.particles.diameter.max());
        N = self.system.N_total;
        pos = snap.particles.position[0:N,:];
        box = freud.box.Box(snap.box.Lx, snap.box.Ly,snap.box.Lz);
        cluster_g = freud.cluster.Cluster(box,rcut=1.4);
        cluster_g.computeClusters(pos)
        cluster_idx = cluster_g.getClusterIdx();
        cluster_prop = freud.cluster.ClusterProperties(box);
        cluster_prop.computeProperties(pos, cluster_idx)
        a = cluster_prop.getClusterSizes();
        a_index = numpy.where(a==max(a))
        rg = cluster_prop.getClusterG();
        MI_rel = rg[a_index,:,:];
#         print(MI_rel)
        eig, eig_val = numpy.linalg.eig(MI_rel);
#         print(numpy.sum(eig))
        return (numpy.sqrt(numpy.sum(eig))+0.5*r_c)

class density_estimate:
    """
        Observable call-back register function to estimate mean
        density of local neighbourhood of particles
    """

    def __init__(self,system,r_cut=3.0):
        self.system = system;
        self.r_cut = r_cut;
    
    def __call__(self,step):
        snap = self.system.take_snapshot();
        N = self.system.N_total;
        print("total number of particles is {}".format(N))
        local_dens = 0.0;
        l_pos = snap.particles.position[0:N,:];
        type_rel= [];
        tag_list = {};
        mult_factor =1;
        for p in numpy.arange(0,N):
            type_of_particle = snap.particles.typeid[p]
            if type_of_particle not in type_rel:
                type_rel.append(type_of_particle);
                tag_list[type_of_particle] = [];
                tag_list[type_of_particle].append(p)
            
            tag_list[type_of_particle].append(p)
            
        for types in tag_list.keys():
            if snap.particles.types[types] =='A':
                mult_factor = 2;
            p = tag_list[types][0];
            volume = mult_factor* 4/3 * numpy.pi* numpy.power(0.5*snap.particles.diameter[p],3);
            diameter = snap.particles.diameter[p]
            gdensity = freud.density.LocalDensity(self.r_cut, volume, diameter)
            fbox = freud.box.Box(snap.box.Lx, snap.box.Ly,snap.box.Lz);         
            gdensity.compute(fbox, l_pos[tag_list[types],:],l_pos)
            y = gdensity.getDensity();
            local_dens = local_dens + numpy.sum(y,axis=0);
        
        return (local_dens/N);
