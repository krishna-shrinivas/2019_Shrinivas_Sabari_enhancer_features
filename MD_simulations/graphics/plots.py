#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 13 18:57:58 2018

@author: krishna
"""

import matplotlib
matplotlib.use('Agg')
import fresnel
import math
import numpy
import gsd
import gsd.fl
import gsd.hoomd
import PIL
import os


device = fresnel.Device(mode='cpu');
path_tracer = fresnel.tracer.Path(device, 300,300)

blue = fresnel.color.linear([0.25,0.5,1])*0.9;
orange = fresnel.color.linear([1.0,0.714,0.169])*0.9


## Function that takes in snapshot, and returns an image 
## with the fundamental units being represented as spheres
def render_sphere_frame(frame, height=None):

    """
        render_sphere_frame takes in a HOOMD system snapshot
        and generates scene of the box, which it returns
        back through a path_tracer sample
        
        This function is an extension of the function provided 
        in the HOOMD tool kit
    
    """
    
    
    if height is None:
        if hasattr(frame, 'configuration'):
            Ly = frame.configuration.box[1]
            height = Ly * math.sqrt(3)
        else:
            Ly = frame.box.Ly;
            height = Ly * math.sqrt(3)

    scene = fresnel.Scene(device)
    scene.lights = fresnel.light.cloudy();
    g = fresnel.geometry.Sphere(scene, position=frame.particles.position, radius=frame.particles.diameter*0.5)
    g.material = fresnel.material.Material(solid=1, color=blue, primitive_color_mix=1.0, specular=1.0, roughness=0.2)
    g.outline_width = 0.07
    scene.camera = fresnel.camera.orthographic(position=(height, height, height), look_at=(0,0,0), up=(0.2,1,0), height=height)

    g.color[frame.particles.typeid == 0] = blue;
    g.color[frame.particles.typeid == 1] = orange;
    g.color[frame.particles.typeid == 2] = [1,0,0];
    g.color[frame.particles.typeid == 3] = [0.5,0.4,0.2];
    g.color[frame.particles.typeid == 4] = [0.2,0.7,0.7];

    scene.background_color = (1,1,1)

    return path_tracer.sample(scene, samples=64, light_samples=20)


## Function that takes input to trajectory file, and outputfile name 
## and saves a high quality gif
def save_movie(frame_gen, gsd_file):
    
    
    """
        save_movie takes in a trajectory file (gsd_file) and a 
        frame method which is typically render_sphere_frame (above)
        
        It generates a movie using the system configuration, which is set by 
        default to 300*300, and saves the movie file for each trajectory under
        a folder called Movies/* in the same folder that the trajectory files
        live. 
    
    """
    f = gsd.fl.GSDFile(gsd_file, 'rb')
    t = gsd.hoomd.HOOMDTrajectory(f)

    a = frame_gen(t[0]);

    if tuple(map(int, (PIL.__version__.split(".")))) < (3,4,0):
        print("Warning! Movie display output requires pillow 3.4.0 or newer.")
        print("Older versions of pillow may only display the first frame.")

    im0 = PIL.Image.fromarray(a[:,:, 0:3], mode='RGB').convert("P", palette=PIL.Image.ADAPTIVE);
    ims = [];
    for f in t[1:]:
        a = frame_gen(f);
        im = PIL.Image.fromarray(a[:,:, 0:3], mode='RGB')
        im_p = im.quantize(palette=im0);
        ims.append(im_p)

    file_full_name = os.path.splitext(gsd_file)[0];
    path_name,file_name = os.path.split(file_full_name);
    f =path_name+'/Movies/'+file_name +'.gif';
    if not os.path.exists(os.path.dirname(f)):
        os.makedirs(os.path.dirname(f),exist_ok=True)

    im0.save(f, 'gif', save_all=True, append_images=ims, duration=1500, loop=0)

    return (f)


## After simulations, generate_plots will plot various logged quantities
## that include temperature, size of largest cluster, RG, PE, KE   
def generate_plots(file_output,f_cluster,file_PE_T_KE):
    
    
    """
        generate_plots takes in the output name of the desired images,
        the data files which are f_cluster and file_PE_T_KE. From these files,
        the function plots the cluster size, cluster RG, fraction of largest cluster,
        the PE, Temp, and the KE, and stores them with uniq identifiers
        for each trajectory.
    
    """
    
    from matplotlib import pyplot
    #get_ipython().magic('matplotlib inline')
    pyplot.figure(figsize=(4,2.2), dpi=140);
    data = numpy.genfromtxt(fname=f_cluster, skip_header=True);
    pyplot.plot(data[:,0]/1e5, data[:,1]);
    pyplot.xlabel('time step (* 10^5)');
    pyplot.ylabel('S(t)');
    pyplot.savefig(file_output+'_S.png',bbox_inches='tight')


    pyplot.figure(figsize=(4,2.2), dpi=140);
    pyplot.plot(data[:,0]/1e5, (data[:,2]) );
    pyplot.xlabel('time step (* 10^5)');
    pyplot.ylabel('Rg(t)');
    pyplot.savefig(file_output+'_Rg.png',bbox_inches='tight')



    pyplot.figure(figsize=(4,2.2), dpi=140);
    vol = 4/3 * numpy.pi * pow((data[:,2]),3)
    pyplot.plot(data[:,0]/1e5, numpy.divide(data[:,1],vol) );
    pyplot.xlabel('time step (* 10^5)');
    pyplot.ylabel('Rg(t)');
    pyplot.savefig(file_output+'_frac.png',bbox_inches='tight')


    # Use matplotlib to plot the potential energy vs time step.

    # In[15]:


    # #get_ipython().magic('matplotlib inline')
    # pyplot.figure(figsize=(4,2.2), dpi=140);
    # data = numpy.genfromtxt(fname=file_g_r_output, skip_header=False);
    # pyplot.plot(data[:,0], data[:,1]);
    # pyplot.xlabel('time step');
    # pyplot.ylabel('g(r)');


    # In[16]:


    data = numpy.genfromtxt(fname=file_PE_T_KE, skip_header=True);
    pyplot.figure(figsize=(4,2.2), dpi=140);
    pyplot.plot(data[:,0]/1e5, data[:,1]);
    pyplot.xlabel('time step (* 10^5)');
    pyplot.ylabel('potential_energy');
    pyplot.savefig(file_output+'_PE.png',bbox_inches='tight')


    # In[17]:


    pyplot.figure(figsize=(4,2.2), dpi=140);
    pyplot.plot(data[:,0]/1e5, data[:,2]);
    pyplot.xlabel('time step (* 10^5)');
    pyplot.ylabel('temperature');
    pyplot.savefig(file_output+'_T.png',bbox_inches='tight')
    # In[17]:


    pyplot.figure(figsize=(4,2.2), dpi=140);
    pyplot.plot(data[:,0]/1e5, data[:,3]);
    pyplot.xlabel('time step (* 10^5)');
    pyplot.ylabel('Kinetic energy');
    pyplot.savefig(file_output+'_KE.png',bbox_inches='tight')
