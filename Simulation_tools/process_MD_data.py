#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Created on Sun Jul  1 12:59:38 2018
#@author: krishna-shrinivas

"""
    Module contains functions that aid in peforming post-processing analysis,
    computations, and generating visualization from output of simulation data.
"""

import PIL
import json
import numpy
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from os.path import isfile, join
import csv
import matplotlib; mpl =matplotlib
import fresnel
import gsd
import gsd.fl
import gsd.hoomd
import freud

mpl.rcParams['figure.figsize'] = (12.0,8.0) # default = (6.0, 4.0)
mpl.rcParams['font.size']      = 24         # default = 10

mpl.rcParams['axes.linewidth']    = 2.0 # default = 1.0
mpl.rcParams['lines.linewidth']   = 2.0 # default = 1.0
mpl.rcParams['patch.linewidth']   = 1.0 # default = 1.0
mpl.rcParams['grid.linewidth']    = 1.0 # default = 0.5
mpl.rcParams['xtick.major.width'] = 1.0 # default = 0.5
mpl.rcParams['xtick.minor.width'] = 1.0 # default = 0.5
mpl.rcParams['ytick.major.width'] = 1.0 # default = 0.5
mpl.rcParams['ytick.minor.width'] = 1.0 # default = 0.5


class NumpyEncoder(json.JSONEncoder):
    """ Special json encoder for numpy types """
    def default(self, obj):
        if isinstance(obj, (numpy.int_, numpy.intc, numpy.intp, numpy.int8,
            numpy.int16, numpy.int32, numpy.int64, numpy.uint8,
            numpy.uint16, numpy.uint32, numpy.uint64)):
            return int(obj)
        elif isinstance(obj, (numpy.float_, numpy.float16, numpy.float32,
            numpy.float64)):
            return float(obj)
        elif isinstance(obj,(numpy.ndarray,)): #### This is the fix
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


def make_nice_axis(ax):
    """ Function to beautify axis, based on version written by D. Murakowski"""

    ax.spines['top'].set_visible(False) # hide top axs
    #ax.spines['right'].set_position(('outward', 30))
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_position(('outward', 20))
    ax.spines['left'].set_position(('outward', 30))
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.xaxis.set_tick_params(pad=10)
    ax.yaxis.set_tick_params(pad=10)
    ax.xaxis.labelpad = 10
    ax.yaxis.labelpad = 20


def generate_average_plot(mypath,observable,pos,show_plot=False,return_ind_traj=False):

    """
    Function takes in path to observables data, and meta-data to describe it
    and uses it generate visualization. This function largely is called
    by plot_observable.

    **Input variables**

    -   mypath = Path to Data/Observables/,
    -   observable = the unique substring observable to be queried. Options are 'PE' or 'cluster'
    -   pos = Column position to aggregate over trajectory repeats. Refer standard output formatting below to identify.
        Columns (1,2,3) in cluster files are

        1. Cluster size
        2. Cluster radius
        3. Normalized density

        Columns (1,2,3) in PE files are

        1. Potential energy
        2. Kinetic energy
        3. Temperature

        Column 0 in all files are the time-step

    -   show_plot = Flag to visualize mean of observable (Default is FALSE)
    -   return_ind_traj =  Flag to return observable data for individual trajectories

    **Function returns**:

    -   cluster_mean =  Numpy matrix of N_timepoints containing the mean observable size

    -   time =  Numpy array of size N_timepoints*1 containing time-step values

    -   header = String of header of observable file

    -   cs =    Numpy matrix of size N_trajectories * N_timepoints with observable data for each trajectory.

    """
    onlyfiles = sorted([f for f in os.listdir(mypath) if isfile(join(mypath, f))])

    cluster_size = [];
    count=0;
    t= [];

    """ New check added to ensure that only trajectories that all converged
    will contribute to averaging """
    size_of_traj = [];

    for f in onlyfiles:
        if observable in f:
            data = numpy.genfromtxt(fname=join(mypath,f), skip_header=True);
            if data.shape[0]>=10:
                t.append(data[:,0]);
                cluster_size.append(data[:,pos]);
                size_of_traj.append(len(data[:,pos]));
                count = count+1;
                if count>=1:
                    header_file = open(join(mypath,f), 'r');
                    reader = csv.reader(header_file,delimiter='\t');
                    header = next(reader);

    longest_run = max(size_of_traj);
    cs = [];
    for d in numpy.arange(len(cluster_size)):
        if len(cluster_size[d]) == longest_run:
            cs.append(cluster_size[d]);
            if len(t[d]) == longest_run:
                time = t[d];


    cs = numpy.reshape(cs,[len(cs),len(cs[0])]);
    cluster_mean = numpy.mean(cs,0);

    if show_plot:
        plt.figure(figsize=(4,2.2), dpi=140);

        traj_color = '#d3d3d3'
        plt.plot(data[:,0]/1e5, cs.T,c=traj_color,lw=1);
        plt.plot(data[:,0]/1e5, cluster_mean,c='green',lw=2);
        plt.xlabel('time step (* 10^5)');
        plt.ylabel(header[pos]);
        plt.show()

#    print(len(time), cs.shape)
    return (cluster_mean,time,header,cs) if return_ind_traj else (cluster_mean,time,header)


def plot_observable(param_name,param_range,unique_paths,input_params,observable,pos,loglog=True,file_save=False,return_ind_traj=False,scaled=False):

    """
    Function takes in input variables related to set of conditions of MD
    data (obtained typically from parse_MD_output) and generates
    beautified plots of mean_observable vs time over a range of conditions.

    **Input parameters**

    -   param_name =    string name of the variable that is swept over.

    -   param_range ==  sorted list of the different values of parameters

    -   unique_paths == list of paths to each observable folder for different conditions in the output data-set (e.g. list of folder/condition/Data/Observables
        for each condition)

    -   input_params == list of input_params for each condition

    -   observable = the unique substring observable to be queried
        Options are PE or cluster

    -   pos = Column position to aggregate over trajectory repeats.
        Refer standard output formatting below to indentify.

        Columns (1,2,3) in cluster files are

        1. Cluster size
        2. Cluster radius
        3. Normalized density

        Columns (1,2,3) in PE files are

        1. Potential energy
        2. Kinetic energy
        3. Temperature

        Column 0 in all files are the time-step

    -   loglog = Flag to generate axes in log-log scales (Default is True)

    -   file_save = Deprecated flag. Has no function!
        (To-do: Remove variable)

    -   scaled =    Flag to return observable data scaled by number of binding
        sites. Only set to TRUE when calculating for explicit-IDR
        systems with observable = 'cluster' and pos ='1'

    -   return_ind_traj =  Flag to return observable data for individual trajectories

    **Output variables returned**

    -   mean_data = Numpy matrix of N_params *N_timepoints containing the mean
        observable size for each parameter condition

    -   time =  Numpy array of size N_timepoints*1 containing time-step values

    -   all_traj = List of numpy matrices, each of size Ntraj*N_timepoints

    -   axes =  Axes variable/plot handle for generated graph
    """


    dt = float(input_params[0]['dt']);

    if param_name.find('seq') > -1:
        a = [len(x) for x in param_range];
        p = numpy.argsort(a)
    else:
        p = numpy.argsort(param_range)



    binding_size = [];
    mean_data = [];
    legend_labels = [];
    if return_ind_traj:
        all_traj = {};

    for count in numpy.arange(len(param_range)):

        if ((input_params[p[count]]['N_A']) and (scaled)):
            if input_params[p[count]]['seq_A'].count('A'):
                binding_size.append(input_params[p[count]]['N_A']*(float((input_params[p[count]]['seq_A'].count('A'))))*(1+input_params[p[count]]['N_bs_AB']+input_params[p[count]]['N_bs_AC']));
            else:
                binding_size.append(1.0)
        else:
            binding_size.append(1.0);

        mypath = unique_paths[p[count]];
        if not return_ind_traj:
            (cl_mean,time,header) = generate_average_plot(mypath,observable,pos);
        else:
            (cl_mean,time,header,cl_all) = generate_average_plot(mypath,observable,pos,return_ind_traj=True);
            all_traj[count] = cl_all/binding_size[count];


        mean_data.append(cl_mean/binding_size[count]);
        if param_name.find('seq') == -1:
            legend_labels.append( param_name + ' = '+ str(round(param_range[p[count]],2)));
        else:
            legend_labels.append( param_name + ' $ _{l} $ = '+ str(len(param_range[p[count]])));

    mean_data = numpy.reshape(mean_data,(len(param_range),len(time)))
    colors = iter(cm.afmhot(numpy.linspace(0, 0.5, len(param_range))))

    fig, axes = plt.subplots(1,1,sharex=True)
    make_nice_axis(axes)
    L_interest = len(time)- 1;
    for count in numpy.arange(len(param_range)):
        if loglog:
            axes.loglog(dt*time[0:L_interest],mean_data[count,0:L_interest],color=next(colors),lw=4);
        else:
            axes.plot(dt*time[0:L_interest],mean_data[count,0:L_interest],color=next(colors),lw=4);

    axes.set_ylabel(header[pos])
    axes.set_xlabel('Time')

    # axes.vlines(5000,min(cluster_mean[0,:]/binding_size),max(cluster_mean[0,:]/binding_size),lw=2,linestyle='--',color="Grey");
    # axes.hlines(1,[0],max(dt*time),lw=2,linestyle='--',color="Grey");

    axes.legend(legend_labels,bbox_to_anchor=(1.05, 1),fontsize =20)

    return (time,mean_data,axes) if not return_ind_traj else (time,mean_data,axes,all_traj)




def plot_window_observable(param_name,param_range,unique_paths,input_params,observable,pos,window_start,window_end,return_ind_traj=False,scaled=False,account_DNA=False):

    """
    Function generates plots of observable over specific time-windows on
    MD output data over a set of conditions.

    **Input parameters**

    -   param_name =   String name of the variable that is swept over

    -   param_range ==  sorted list of the different values of parameters

    -   unique_paths == list of paths to each observable folder for different
        conditions in the output data-set
        (e.g. list of folder/condition/Data/Observables for each condition)

    -   input_params == list of input_params for each condition

    -   observable = the unique substring observable to be queried
        Options are PE or cluster

    -   pos = Column position to aggregate over trajectory repeats.
        Refer standard output formatting below to indentify.

        Columns (1,2,3) in cluster files are

        1. Cluster size
        2. Cluster radius
        3. Normalized density

        Columns (1,2,3) in PE files are

        1. Potential energy
        2. Kinetic energy
        3. Temperature

        Column 0 in all files are the time-step

    -   window_start = Index of time_point to start calculation

    -   window_end =    Index of time_point to end calculation

    -   scaled =    Flag to return observable data scaled by number of binding
        sites. Only set to TRUE when calculating for explicit-IDR
        systems with observable = 'cluster' and pos ='1'

        DEFAULT: FALSE

    -   return_ind_traj =   Flag to return variance of observable data for individual trajectories
        DEFAULT: FALSE

    -   account_DNA =   Variable that calculates an alternate metric of scaled
        size where the size of DNA is explicitly subtracted out.

        DEFAULT: FALSE

    **Output variables returned**

    -   mean_all =list of size N_params containing the mean
        observable in the window size for each parameter condition

    -   axes =  Axes variable/plot handle for generated graph

    -   std_all =   list of standard deviance in observable value over
        window period

    """


    dt = float(input_params[0]['dt']);

    if param_name.find('seq') > -1:
        a = [len(x) for x in param_range];
        p = numpy.argsort(a)
    else:
        p = numpy.argsort(param_range)



    binding_size = [];
    binding_sites = [];

    mean_data = [];
    legend_labels = [];
    if return_ind_traj:
        all_traj = {};

    for count in numpy.arange(len(param_range)):

        if ((input_params[p[count]]['N_A']) and (scaled)):
            if input_params[p[count]]['seq_A'].count('A'):
                binding_size.append(input_params[p[count]]['N_A']*(float((input_params[p[count]]['seq_A'].count('A'))))*(1+input_params[p[count]]['N_bs_AB']+input_params[p[count]]['N_bs_AC']));
                binding_sites.append(input_params[p[count]]['N_A']*(float((input_params[p[count]]['seq_A'].count('A'))))*(input_params[p[count]]['N_bs_AB']+input_params[p[count]]['N_bs_AC']));
            else:
                binding_size.append(1.0)
                binding_sites.append(1.0)

        else:
            binding_size.append(1.0);
            binding_sites.append(1.0)

        if (binding_sites[count] < 1.0):
            binding_sites[count] = input_params[p[count]]['seq_A'].count('A');

        len_A =len(input_params[p[count]]['seq_A']);
        mypath = unique_paths[p[count]];
        if not return_ind_traj:
            (cl_mean,time,header) = generate_average_plot(mypath,observable,pos);
        else:
            (cl_mean,time,header,cl_all) = generate_average_plot(mypath,observable,pos,return_ind_traj=True);
            if not account_DNA:
                all_traj[count] = cl_all/binding_size[count];
            else:
                all_traj[count] = (cl_all-len_A)/binding_sites[count];

        if not account_DNA:
            mean_data.append(cl_mean/binding_size[count]);
        else:
            mean_data.append((cl_mean-len_A)/binding_sites[count]);

    mean_data = numpy.reshape(mean_data,(len(param_range),len(time)))

    fig, axes = plt.subplots(1,1,sharex=True)
    make_nice_axis(axes)
    L_interest = len(time)- 1;
    mean_all = [];
    std_all = [];
    for count in numpy.arange(len(param_range)):
        mean = numpy.mean(mean_data[count,window_start:window_end]);
        mean_all.append(mean);
        if return_ind_traj:
            std = numpy.std((all_traj[count][:,window_start:window_end]));
            std_all.append(std);

    mean_all = numpy.array(mean_all)
    std_all = numpy.array(std_all)
    if return_ind_traj:
        axes.plot(param_range,mean_all,lw =3,color='black');
        axes.fill_between(param_range,mean_all-std_all,mean_all+std_all,facecolor='grey',alpha=0.1)

    else:
        axes.scatter(param_range,mean_all);
    axes.set_ylabel(header[pos])
    axes.set_xlabel(param_name)

    # axes.vlines(5000,min(cluster_mean[0,:]/binding_size),max(cluster_mean[0,:]/binding_size),lw=2,linestyle='--',color="Grey");
    # axes.hlines(1,[0],max(dt*time),lw=2,linestyle='--',color="Grey");

    return (mean_all,axes) if not return_ind_traj else (mean_all,std_all,axes)





def estimate_chain_volume(seq,input_params):

    """
    Estimate the total volume occupied by the monomers in any given
    sequence.

    Takes in sequence string, input_params and returns volume.
    """
    volume = 0;
    for aa in seq:
        if aa is not 'A':
            volume = volume + 4/3*numpy.pi*numpy.power(input_params['diameter'+aa]/2,3);
        else:
            volume = volume + numpy.power(input_params['diameter'+aa],3);

    return volume;




def get_volume_fraction(input_params):
    """
    Estimate the volume fraction of DNA (type A), TF (type B), and coa (type C)

    Takes in input_parameters and returns volume fractions.
    """

    vol_A =input_params['N_A']*(estimate_chain_volume(input_params['seq_A'],input_params))/float(pow(input_params['L'],3));
    vol_B = input_params['N_B']*(estimate_chain_volume(input_params['seq_B'],input_params))/float(pow(input_params['L'],3));
    vol_C = input_params['N_C']*(estimate_chain_volume(input_params['seq_C'],input_params))/float(pow(input_params['L'],3));

    return (vol_A,vol_B,vol_C)



def calculate_entropy_enthalpy(param_name,param_range,unique_paths,input_params):

    """
    Function generates estimates of entropy and enthalpy over list of different
    parameter conditions of input MD data.

    **Input parameters**

    -   param_name =    The string name of the variable that is swept over.

    -   param_range =   sorted list of the different values of parameters

    -   unique_paths = list of paths to each observable folder for different
        conditions in the output data-set
        (e.g. list of folder/condition/Data/Observables for each condition)

    -   input_params = list of input_params for each condition

    **Output variables returned**

    -   PE_all_traj =   list of size N_params containing the Ntraj*Ntimepoints
        matrices of potential energy for all trajectories

    -   entropy_all_traj =  list of size N_params containing the Ntraj*Ntimepoints
        matrices of entropy for all trajectories

    -   scaled_cluster_size_all_traj =  list of size N_params containing the Ntraj*Ntimepoints
        matrices of scaled cluster size for all trajectories

    -   time*dt =  Numpy array of size N_timepoints*1 containing times

    -   PE_store_mean = list of size N_params *Ntimepoints matrices of <potential energy>

    -   entropy_mean = list of size N_params *Ntimepoints matrices of <entropy>

    -   scaled_cs_mean = list of size N_params *Ntimepoints matrices of <cluster size>

    """

    dt = float(input_params[0]['dt']);

    binding_size = [];
    volume_fraction = numpy.zeros((len(param_range),3))
    PE_all_traj ={};
    entropy_all_traj ={};
    entropy_mean = [];
    scaled_cs_mean = [];
    PE_store_mean = [];
    scaled_cluster_size_all_traj = {};
    count = 0;
    for param in numpy.arange(len(param_range)):
        mypath = unique_paths[param];

        (cl_mean,time,header,cl_all) = generate_average_plot(mypath,'cluster',1,return_ind_traj=True);
        (rg_mean,time,header,rg_all) = generate_average_plot(mypath,'cluster',2,return_ind_traj=True);
        (PE_mean,time,header,PE_all) = generate_average_plot(mypath,'PE',1,return_ind_traj=True);
        L_interest = int( numpy.floor(len(time)/1));

        (vol_A,vol_B,vol_C) = get_volume_fraction(input_params[param]);
        volume_fraction[param,0] = vol_A;
        volume_fraction[param,1] = vol_B;
        volume_fraction[param,2] = vol_C;
        if input_params[param]['N_A']:
            if input_params[param]['seq_A'].count('A'):
                binding_size.append(input_params[param]['N_A']*(float((input_params[param]['seq_A'].count('A'))))*(1+input_params[param]['N_bs_AB']+input_params[param]['N_bs_AC']));
            else:
                binding_size.append(1.0)
            size_A = int(len(input_params[param]['seq_A']))
        else:
            binding_size.append(1.0);
            size_A = 1;

        volume_rest = float(pow(input_params[count]['L'],3)) * (1-numpy.sum(volume_fraction[count,:])) ;

        entropy = numpy.multiply(cl_all[:,0:L_interest]-size_A+1,numpy.log(numpy.divide(4/3*numpy.pi*numpy.power(rg_all[:,0:L_interest],3),volume_rest)));
        entropy_mean.append(numpy.mean(entropy,0));
        scaled_cs_mean.append(numpy.mean(numpy.divide(cl_all[:,0:L_interest],binding_size[count]),axis=0));
        PE_store_mean.append(PE_mean)
        PE_all_traj[count] = PE_all;
        entropy_all_traj[count] = entropy;
        scaled_cluster_size_all_traj[count] = numpy.divide(cl_all[:,0:L_interest],binding_size[count])
        count =count+1;

    PE_store_mean = numpy.reshape(PE_store_mean,(len(param_range),len(time)))

    entropy_mean = numpy.reshape(entropy_mean,(len(param_range),len(time)))
    scaled_cs_mean = numpy.reshape(scaled_cs_mean,(len(param_range),len(time)))

    return(PE_all_traj,entropy_all_traj,scaled_cluster_size_all_traj,time*dt,entropy_mean,PE_store_mean,scaled_cs_mean)

def plot_entropy_enthalpy(param_name,param_range,unique_paths,input_params):
    """
    Function generates plots of entropy and enthalpy over list of different
    parameter conditions of input MD data

    **Input parameters**

    -   param_name =    The string name of the variable that is swept over.

    -   param_range =  sorted list of the different values of parameters

    -   unique_paths = list of paths to each observable folder for different
        conditions in the output data-set
        (e.g. list of folder/condition/Data/Observables for each condition)

    -   input_params == list of input_params for each condition

    **Output variables returned**

    -   PE_mean = list of size N_params *Ntimepoints
        matrices of <potential energy>

    -   entropy_store = list of size N_params *Ntimepoints
        matrices of <entropy>

    -   cluster_mean = list of size N_params *Ntimepoints
        matrices of <cluster size>

    -   rg_mean = list of size N_params *Ntimepoints
        matrices of <cluster radius>

    -   axes =  Axes variable/plot handle for generated graph

        NOTE: The entropy calculation is performed for the largest cluster
        by estimating the loss in "Free volume" for all individual
        monomers.

        For polymers, this is in general an over-estimate of the entropy
        loss. Rather, a unit of kT log(V1/V2) should be added only
        for each chain that joins the largest cluster. These calculations
        are done separately now, employing the cluster_composition
        model.

    """
    dt = float(input_params[0]['dt']);

    cluster_mean = [];
    binding_size =  [];

    rg_mean = [];
    PE_mean = [];
    volume_fraction = numpy.zeros((len(param_range),3))
    for param in numpy.arange(len(param_range)):
        mypath = unique_paths[param];
        (cl_mean,time,header) = generate_average_plot(mypath,'cluster',1);
        cluster_mean.append(cl_mean);
        (rg,time,header) = generate_average_plot(mypath,'cluster',2);
        rg_mean.append(rg);
        (PE,time,header) = generate_average_plot(mypath,'PE',1);
        PE_mean.append(PE)
        (vol_A,vol_B,vol_C) = get_volume_fraction(input_params[param]);
        volume_fraction[param,0] = vol_A;
        volume_fraction[param,1] = vol_B;
        volume_fraction[param,2] = vol_C;
        if input_params[param]['N_A']:
            if input_params[param]['seq_A'].count('A'):
                binding_size.append(input_params[param]['N_A']*(float((input_params[param]['seq_A'].count('A'))))*(1+input_params[param]['N_bs_AB']+input_params[param]['N_bs_AC']));
            else:
                binding_size.append(1.0)
            size_A = int(len(input_params[param]['seq_A']))
        else:
            binding_size.append(1.0);
            size_A = 1;
    cluster_mean = numpy.reshape(cluster_mean,(len(param_range),len(time)))
    rg_mean = numpy.reshape(rg_mean,(len(param_range),len(time)))
    PE_mean = numpy.reshape(PE_mean,(len(param_range),len(time)))

    L_interest = int( numpy.floor(len(time)/1));
    fig, axes = plt.subplots(len(param_range),2,sharex=True,figsize=(12*len(param_range),9*len(param_range)))
    legend_labels = [];
    entropy_store = [];
    for count in numpy.arange(len(param_range)):
        make_nice_axis(axes[count,1]);
        make_nice_axis(axes[count,0]);

        axes[count,1].plot(dt*time[0:L_interest],PE_mean[count,0:L_interest],color='black',lw=4,label = 'H');

        volume_rest = float(pow(input_params[count]['L'],3)) * (1-numpy.sum(volume_fraction[count,:])) ;
        entropy = numpy.multiply(cluster_mean[count,0:L_interest],numpy.log(numpy.divide(4/3*numpy.pi*numpy.power(rg_mean[count,0:L_interest],3),volume_rest)));
        entropy_store.append(entropy)
        axes[count,1].plot(dt*time[0:L_interest],entropy,color='green',lw=4, label = ' $\delta$S, ');
        if param_name.find('seq') == -1:
            legend_labels.append( param_name + ' = '+ str(round(param_range[count],2)));
        else:
            legend_labels.append( param_name + ' $ _{l} $ = '+ str(len(param_range[count])));

        axes[count,1].legend(fontsize =20)



        axes[count,0].plot(dt*time[0:L_interest],cluster_mean[count,0:L_interest]/binding_size[count],color='black',lw=4,label = 'Scaled cluster size');
        axes[count,0].set_ylabel('Scaled size')
        axes[count,0].set_ylim(0,max(cluster_mean[count,0:L_interest]/binding_size[count]))
        axes[count,0].hlines(1,0,max(dt*time[0:L_interest]),lw=2,linestyle='--',color="Grey")
        axes[count,0].legend([legend_labels[count],'Stoich'],fontsize =20)

        axes[count,1].set_ylabel('Energy (kT)')


    axes[count,1].set_xlabel('Time')
    axes[count,0].set_xlabel('Time')

    entropy_store = numpy.reshape(entropy_store,(len(param_range),len(time)))

    return(cluster_mean,rg_mean,PE_mean,entropy_store,axes)



class cluster_composition():

    """

    Module computes detailed statistics on cluster compositions of different
    types of molecules at all times + conditions.

    """
    def __init__(self):
        self.cluster_distribution = {};
        self.N_tot_mean = {};
    def convert_snap_cluster_dist(self,snap,input_params):
        """
        Function generates detailed cluster statistics from a MD snapshot

        **Input parameters**

        -   snap =  MD frame from simulation

        -   input_params =  Input parameters for simulations

        **Output variables returned**

        -   c_max =     Size of the largest cluster

        -   rg_max =    Radius of the largest cluster

        -   count_types =   Dictionary of types of chains (DNA,TF,coa - A,B,C)
            in largest cluster

        -   chain_list  = IDs of chains in largest cluster

        -   monomers_in_largest_cluster[0] = IDs of monomers in largest cluster

        """
        N = input_params['N_monomers']
        monomer_to_chain_map = input_params['MC_map'];
        chain_type = input_params['CT_map'];


        pos =snap.particles.position[0:N,:];
        box = freud.box.Box(snap.configuration.box[0],snap.configuration.box[1],snap.configuration.box[2]);
        cluster_g = freud.cluster.Cluster(box,rcut=1.4);
        cluster_g.computeClusters(pos)
        cluster_idx = cluster_g.getClusterIdx();
        cluster_prop = freud.cluster.ClusterProperties(box);
        cluster_prop.computeProperties(pos, cluster_idx)
        a = cluster_prop.getClusterSizes();
        rg = cluster_prop.getClusterG();
        all_a_index = numpy.where(a==max(a));
        a_index = all_a_index[0][0];
        monomers_in_largest_cluster = numpy.where(cluster_idx==a_index)
    #    print(monomers_in_largest_cluster[0])
        count_types = {};
        chain_list = {};
        count_types['A'] =0;
        count_types['B'] =0;
        count_types['C'] =0;

        for monomer in monomers_in_largest_cluster[0]:
            chain_id = monomer_to_chain_map[monomer]

            if str(chain_id) not in chain_list.keys():
                ctype = chain_type[chain_id];
                chain_list[str(chain_id)] = 1;
                count_types[ctype] +=1;

        MI_rel = rg[a_index,:,:];
        eig, eig_val = numpy.linalg.eig(MI_rel);
        rg_max = numpy.sqrt(numpy.sum(eig))+0.5*max(snap.particles.diameter);
        c_max = a[a_index];
        return (c_max,rg_max,count_types,chain_list,monomers_in_largest_cluster[0])

    def generate_cluster_distribution(self,param_name,param_range,unique_paths,input_params,file_save = False):
        """
        Function generates detailed cluster statistics for MD simulations across
        all conditions and times under output Folder/* , and writes output
        to Folder../Data/Analysis/ for each condition.

        Each trajectory will have a text file written with detailed
        statistics of output.

        **Input parameters**

        -   param_name =   The string name of the variable that is swept over.

        -   param_range =  sorted list of the different values of parameters

        -   unique_paths = list of paths to each observable folder for different
            conditions in the output data-set
            (e.g. list of folder/condition/Data/Observables for each condition)

        -   input_params = list of input_params for each condition

        **Output variables returned**

        -   cluster_distribution = Dictionary of entire cluster distribution, with the following
            nested key structure.

            1. Path to condition trajectories

            2. Unique replicate traj ID

            3. Time step of cluster statistics

            4. Cluster statistics for that relevant time-step

        """
        count = 0;
#        cluster_distribution = {};
        for path in unique_paths:
            traj_path = path.replace('Observables','Trajectory');
            ana_path  = path.replace('Observables','Analysis');
            ll = ana_path.split('/');
            for p in numpy.arange(1,len(ll)):
                try: os.makedirs('/'.join(ll[0:p]))
                except OSError: pass

            p =sorted(os.listdir(traj_path));
            self.cluster_distribution[traj_path] = {};
            if p[0] == 'Movies':
                del(p[0])
            if p[-1] == 'Movies':
                del(p[-1])
            for file in p:
                self.cluster_distribution[traj_path][file]= {};
                gsd_file =traj_path + file;
                f = gsd.fl.GSDFile(gsd_file, 'rb')
                t = gsd.hoomd.HOOMDTrajectory(f)
                for snap in t:
                    self.cluster_distribution[traj_path][file][str(snap.configuration.step)] = {};
                    (self.cluster_distribution[traj_path][file][str(snap.configuration.step)]['c_max'],self.cluster_distribution[traj_path][file][str(snap.configuration.step)]['rg_max'],self.cluster_distribution[traj_path][file][str(snap.configuration.step)]['count_types'],self.cluster_distribution[traj_path][file][str(snap.configuration.step)]['chain_list'],self.cluster_distribution[traj_path][file][str(snap.configuration.step)]['monomer_list']) = self.convert_snap_cluster_dist(snap,input_params[count]);

                output_file = gsd_file.replace('Trajectory','Analysis');
                output_file = output_file.replace('gsd','text');
                with open(output_file,'w') as f:
                    f.write(json.dumps(self.cluster_distribution[traj_path][file],cls=NumpyEncoder));

            count = count+1;
        self.param_range = param_range;
        self.param_name = param_name;
        self.unique_paths = unique_paths;
        self.input_params = input_params;
        return (self.cluster_distribution)


    def parse_cluster_composition(self,param_name,param_range,unique_paths,input_params):

        """
        Function reads detailed statistics that are written from generate_cluster_distribution

        Each trajectory will have a text file written with detailed
        statistics of output.

        **Input parameters**

        -   param_name =   The string name of the variable that is swept over.

        -   param_range =  sorted list of the different values of parameters

        -   unique_paths = list of paths to each observable folder for different
            conditions in the output data-set
            (e.g. list of folder/condition/Data/Observables for each condition)

        -   input_params = list of input_params for each condition

        **Output variables returned**

        -   cluster_distribution = Dictionary of entire cluster distribution, with the following
            nested key structure.

            1. Path to condition trajectories

            2. Unique replicate traj ID

            3. Time step of cluster statistics

            4. Cluster statistics for that relevant time-step

        """

        count = 0;
#        cluster_distribution = {};
        for path in unique_paths:
            traj_path = path.replace('Observables','Trajectory');
            p =sorted(os.listdir(traj_path));
            self.cluster_distribution[traj_path] = {};
            if p[0] == 'Movies':
                del(p[0])
            if p[-1] == 'Movies':
                del(p[-1])
            for file in p:
                self.cluster_distribution[traj_path][file]= {};
                gsd_file =traj_path + file;
                analysis_file = gsd_file.replace('Trajectory','Analysis');
                analysis_file = analysis_file.replace('gsd','text');
                with open(analysis_file,'r') as f:
                    self.cluster_distribution[traj_path][file] = json.load(f);

            count = count+1;
        self.param_range = param_range;
        self.param_name = param_name;
        self.unique_paths = unique_paths;
        self.input_params = input_params;
        return(self.cluster_distribution)

    def plot_cluster_composition(self,scaled=False,no_legend=False,ind_plot=False):

        """
        Function reads detailed statistics that are written from generate_cluster_distribution

        Each trajectory will have a text file written with detailed
        statistics of output.

        **Input parameters**

        -   no_legend = Flag to toggle legend ON/OFF

        -   scaled =    Flag to return observable data scaled by number of binding sites.

        -   ind_plot = Deprecated flag - Needs to be removed. ONLY DEFAULT VALUE OF FALSE.

        **Output variables returned**

        -   N_tot_mean =    Numpy array of <cluster_size> of size timepoints for last computed condition

        -   time_points =   Numpy array of time steps

        -   axes =          Handle for generated plot

        -   store_all_means=    Numpy matrix of size N_params * N_timepoints with
            <cluster size> for each trajectory.

        """

        colors = list(cm.afmhot(numpy.linspace(0, 0.5, len(self.param_range))));

        if ~ind_plot:
            fig, axes = plt.subplots(1,1,sharex=True)
            make_nice_axis(axes);
        else:
            fig, axes = plt.subplots(2,2,sharex=True)
            make_nice_axis(axes[0]) ;
            make_nice_axis(axes[1]);
            make_nice_axis(axes[2]);
            make_nice_axis(axes[3]);

        files = list(self.cluster_distribution.keys());
        count_file = 0;
        binding_size = [];
        store_all_means =[];
        self.N_tot = {};
        for file in files:
            if ((self.input_params[count_file]['N_A']) and (scaled)):
                if self.input_params[count_file]['seq_A'].count('A'):
                    binding_size.append(self.input_params[count_file]['N_A']*(float((self.input_params[count_file]['seq_A'].count('A'))))*(self.input_params[count_file]['N_bs_AB']+self.input_params[count_file]['N_bs_AC']));
                else:
                    binding_size.append(1.0)
            else:
                binding_size.append(1.0);

            self.N_tot[count_file] = {};
            traj = list(self.cluster_distribution[file].keys());

            size_of_traj = [];
            for tr in traj:
                size_of_traj.append((len(self.cluster_distribution[file][tr].keys())))
            longest_size = max(size_of_traj);

            Na = [];
            Nb = [];
            Nc = [];
            Ntot = [];
            count = 0;
            for tr in traj:
                if (len(self.cluster_distribution[file][tr].keys()) == longest_size):

                    self.N_tot[count_file][count] = [];

                    time_points = [int(x) for x in self.cluster_distribution[file][tr].keys()];

                    for t in self.cluster_distribution[file][tr].keys():
                        Nc.append([]);
                        Nb.append([]);
                        Na.append([]);
                        Ntot.append([]);
                        Nc[count].append((self.cluster_distribution[file][tr][t]['count_types']['C']))
                        Nb[count].append((self.cluster_distribution[file][tr][t]['count_types']['B']))
                        Na[count].append((self.cluster_distribution[file][tr][t]['count_types']['A']))
                        Ntot[count].append(sum(self.cluster_distribution[file][tr][t]['count_types'].values()));

                    if ~ind_plot:
                        axes.plot(time_points,numpy.array(Ntot[count])/float(binding_size[count_file]),color='grey',lw=0.1);

                    self.N_tot[count_file][count] = numpy.array(Ntot[count])/float(binding_size[count_file]);
                    count = count+1;
            if ~ind_plot:
                N_tot_mean = numpy.array(Ntot[0]);
                for p in numpy.arange(1,len(self.N_tot[count_file])):
                    N_tot_mean=N_tot_mean+numpy.array(Ntot[p]);
                N_tot_mean = (N_tot_mean)/float(len(self.N_tot[count_file]));
                store_all_means.append([]);
                store_all_means[count_file] = N_tot_mean/binding_size[count_file] ;
                axes.plot(time_points,N_tot_mean/binding_size[count_file],color=colors[count_file],lw=2);
            count_file = count_file +1;

        if ~ind_plot:
            if ~no_legend:
                axes.legend(bbox_to_anchor=(1.05,1.1));
            axes.set_xlabel('Time');
            axes.set_ylabel('Number of chains');
        return(N_tot_mean,time_points,axes,store_all_means)



    """
    Function below is deprecated, and not used.
    """
#    def calculate_entropy_enthalpy_chains(self,scaled=True):
#        """
#            Function within the cluster_composition module
#            that calculates the entropy of the system via an approximation
#            that entropy lost in cluster is that "per" chain!
#        """
#
#        dt = float(self.input_params[0]['dt']);
#        binding_size = [];
#        volume_fraction = numpy.zeros((len(self.param_range),3))
#        PE_all_traj ={};
#        entropy_all_traj ={};
#        entropy_mean = [];
#        scaled_cs_mean = [];
#        PE_store_mean = [];
#        scaled_cluster_size_all_traj = {};
#        count = 0;
#
#        for param in numpy.arange(len(self.param_range)):
#            mypath = self.unique_paths[param];
#
#
#            (PE_mean,time,header,PE_all) = generate_average_plot(mypath,'PE',1,return_ind_traj=True);
#            L_interest = int( numpy.floor(len(time)/1));
#
#            (vol_A,vol_B,vol_C) = get_volume_fraction(self.input_params[param]);
#            volume_fraction[param,0] = vol_A;
#            volume_fraction[param,1] = vol_B;
#            volume_fraction[param,2] = vol_C;
#            if ((self.input_params[count_file]['N_A']) and (scaled)):
#                if self.input_params[count_file]['seq_A'].count('A'):
#                    binding_size.append(self.input_params[count_file]['N_A']*(float((self.input_params[count_file]['seq_A'].count('A'))))*(self.input_params[count_file]['N_bs_AB']+self.input_params[count_file]['N_bs_AC']));
#                else:
#                    binding_size.append(1.0)
#            else:
#                binding_size.append(1.0);
#
#            volume_rest = float(pow(self.input_params[count]['L'],3)) * (1-numpy.sum(volume_fraction[count,:])) ;
#
#            N_tot = numpy.reshape(list(self.N_tot[count].values()),(len(self.N_tot[count]),len(self.N_tot[count][0])));
#
#            entropy = numpy.multiply(self.N_tot[count_file][count],numpy.log(numpy.divide(4/3*numpy.pi*numpy.power(rg_all[:,0:L_interest],3),volume_rest)));
#            entropy_mean.append(numpy.mean(entropy,0));
#            scaled_cs_mean.append(numpy.mean(numpy.divide(cl_all[:,0:L_interest],binding_size[count]),axis=0));
#            PE_store_mean.append(PE_mean)
#            PE_all_traj[count] = PE_all;
#            entropy_all_traj[count] = entropy;
#            scaled_cluster_size_all_traj[count] = numpy.divide(cl_all[:,0:L_interest],binding_size[count])
#            count =count+1;
#
#        PE_store_mean = numpy.reshape(PE_store_mean,(len(param_range),len(time)))
#
#        entropy_mean = numpy.reshape(entropy_mean,(len(param_range),len(time)))
#        scaled_cs_mean = numpy.reshape(scaled_cs_mean,(len(param_range),len(time)))
#
#        return(PE_all_traj,entropy_all_traj,scaled_cluster_size_all_traj,time*dt,entropy_mean,PE_store_mean,scaled_cs_mean)
#



## FRAP analysis
class FRAP_analysis():

    """
        Module that performs computational FRAP analysis on MD simulation data.

        Class is defined by taking in the unique paths, as well as a time window (Start
        frame and end frama) and a fluorophore to FRAP (label of monomer), and returns the mean
        normalized intensity profile over the corresponding time window.
    """

    def __init__(self):
        self.int_total = [];

    def calc_distance_2points(self,pos1,pos2,L):
        """
        Function calculates distance between two points, accounting for
        periodic boundary conditions.

        **Input**

        -   pos1,pos2 - Positions of two monomers
        -   L - box size

        **Output**

        -   numpy.sqrt(y)) = Distance between two monomers

        """
        y =0;
        for count in numpy.arange(len(pos1)):
            if abs(pos1[count]-pos2[count]) > float(L)/2:
                y = y + numpy.power(L -abs(pos1[count]-pos2[count]),2);
            else:
                y = y + numpy.power(pos1[count]-pos2[count],2);
        return (numpy.sqrt(y));

    def calculate_intensity_FRAP(self,position,FRAP_centroid,FRAP_radius,intensity,L):
        """
        Function calculates FRAP intensity in control volume

        **Input**

        -   position - Positions of all flurophore monomers
        -   FRAP_centroid - center of control FRAP volume
        -   FRAP_radius - Radius of bleach volume
        -   intensity - intensity matrix of flurophores (+1 if on, +0  bleached)
        -   L - box size


        **Output**

        -   intensity_volume = Total intensity in bleach volume

        """
        intensity_volume = 0;
        count = 0;
        for pos in position:
            distance_between_points =self.calc_distance_2points(FRAP_centroid,pos,L)
            if distance_between_points <= FRAP_radius:
                intensity_volume += intensity[count]

            count = count+1;
        return(intensity_volume)

    def return_center_of_DNA(self,snap):
        """
        Function calculates centroid of DNA

        **Input**

        -   snap - MD snapshot at a given time

        **Output**

        -   pos_center - centroid of DNA chain

        """
        pos_A = [x for x in numpy.arange(snap.particles.N) if snap.particles.typeid[x]==0];

        if not pos_A:
            pos_center = numpy.zeros((1,3))
        else:
            p = int(numpy.floor(len(pos_A)/2));
            pos_center = snap.particles.position[pos_A[p],:]

        return (pos_center)

    def calculate_FRAP_curve(self,t,N_start,N_end,fluoro_type,FRAP_radius,moving_frame):
        """
        Function reconstructs a simulated FRAP window over a particular
        time window for a single trajectory. Typically called from return_mean_FRAP_data

        **Input**

        -   t - Matrix of all snapshots from .gsd file
        -   N_start -       Index of time_point to start calculation (step-index)
        -   N_end -         Index of time_point to end calculation (step-index)
        -   fluoro_type -   Label of fluorophore
            'B' for TF, 'C' for coactivator
        -   FRAP_radius -   size of bleached volume
        -   moving_frame -  Track the center of the DNA molecule to
            calculate FRAP around cluster TRUE enables this FALSE only uses fixed control volume

            step-index takes values from 0 to N_time_points-1

        **Output**

        -   int_store - List of FRAP intensities in control volume over time window
        -   len(parti_id)= Total number of fluorophores of type

        """
        init_snap = t[N_start];
        L = init_snap.configuration.box[0];

        typeid_rel=(init_snap.particles.types.index(fluoro_type));

        parti_id = [x for x in numpy.arange(init_snap.particles.N) if init_snap.particles.typeid[x]==typeid_rel];
    #     print("Number of particles of type {} is {}".format(fluoro_type,len(parti_id)))
        intensity = numpy.ones((len(parti_id),1));
        pos_center = self.return_center_of_DNA(init_snap);
        int_store = [];
        FRAP_radius =5;

        for snap in t[N_start-5:N_start]:
            a= snap.particles.position[parti_id,:];
            if moving_frame:
                int_store.append(self.calculate_intensity_FRAP(a,self.return_center_of_DNA(snap),FRAP_radius,intensity,L));
            else:
                int_store.append(self.calculate_intensity_FRAP(a,pos_center,FRAP_radius,intensity,L));


        #    Assume that each particle inside FRAP volume loses its intensity.
        count = 0;
        init_pos = t[N_start].particles.position[parti_id,:];
        for pos in init_pos:

            if moving_frame:
                distance_between_points =self.calc_distance_2points(self.return_center_of_DNA(init_snap),pos,L)
            else:
                distance_between_points =self.calc_distance_2points(pos_center,pos,L)

            if distance_between_points <= FRAP_radius:
                intensity[count] = 0;

            count = count+1;

        for snap in t[N_start:N_end]:
            a= snap.particles.position[parti_id,:];
            if moving_frame:
                int_store.append(self.calculate_intensity_FRAP(a,self.return_center_of_DNA(snap),FRAP_radius,intensity,L));
            else:
                int_store.append(self.calculate_intensity_FRAP(a,pos_center,FRAP_radius,intensity,L));
        return(int_store,len(parti_id))

    def return_mean_FRAP_data(self,mypath,N_start,N_end,fluoro_type,FRAP_radius,moving_frame):
        """
        Function reconstructs a simulated FRAP window over a particular
        time window over multiple trajectories.

        **Input**

        -   mypath -        Path to /.../Data/Trajectory/ of output data
        -   N_start -       Index of time_point to start calculation (step-index)
        -   N_end -         Index of time_point to end calculation (step-index)
        -   fluoro_type -   Label of fluorophore
            'B' for TF, 'C' for coactivator
        -   FRAP_radius -   size of bleached volume
        -   moving_frame -  Track the center of the DNA molecule to
            calculate FRAP around cluster
            TRUE enables this
            FALSE only uses fixed control volume

            step-index takes values from 0 to N_time_points-1

        **Output**

        -   int_store_mean -    List of FRAP intensities in control volume over time window
            averaged over multiple trajectories

        -   N_mol_mean-         Total number of fluorophores of type

        """

        int_store_mean = [];
        N_mol_mean = [];
        for path in mypath:
            onlyfiles = [f for f in os.listdir(path) if isfile(join(path, f))]
            Ntraj = len(onlyfiles);
            N_mol = [];
            intensity = [];
            for file in onlyfiles:
                gsd_file =join(path,file);
                f = gsd.fl.GSDFile(gsd_file, 'rb')
                t = gsd.hoomd.HOOMDTrajectory(f)

                intensity_temp,N_part = self.calculate_FRAP_curve(t,N_start,N_end,fluoro_type,FRAP_radius,moving_frame)
                intensity.append(intensity_temp);
                N_mol.append(N_part)

            int_total = numpy.reshape(intensity,(Ntraj,N_end-N_start+5))
            int_store_mean.append(numpy.mean(int_total,axis=0));
            N_mol_mean.append(numpy.mean(N_mol))
        int_store_mean = numpy.reshape(int_store_mean,(len(mypath),N_end-N_start+5));
        return(int_store_mean,N_mol_mean);

    def plot_FRAP_data(self,int_store_mean,N_mol_mean,param_range,param_name,input_parameters,mol_type,normalized=False,bleach_norm=False):

        """
        Function that plots FRAP data.

        FRAP experiments over single conditions are obtained and a list of
        mean intensities are constructed to pass into the system

        **Input**

        -   int_store_mean -Matrix of <FRAP> intensities of size N_params * N_FRAP_time_window

        -   N_mol_mean -    Number of fluorophores (used to make correction for FRAP after bleaching)

        -   N_end -         Index of time_point to end calculation (step-index)

        -   param_name-     The string name of the variable that is swept over.

        -   param_range -   Sorted list of the different values of parameters

        -   input_parameters - List of input_params for each condition

        -   mol_type -      Label of fluorophore
            'B' for TF, 'C' for coactivator

        -   normalized -      Normalize intensity by pre-bleach <intensity> in volume
            if set to TRUE. Default is FALSE.

        -   bleach_norm -   Normalize intensity for finite-size loss in
            bleached molecules. Correction takes the form

            I_{cor} = I * (1 - N_{bleached}/N_{total})

            step-index takes values from 0 to N_time_points-1

        **Output**

        -   int_store_mean -    List of FRAP intensities in control volume over time window
            averaged over multiple trajectories

        -   t -                  Time-steps

        """


        colors = iter(cm.afmhot(numpy.linspace(0, 0.5, len(param_range))))
        legend_labels = [];
        fig, axes = plt.subplots(1,1,sharex=True)
        make_nice_axis(axes)
        fig.suptitle('FRAP analysis');
        N_time_steps = len(int_store_mean[0,:]);
        size_step = input_parameters[0]['dt']*input_parameters[0]['log_step'];
        t = numpy.arange(N_time_steps)*size_step;
        for count in numpy.arange(len(param_range)):
            if normalized:
                if bleach_norm:
                    n_bleached = numpy.mean(int_store_mean[count,0:5]);
                    n_total = N_mol_mean[count];
                    int_store_mean[count,:]= int_store_mean[count,:]/numpy.mean(int_store_mean[count,0:5]);
                    int_store_mean[count,5:]= (int_store_mean[count,5:])/(1-n_bleached/n_total);

                else:
                    int_store_mean[count,:]= int_store_mean[count,:]/numpy.mean(int_store_mean[count,0:5]);
                axes.plot(t,int_store_mean[count,:],color=next(colors),lw=4);
            else:
                axes.plot(t,int_store_mean[count,:],color=next(colors),lw=4);
            if param_name.find('seq') == -1:
                legend_labels.append( param_name + ' = '+ str(round(param_range[count],2)));
            else:
                legend_labels.append( param_name + ' $ _{l} $ = '+ str(len(param_range[count])));

        axes.set_ylabel(str(mol_type));
        axes.set_xlabel('Time')

        axes.legend(legend_labels,bbox_to_anchor=(1.05, 1),fontsize =20)

        return (int_store_mean,t)








"""
    Making movies from simulation data!
"""
device = fresnel.Device(mode='cpu');

blue = fresnel.color.linear([0.25,0.5,1])*0.9;
orange = fresnel.color.linear([1.0,0.714,0.169])*0.9



def render_sphere_frame(frame,path_tracer, height=None):

    """
        render_sphere_frame takes in a HOOMD system snapshot
        and generates scene of the box, which it returns
        back through a path_tracer sample

        This function is an extension of the function provided
        in the HOOMD tool kit

        **Output** Returns frame render

    """


    if height is None:
        if hasattr(frame, 'configuration'):
            Ly = frame.configuration.box[1]
            height = Ly * numpy.sqrt(3)
        else:
            Ly = frame.box.Ly;
            height = Ly * numpy.sqrt(3)

    scene = fresnel.Scene(device)
    scene.lights = fresnel.light.cloudy();
    g = fresnel.geometry.Sphere(scene, position=frame.particles.position, radius=frame.particles.diameter*0.5)
    g.material = fresnel.material.Material(solid=1, color=blue, primitive_color_mix=1.0, specular=1.0, roughness=0.2)
    g.outline_width = 0.07
    scene.camera = fresnel.camera.orthographic(position=(height, height, height), look_at=(0,0,0), up=(0.2,1,0), height=height)

    g.color[frame.particles.typeid == 0] =[0.9,0,0];
    g.color[frame.particles.typeid == 1] = [228/256,190/256,108/256];
    g.color[frame.particles.typeid == 2] = [89/256,96/256,174/256];
    g.color[frame.particles.typeid == 3] = [0.5,0.5,0.5];
    g.color[frame.particles.typeid == 4] = [0.5,0.5,0.5];

    scene.background_color = (1,1,1)

    return path_tracer.sample(scene, samples=64, light_samples=20)


## Function that takes input to trajectory file, and outputfile name
## and saves a high quality gif
def save_movie( gsd_file,output_file,resolution,file_save=False,down_sample=1):


    """
        save_movie takes in a trajectory file (gsd_file) and a
        frame method which is typically render_sphere_frame (above)

        It needs the resolution list (2 values), output_file name,
        and down_sampling frequency to get the right frames

        The resolution list sets the quality of the image. Typical values
        include 300*300 dpi.

        Higher values of down_sample constitutes movies with lesser frames,
        but equally sampled from the trajectory data.

        if file_save is set to True --> the movies are outputted to the output_file

        **Output of function**

        -   Returns last snap that was framed.

    """
    path_tracer = fresnel.tracer.Path(device,resolution[0],resolution[1])

    f = gsd.fl.GSDFile(gsd_file, 'rb')
    t = gsd.hoomd.HOOMDTrajectory(f)

    a = render_sphere_frame(frame=t[0],path_tracer=path_tracer);

    if tuple(map(int, (PIL.__version__.split(".")))) < (3,4,0):
        print("Warning! Movie display output requires pillow 3.4.0 or newer.")
        print("Older versions of pillow may only display the first frame.")

    im0 = PIL.Image.fromarray(a[:,:, 0:3], mode='RGB').convert("P", palette=PIL.Image.ADAPTIVE);
    ims = [];
    points = numpy.linspace(1,len(t)-1,(len(t)-1)/down_sample);
    print(points)
    for point in points:
        f = t[int(numpy.floor(point))];
        a = render_sphere_frame(frame=f,path_tracer=path_tracer);
        im = PIL.Image.fromarray(a[:,:, 0:3], mode='RGB')
        im_p = im.quantize(palette=im0);
        ims.append(im_p)
    if file_save:
        if not os.path.exists(os.path.dirname(output_file)):
            os.makedirs(os.path.dirname(output_file),exist_ok=True);

        im0.save(output_file, 'gif', save_all=True, append_images=ims, duration=1500, loop=0)

    return (f)

class contact_frequency_matrix():

    """
        Module that helps estimate the contact frequency matrix of the DNA
        chain for each param_index, and this average is returned.

    """

    def __init__(self):
        """
        Intializes an empty list to variables of this class.

        """
        self.contact_freq = [];

    def return_pos_of_DNA(self,snap,input_params):
        """
        Function calculates positions of DNA monomers.

        **Input**

        -   snap - MD snapshot at a given time
        -   input_params - Input parameters

        **Output**

        -   pos_A - list of ID's of DNA monomers

        """
        monomer_to_chain_map = input_params['MC_map'];
        chain_type = input_params['CT_map'];
        N_monomers = len(monomer_to_chain_map.keys())
        pos_A = [x for x in numpy.arange(N_monomers) if chain_type[monomer_to_chain_map[x]]=='A'];

        return (pos_A);

    def calc_distance_2points(self,pos1,pos2,L):

        """
        Function calculates distance between two points, accounting for
        periodic boundary conditions.

        **Input**

        -   pos1,pos2 - Positions of two monomers
        -   L - box size

        **Output**

        -   numpy.sqrt(y)) = Distance between two monomers

        """

        y =0;
#        print((pos1),pos2)
        for count in numpy.arange(len(pos1)):
            if abs(pos1[count]-pos2[count]) > float(L)/2:
                y = y + numpy.power(L -abs(pos1[count]-pos2[count]),2);
            else:
                y = y + numpy.power(pos1[count]-pos2[count],2);

        return (numpy.sqrt(y));

    def estimate_snap_contact_frequency(self,snap,in_params):

        """
        Function calculates the contact frequency for an individual
        snapshot from the MD simulation. The input params are passed from
        the overall function.

        Typically called and employed from calculate_contact_frequency

        **Returns**

        -   c_map = Contact frequency matrix of that snapshot

        """
        pos = self.return_pos_of_DNA(snap,in_params);
        L = snap.configuration.box[0]
#        print(L)
        n_rows = len(pos);
        c_map = numpy.zeros((n_rows,n_rows));
        for i in numpy.arange(0,n_rows-1):
            for j in numpy.arange(1,n_rows):
                c_map[i,j] = self.calc_distance_2points(snap.particles.position[pos[i],:],snap.particles.position[pos[j],:],L);
                c_map[j,i] = c_map[i,j];

        return c_map;

    def calculate_contact_frequency(self,param_name,param_range,unique_paths,input_params,N_start,N_end,threshold=2.5):

        """
        Module that helps estimate the contact frequency matrix of the DNA
        chain for each param_index, and this average is returned.

        For pairs of monomers < threshold distance, their contact is counted
        as +1, or else 0. This is used to estimate the contact frequency
        matrix.

        **Input to function**

        -   param_name =   The string name of the variable that is swept over.

        -   param_range =  sorted list of the different values of parameters

        -   unique_paths = list of paths to each observable folder for different
            conditions in the output data-set
            (e.g. list of folder/condition/Data/Observables
            for each condition)

        -   input_params = list of input_params for each condition

        -   N_start = Index of time_point to start calculation

        -   N_end =    Index of time_point to end calculation

        -   threshold = radius threshold to estimate the contact frequency
            (threshold sets when the pairs of monomers are
            cross-linked)


        **Output variables**

        -   *self*.contact_freq_all - Averaged contact frequency matrix (after
            cross-linking)

        -   *self*.contact_dist_all - Averaged contact distance matrix (after
            cross-linking)

        """
        unique_paths = [uni_path.replace('Observables','Trajectory') for uni_path in unique_paths]
        self.contact_freq_all = [];
        self.contact_dist_all = [];
        count =0;

        for path in unique_paths:
            in_params = input_params[count];
            onlyfiles = [f for f in os.listdir(path) if isfile(join(path, f))]
            traj =0;

            contact_dist_traj = [];
            contact_frequency_traj = [];

            for file in onlyfiles:
                gsd_file =join(path,file);
                f = gsd.fl.GSDFile(gsd_file, 'rb')
                t = gsd.hoomd.HOOMDTrajectory(f)

                contact_distance = [];
                contact_frequency = [];

                time_point =0;
                for snap in t[N_start:N_end+1]:
#                    print(snap.particles.N)
                    contact_distance_temp = self.estimate_snap_contact_frequency(snap,in_params)
#                    print(contact_frequency_temp)
                    if not time_point:

                        contact_distance = contact_distance_temp;
                        contact_frequency = (contact_distance_temp<threshold).astype(int)
                        time_point +=1;
                    else:
                        contact_distance = contact_distance_temp + contact_distance;
                        contact_frequency = contact_frequency + (contact_distance_temp<threshold).astype(int);
                        time_point +=1;

                contact_frequency_traj.append([])
                contact_dist_traj.append([]);

                contact_dist_traj[traj] = contact_distance/(N_end-N_start+1);
                contact_frequency_traj[traj] = contact_frequency.astype(float)/(N_end-N_start+1);
#                print(contact_dist_traj[traj],contact_frequency_traj[traj])
                traj = traj+1;

            self.contact_freq_all.append([])
            self.contact_dist_all.append([])
            self.contact_dist_all[count] = sum(contact_dist_traj)/traj;
            self.contact_freq_all[count] = sum(contact_frequency_traj)/traj;

            count = count+1;

        return(self.contact_freq_all,self.contact_dist_all);
