#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 11:23:20 2018
This file is primarily for functions that aid post-processing of droplet characterization analysis
from MATLAB files
@author: krishna
"""

import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile
import numpy
import matplotlib.pyplot as plt
import os
import matplotlib.cm as cm
import matplotlib; mpl =matplotlib
from scipy import stats;


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

#colors = ['#ccebc5','#b3cde3','#fbb4ae','#decbe4','#fed9a6'];

def make_nice_axis(ax):
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



def parse_data(folder_base,conditions,sheet_name,clean_xdata = False,pos_clean=-1):
    
    """
        The base path to directory containing Output directories is passed as
        folder_base. Different conditions represent distinct lines which are usually
        also distinct output_condition folders. sheet_name has meta_data on
        image_analysis parameters!
    """
    folder_base = folder_base + 'Output_'  
    ydata = {};
    xdata = {};
    for condition in conditions:
        p = os.listdir(folder_base+condition+'/');
        files = [x for x in p if x.find('xlsx') >-1]
        file_name = folder_base +condition + '/' + files[0]    
        df = pd.read_excel(file_name, sheetname=sheet_name);
        headers = (df.columns);
        df=df.dropna();
#        print(df,condition.split('-')[0])
        ydata[condition] = numpy.array(df.values[:,0:-1],dtype='float');
        if not clean_xdata:
            xdata[condition] = list(df[headers[-1]].values);
            xdata['data_type'] = type(xdata[condition][0])
        else:
            xdata[condition] = [float((x.split('_')[pos_clean]).replace('-','')) for x in df[headers[-1]].values]
            xdata['data_type'] = type(xdata[condition][0]);
         
        xdata[condition+'_header'] = headers;
    return (xdata,ydata);



def parse_data_individual(folder_base,conditions,sheet_name,sub_folder,clean_xdata = False):
    
    """
        The base path to directory containing Output directories is passed as
        folder_base. Different conditions represent distinct lines which are usually
        also distinct output_condition folders. sheet_name has meta_data on
        image_analysis parameters!
    """
    folder_base = folder_base + 'Output_'  
    ydata = {};
    xdata = {};
    for condition in conditions:
        p = os.listdir(folder_base+condition+'/' + sub_folder + '/');
        files = [x for x in p if x.find('xlsx') >-1];
        for file in files:
            file_name = folder_base +condition + '/' + sub_folder + '/' + file;    
            df = pd.read_excel(file_name, sheetname=sheet_name);
            headers = (df.columns);
            df=df.dropna();
            #        print(df,condition.split('-')[0])
            file = file.strip('.xlsx');
            
            if not clean_xdata:
                xdata[condition + '_' + file] = list(df[headers[-1]].values);
                ydata[condition+'_'+file] = numpy.array(df.values[:,0:-1],dtype='float');
            else:    
                value = float(file.split('MED1_')[1].split('_')[0]);
                xdata[value] = value; 
                ydata[value] = numpy.array(df.values[:,0:-1],dtype='float');
    return (xdata,ydata);




def plot_all_data(xdata,ydata,error_type='shaded',file_save = False,figure_name='test',colors = None, labels = None,data_pos=0,line_plot=True):

    """
        xdata and ydata are typically dictionaries with the same key linking a 
        particulary xdata set with a y data-set. Colors are passed for each "unique"
        condition i.e. N of conditions, and so are labels.
    """
    N_conditions = len(ydata.keys());
    
    fig,axes = plt.subplots(1,1)
    make_nice_axis(axes)

    if colors is None:
        colors = list(cm.rainbow(numpy.linspace(0,1,N_conditions)));

    if labels is None:
        labels = ['data_'+ str(x+1) for x in numpy.arange(N_conditions)];

    y_count = 0;
    for condition in ydata.keys():
        x_d = xdata[condition];
        y_d = ydata[condition];

        if error_type is 'bars':
            if line_plot:
                axes.errorbar(x_d,y_d[:,data_pos],y_d[:,data_pos+1],label=labels[y_count],lw =3,color = colors[y_count],elinewidth=1,capsize=5);
            else:
                axes.errorbar(x_d,y_d[:,data_pos],y_d[:,data_pos+1],label=labels[y_count],fmt='o',color = colors[y_count],elinewidth=1,capsize=5);
        elif error_type is 'shaded':
            axes.plot(x_d,y_d[:,data_pos],lw=3,color=colors[y_count],label=labels[y_count]);
            axes.fill_between(x_d,y_d[:,data_pos]-y_d[:,data_pos+1],y_d[:,data_pos]+y_d[:,data_pos+1],facecolor=colors[y_count],alpha=0.1)

        y_count = y_count+1;

    axes.legend(bbox_to_anchor=(1.05,1));

    if file_save:
        if not os.path.exists(os.path.dirname(figure_name)):
            os.makedirs(os.path.dirname(figure_name));

        plt.savefig(fname = figure_name + '.svg',format='svg',dpi=600,bbox_inches='tight')
        plt.savefig(fname = figure_name + '.png',format='png',dpi=600,bbox_inches='tight')

    return (fig,axes);


def plot_across_channels(xdata,ydata,colors = None, labels = None,error_type='shaded',line_plot=True,data_init_pos=2):

    """
        Plots a given quantity (like condensed fraction) across all conditions
        and all relevant imaged channels!
    """
    N_conditions = len(ydata.keys());
    conditions = list(ydata.keys());
    N_channels = int(ydata[conditions[0]].shape[1]/6);
    if colors is None:
        colors = list(cm.rainbow(numpy.linspace(0,1,N_conditions)));

    if labels is None:
        labels = ['data_'+ str(x+1) for x in numpy.arange(N_conditions)];
    
    axes = [];
    figures=  [];
    channels = [488,561,640];
    for channel in numpy.arange(N_channels):    
        (temp_fig,temp_axes) = plot_all_data(xdata,ydata,error_type=error_type,file_save = False,figure_name='test',colors = colors, labels = labels,data_pos=channel*6+data_init_pos,line_plot=line_plot)    
        temp_axes.legend(bbox_to_anchor=(1.05,1));
        temp_axes.set_title('Channel = ' + str(channels[channel]) )
        axes.append([]);
        axes[channel]= temp_axes;
        figures.append([]);
        figures[channel] = temp_fig;
    return (figures,axes);
    
def production_figures(xdata,ydata,colors = None, labels = None,error_type='shaded',line_plot=True,output_folder=None,xscale='linear',yscale='linear',xlabel='Concentration'):
    
    """
        Plots all relevant quantities (like condensed fraction) across all conditions
        and all relevant imaged channels!
    """

    if output_folder is not None:
        os.makedirs(output_folder,exist_ok=True)
        os.makedirs(output_folder,exist_ok=True)
        os.makedirs(output_folder,exist_ok=True)

    channels = ['488','561','640'];

    conditions = list(ydata.keys());
 
    if labels is None:
        labels = conditions;
    axes = [];
    figures = [];
    temp_figures,temp_axes = (plot_across_channels(xdata,ydata,colors = colors,labels = labels,error_type=error_type,line_plot=line_plot,data_init_pos=0));
    for ax in temp_axes:
        ax.set_ylabel('Partition ratio');
        if xdata['data_type'] is not str:
            ax.set_yscale(yscale);
            ax.set_xscale(xscale);
        ax.set_xlabel(xlabel);

    if output_folder is not None:
        count = 0;
        for fig in temp_figures:
            fig.savefig(output_folder+'/partition_' + channels[count]+ '.svg',dpi=600,bbox_inches='tight')
            count = count+1;
    
    axes.append(temp_axes);   
    figures.append(temp_figures)
    temp_figures,temp_axes = (plot_across_channels(xdata,ydata,colors = colors,labels = labels,error_type=error_type,line_plot=line_plot,data_init_pos=2));
    for ax in temp_axes:
        ax.set_ylabel('Condensed fraction');
        ax.set_xlabel(xlabel);
        if xdata['data_type'] is not str:
            ax.set_yscale(yscale);
            ax.set_xscale(xscale);

        
    if output_folder is not None:
        count = 0;
        for fig in temp_figures:
            fig.savefig(output_folder+'/cf_' + channels[count]+ '.svg',dpi=600,bbox_inches='tight')
            count = count+1;

    
    axes.append(temp_axes); 
    figures.append(temp_figures)

    temp_figures,temp_axes =(plot_across_channels(xdata,ydata,colors = colors,labels = labels,error_type=error_type,line_plot=line_plot,data_init_pos=4));
    for ax in temp_axes:
        ax.set_ylabel('Total intensity');
        if xdata['data_type'] is not str:
            ax.set_yscale(yscale);
            ax.set_xscale(xscale);

        ax.set_xlabel(xlabel);
       
    if output_folder is not None:
        count = 0;
        for fig in temp_figures:
            fig.savefig(output_folder+'/totali_' + channels[count]+ '.svg',dpi=600,bbox_inches='tight')
            count = count+1;

    axes.append(temp_axes); 
    figures.append(temp_figures)

    return(figures,axes);

def infer_saturation_concentration(xdata,ydata,threshold,pos):
    
    """
        Infers the saturation concentration as the value at which the condensed
        fraction in the "scaffold" channel crosses the threshold value that suggests.
        
        The precise calculation is carried out as follows: The concentrations
        before and after the 0.5 values are identified, and a linear fit is 
        obtained between the condensed fraction means and the log of the concentrations.
        Subsequently, the c.f. 0.5 value is inferred back to a specifc concentration
        value. The std deviation of this bound is obtained by fitting lines
        to the mean1+std1 --> mean2+std2 and mean1-std1 --->mean2-std2
    """

    conditions = list(ydata.keys());

    saturation_concentrations = {};
    saturation_concentrations_error = {};
    
    for cond in conditions:
        y = ydata[cond][:,pos]*100;
        yerr = ydata[cond][:,pos+1]*100;
        x = xdata[cond];
    
        
        idx = (numpy.abs(y - threshold)).argmin()
        if y[idx] > threshold:
            idx = idx-1;
#        print(idx)
        x1 = numpy.log10(x[idx]);
        x2 = numpy.log10(x[idx+1]);
        y1 = y[idx];
        y2 = y[idx+1];
        y1_std = yerr[idx];
        y2_std = yerr[idx+1];
        
        # inferring slopes from the error curves
        inferred_slope = (y2-y1)/(x2-x1);
        inferred_slope_top = (y2+y2_std-y1-y1_std)/(x2-x1);
        inferred_slope_bottom = (y2-y2_std-y1+y1_std)/(x2-x1);
        
        
        c_sat = numpy.power(10,(0.5-y1)/inferred_slope + x1);
        c_sat_std = 0.5*(-numpy.power(10,(0.5-y1-y1_std)/inferred_slope_top + x1) + numpy.power(10,(0.5-y1+y1_std)/inferred_slope_bottom + x1)) ;

    #     print(conditions[count],c_sat,c_sat_std)
    
        saturation_concentrations[cond] = (c_sat)
        saturation_concentrations_error[cond]=(c_sat_std)
        
    return(saturation_concentrations,saturation_concentrations_error)
    
def infer_individual_saturation_concentration(xdata,ydata,threshold,pos,all_pairs=False):
    a = sorted(list((xdata.keys())));
    conc = [];
    cf = [];

    before = [];
    after = [];
    mean_data = [];
    for key in sorted(a):
        conc.append(xdata[key]);
        cf.append(ydata[key][:,pos]*100)
        mean_data.append(float(ydata[key][-2,pos]*100));
        
    """closest value to threshold"""
#    idx = (numpy.abs(numpy.array(mean_data) - threshold)).argmin()

    """ first value greater than threshold"""
    idx = (numpy.argmax(numpy.array(mean_data) > threshold))

    if (mean_data[idx]) >threshold:
        idx = idx-1;
    
    before = numpy.array(cf[idx]);
    after = numpy.array(cf[idx+1]);
    data = numpy.concatenate((before,after));
    xcloud= numpy.concatenate((numpy.ones((len(before)))*conc[idx],numpy.ones((len(after)))*conc[idx+1])); 

    slope, intercept, r_value, p_value, std_err = stats.linregress((data),xcloud)
#    print(slope,intercept,r_value,p_value,std_err)
    if all_pairs:
        c_sat_all = [];
        for i in numpy.arange(len(before)):
            for j in numpy.arange(len(after)):
                slope_temp = (after[j]-before[i])/(conc[idx+1] - conc[idx]);
                c_sat_temp = conc[idx] + (threshold-before[i])/slope_temp;
                c_sat_all.append(c_sat_temp);
    c_sat = threshold*slope + intercept ;

    
    return (c_sat,std_err) if not all_pairs else (c_sat,std_err,c_sat_all)



def return_csat_threshold_sweep(folder_base,conditions,sheet_name,sub_folder,clean_xdata=True,threshold=[0,0.5,1],pos=6):
    
    """
        Calculate the inferred threshold csats at different threshold values
        and return those values.
    """
    csat_all = {};
    csat_all_errors = {};
    N_replicates = {};
    for cond in conditions:
        (xdata_temp,ydata_temp) = parse_data_individual(folder_base,[cond],sheet_name,sub_folder,clean_xdata=True)
        csat_all[cond] ={};  
        csat_all_errors[cond] = {};

        for t in threshold:

            (csat,std_error) =infer_individual_saturation_concentration(xdata_temp,ydata_temp,t,pos,all_pairs=False)
            csat_all[cond][t] = csat;
            csat_all_errors[cond][t]= std_error;
            keys = list(ydata_temp.keys())
            N_replicates[cond] = ((ydata_temp[keys[0]].shape[0]-2));

    return (csat_all,csat_all_errors)
    

def return_csat_sweep_plots(folder_base,conditions,sheet_name,sub_folder,clean_xdata=True,threshold=[0,0.5,1],pos=6):    

    """
        Return summary plots of threshold sweeps on csats
    """

    (csat_sweep,csat_error_sweep) = return_csat_threshold_sweep(folder_base,conditions,sheet_name,sub_folder,clean_xdata=True,threshold=threshold,pos=6);

    colors = ['#7fc97f','#beaed4','#fdc086','#386cb0','#ff7f00','#fb9a99','#a6cee3'];

    fig,axes = plt.subplots(2,1,sharex=True,figsize=(20,16))
    make_nice_axis(axes[0])
    count= 0;

    for cond in conditions:
        axes[0].errorbar(threshold, list(csat_sweep[cond].values()), yerr=list(csat_error_sweep[cond].values()), color=colors[count],label=cond);
        count = count+1;

    axes[0].legend(bbox_to_anchor=(1.05,0.5))

    axes[0].set_ylabel(' $ C_{sat} $ (nM)');
    make_nice_axis(axes[1]);
    axes[1].plot(threshold, numpy.divide(list(csat_sweep[conditions[0]].values()),list(csat_sweep[conditions[1]].values())),color='black');
    axes[1].set_ylabel('Ratio of $ C_{sat} $')
    axes[1].set_xlabel('c.f % threshold ')
    axes[1].set_xlim((0,1))

    return (fig,axes)