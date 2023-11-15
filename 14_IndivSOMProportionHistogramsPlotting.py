# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 13:49:49 2023

@author: emrus2

Uses:
    I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_{numpatterns}NodeAssignMonthlyHistogram.npy
    
to plot the monthly node frequency (for individual nodes)

Instead of the absolute number of days, instead counts the proportion of days assigned to node in each month


For a 9-node SOM

UPDATED 6/12/2023
"""
#%% IMPORT MODULES
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
import scipy.io

#%% IMPORT SOM DATA

# define metvar
metvar = 'IVT'
numpatterns = 12
percentile = 90
clusters = 5
lettered = True

# change directory and import SOM data from .mat file
mat_dir='I:\\Emma\\FIROWatersheds\\Data\\SOMs\\SomOutput'
os.chdir(mat_dir)
soms = scipy.io.loadmat(f'{metvar}_{percentile}_{numpatterns}sompatterns_{clusters}d.mat')
pat_prop = np.squeeze(soms['pat_prop'])
pat_freq = np.squeeze(soms['pat_freq'])
#%% IMPORT HISTOGRAM DATA
data = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_{numpatterns}NodeAssignMonthlyHistogram_{clusters}d.npy')

#%% SUBPLOTS OF NODES
if numpatterns == 9:
    months = ("O", "N", "D", "J", "F", "M")
    #color = 'tomato'
    colors = ('tomato','indianred','gold','lightgreen','mediumseagreen','cornflowerblue','royalblue','plum','darkorchid','magenta',)
    
    
    x = np.arange(len(months))  # the label locations
    width = 0.8  # the width of the bars
    multiplier = 0
    
    fig = plt.figure(figsize=(7.2,5))
    
    if clusters == 1:
        ymax = 50
        yint = 10
    elif clusters == 5:
        ymax = 70
        yint = 15
    
    
    for i,node in enumerate(data):
        perc_assigned = str(round(pat_prop[i]*100,1))
        # precipavg = precip_rounded[i]
        node_rel = (node / pat_freq[i]) * 100
        round_node_rel = [ round(elem, 0) for elem in node_rel]
        offset = width * multiplier
        gbcolor = 1.1-(pat_freq[i]/44)
        patfreq_col = (1,gbcolor,gbcolor)
        ax = fig.add_subplot(3,3,i+1)
        sublabel_loc = mtransforms.ScaledTranslation(4/72, -4/72, fig.dpi_scale_trans)
        ax.text(x=0.0, y=1.0, s=i+1, transform=ax.transAxes + sublabel_loc,
            fontsize=10, fontweight='bold', verticalalignment='top', 
            bbox=dict(facecolor='1', edgecolor='none', pad=1.5),zorder=3)
        freqs = plt.bar(x + offset,height=round_node_rel,width=width,color=colors[i],label=i+1,align='center')
        ax.bar_label(freqs, padding=3)
        ax.set_ylim(0, ymax)
        ax.set_xlim(-0.5,5.5)
        ylabels = np.arange(0,ymax,yint)
        if i == 0 or i == 3:
            ax.set_yticks(ylabels)
            ax.set_xticks(x, [])
            if i == 3:
                ax.set_ylabel('Days (%)',fontweight='bold')
        elif i == 7 or i == 8:
            ax.set_yticks(ylabels,[])
            ax.set_xticks(x, months)
            if i == 7:
                ax.set_xlabel('Month',fontweight='bold')
        elif i == 6:
            ax.set_yticks(ylabels)
            ax.set_xticks(x, months)
        else:
            ax.set_yticks(ylabels,[])
            ax.set_xticks(x, [])
        ax.tick_params(direction='in',which='both',axis='y')
    
    #CUSTOMIZE SUBPLOT SPACING
    fig.subplots_adjust(left=0.065,right=0.985,bottom=0.082, top=0.985,hspace=0.05, wspace=0.05) #bottom colorbar
    
    
    save_dir='I:\\Emma\\FIROWatersheds\\Figures\\NodeHistograms'
    os.chdir(save_dir)
    plt.savefig(f'{percentile}_{numpatterns}_NodeProportionbyNode_{clusters}d.png',dpi=300)
    plt.show()
    #%%
if numpatterns == 12:

    months = ("O", "N", "D", "J", "F", "M")
    #color = 'tomato'
    colors = ('tomato','cornflowerblue','lightgreen','darkorchid','gold','lightblue','plum','mediumseagreen','indianred','royalblue','grey',(.9,0,.9))
    
    
    x = np.arange(len(months))  # the label locations
    width = 0.8  # the width of the bars
    multiplier = 0
    
    if lettered == True:
        fig = plt.figure(figsize=(7.2,6.9))
        fig.suptitle('b)',fontsize=11,fontweight="bold",y=0.997,x=0.07)
    else:
        fig = plt.figure(figsize=(7.2,6.7))
    
    
    if clusters == 1:
        ymax = 50
        yint = 10
    elif clusters == 5:
        ymax = 71
        yint = 15
    
    
    for i,node in enumerate(data):
        perc_assigned = str(round(pat_prop[i]*100,1))
        # precipavg = precip_rounded[i]
        node_rel = (node / pat_freq[i]) * 100
        round_node_rel = [ round(elem, 0) for elem in node_rel]
        offset = width * multiplier
        gbcolor = 1.1-(pat_freq[i]/44)
        patfreq_col = (1,gbcolor,gbcolor)
        ax = fig.add_subplot(4,3,i+1)
        sublabel_loc = mtransforms.ScaledTranslation(4/72, -4/72, fig.dpi_scale_trans)
        ax.text(x=0.0, y=1.0, s=i+1, transform=ax.transAxes + sublabel_loc,
            fontsize=10, fontweight='bold', verticalalignment='top', 
            bbox=dict(facecolor='1', edgecolor='none', pad=1.5),zorder=3)
        freqs = plt.bar(x + offset,height=round_node_rel,width=width,color=colors[i],label=i+1,align='center')
        ax.bar_label(freqs, padding=3)
        ax.set_ylim(0, ymax)
        ax.set_xlim(-0.5,5.5)
        ylabels = np.arange(0,ymax,yint)
        if i in range(0,7,3):
            ax.set_yticks(ylabels)
            ax.set_xticks(x, [])
            if i == 3:
                ax.set_ylabel('Days (%)',fontweight='bold')
                ax.yaxis.set_label_coords(-0.12,-0.1)
        elif i == numpatterns-2 or i == numpatterns-1:
            ax.set_yticks(ylabels,[])
            ax.set_xticks(x, months)
            if i == 10:
                ax.set_xlabel('Month',fontweight='bold')
        elif i == numpatterns-3:
            ax.set_yticks(ylabels)
            ax.set_xticks(x, months)
        else:
            ax.set_yticks(ylabels,[])
            ax.set_xticks(x, [])
        ax.tick_params(direction='in',which='both',axis='y')
    
    #CUSTOMIZE SUBPLOT SPACING
    if lettered == True:
        fig.subplots_adjust(left=0.065,right=0.985,bottom=0.065, top=0.975,hspace=0.05, wspace=0.05) #bottom colorbar
    else:
        fig.subplots_adjust(left=0.065,right=0.985,bottom=0.065, top=0.985,hspace=0.05, wspace=0.05) #bottom colorbar
    
    
    save_dir='I:\\Emma\\FIROWatersheds\\Figures\\NodeHistograms'
    os.chdir(save_dir)
    if lettered == True:
        plt.savefig(f'{percentile}_{numpatterns}_NodeProportionbyNode_{clusters}d_lettered.png',dpi=300)
    else:
        plt.savefig(f'{percentile}_{numpatterns}_NodeProportionbyNode_{clusters}d.png',dpi=300)
    plt.show()