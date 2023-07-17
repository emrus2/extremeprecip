# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 13:49:49 2023

@author: emrus2

Uses:
    I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_{numpatterns}NodeAssignMonthlyHistogram.npy
    
to plot the monthly node frequency (for cummulative nodes)

Instead of the absolute number of days, instead counts the proportion of days assigned to node in each month

UPDATED 6/12/2023
"""
#%% IMPORT MODULES
#import netCDF4 as nc #if this doesn't work, try to reinstall in anaconda prompt using
    #conda install netcdf4
import numpy as np
import os
import matplotlib.pyplot as plt
#from matplotlib.colors import ListedColormap
#import matplotlib.cm as cm
#import matplotlib.transforms as mtransforms
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
#from mpl_toolkits.basemap import Basemap #installed using 
    #conda install -c anaconda basemap
import scipy.io

#%% IMPORT SOM DATA

# define metvar
#metvars = ['Z500', 'SLP', '850T', '300W', 'IVT','Z500Anom','300W']
metvar = 'IVT'
numpatterns = 9
percentile = 90

# change directory and import SOM data from .mat file
mat_dir='I:\\Emma\\FIROWatersheds\\Data\\SOMs\\SomOutput'
os.chdir(mat_dir)
soms = scipy.io.loadmat(f'{metvar}_{percentile}_{numpatterns}_sompatterns.mat')
pat_prop = np.squeeze(soms['pat_prop'])
pat_freq = np.squeeze(soms['pat_freq'])
#%% IMPORT HISTOGRAM DATA
data = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_{numpatterns}NodeAssignMonthlyHistogram.npy')
#data[data == 0] = 'nan'

#%% CLUSTERED BAR CHART - 12 NODE
if numpatterns == 12:
    months = ("October", "November", "December", "January", "February", "March")
    colors = ('orangered','indianred','gold','lightgreen','mediumseagreen','olive','cornflowerblue','royalblue','plum','darkorchid','magenta','slategrey')
    
    if percentile == 95:
        ymax = 10
        yint = 2
    elif percentile == 90:
        ymax = 13
        yint = 3
        
    x = np.arange(len(months))  # the label locations
    width = 0.05  # the width of the bars
    multiplier = 0
    
    fig, ax = plt.subplots(layout='constrained')
        
    for i, node in enumerate(data):
        #print(i,node,pat_freq[i])
        node_rel = (node / pat_freq[i]) * 100
        offset = width * multiplier
        rects = ax.bar(x + offset, node_rel, width, label=i+1,align='center',color=colors[i])
        #ax.bar_label(rects, padding=3)
        multiplier += 1
    
    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Days',fontweight='bold')
    ax.set_xlabel('Month',fontweight='bold')
    #ax.set_title('Node Frequency',fontweight='bold')
    ax.set_xticks(x + width + 0.2, months)
    ax.legend(loc='upper left', ncols=2)
    ax.set_ylim(0, ymax)
    
    save_dir='I:\\Emma\\FIROWatersheds\\Figures\\NodeHistograms'
    os.chdir(save_dir)
    plt.savefig(f'{percentile}_{numpatterns}_NodeFreqMonthly.png',dpi=300,bbox_inches='tight',pad_inches=0.08)
    plt.show()

#%% CLUSTERED BAR CHART - 9 NODE
if numpatterns == 9:
    months = ("October", "November", "December", "January", "February", "March")
    colors = ('tomato','indianred','gold','lightgreen','mediumseagreen','cornflowerblue','royalblue','plum','darkorchid','magenta',)
    
    # if percentile == 95:
    #     ymax = 10
    #     yint = 2
    # elif percentile == 90:
    #     ymax = 18
    #     yint = 3
    
    x = np.arange(len(months))  # the label locations
    width = 0.09  # the width of the bars
    multiplier = 0
    
    fig, ax = plt.subplots(layout='constrained')
        
    for i, node in enumerate(data):
        node_rel = (node / pat_freq[i]) * 100
        print(node_rel)
        offset = width * multiplier
        rects = ax.bar(x + offset, node_rel, width, label=i+1,align='center',color=colors[i])
        #ax.bar_label(rects, padding=3)
        multiplier += 1
    
    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Days (%)',fontweight='bold')
    ax.set_xlabel('Month',fontweight='bold')
    #ax.set_title('Node Frequency',fontweight='bold')
    ax.set_xticks(x + width + 0.26, months)
    #ax.legend(bbox_to_anchor=(1.0, 1.0),ncols=1)
    ax.legend(loc='upper center',ncols=9,columnspacing=1.4,handletextpad=0.4,fontsize=9.7)
    ax.set_ylim(0, 46.5)
    ax.set_xlim(-0.05,5.8)
    ax.set_yticks(np.arange(0,46.5,5))
    
    save_dir='I:\\Emma\\FIROWatersheds\\Figures\\NodeHistograms'
    os.chdir(save_dir)
    plt.savefig(f'{percentile}_{numpatterns}_NodeProportionMonthly.png',dpi=300,bbox_inches='tight',pad_inches=0.08)
    plt.show()