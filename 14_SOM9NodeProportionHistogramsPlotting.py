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
#import netCDF4 as nc #if this doesn't work, try to reinstall in anaconda prompt using
    #conda install netcdf4
import numpy as np
import os
import matplotlib.pyplot as plt
#from matplotlib.colors import ListedColormap
#import matplotlib.cm as cm
import matplotlib.transforms as mtransforms
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

#%% IMPORT AVERAGE PRECIP DATA
precip = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_{numpatterns}NodeAssignAvergagePrecip.npy')
precip_rounded = np.round(a=precip,decimals=1)

#%% SUBPLOTS OF NODES

months = ("O", "N", "D", "J", "F", "M")
#color = 'tomato'
colors = ('tomato','indianred','gold','lightgreen','mediumseagreen','cornflowerblue','royalblue','plum','darkorchid','magenta',)


x = np.arange(len(months))  # the label locations
width = 0.8  # the width of the bars
multiplier = 0

fig = plt.figure(figsize=(7.2,5))
#fig.suptitle('Monthly Node Frequency',fontsize=13,fontweight="bold",y=0.9875)

# if percentile == 95:
#     ymax = 11
#     yint = 2
# elif percentile == 90:
#     ymax = 20.6
#     yint = 3

ymax = 50
yint = 10


for i,node in enumerate(data):
    perc_assigned = str(round(pat_prop[i]*100,1))
    precipavg = precip_rounded[i]
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
    # ax.text(x=0.77, y=1.0, s=f'{perc_assigned}%', transform=ax.transAxes + sublabel_loc,
    #     fontsize=9, fontweight='bold', verticalalignment='top', color = 'red',
    #     bbox=dict(facecolor='1', edgecolor='none', pad=1.5),zorder=3)
    # ax.text(x=0.45, y=1.0, s=precipavg, transform=ax.transAxes + sublabel_loc,
    #     fontsize=9, fontweight='bold', verticalalignment='top', color = 'blue',
    #     bbox=dict(facecolor='1', edgecolor='none', pad=1.5),zorder=3)
    # if len(perc_assigned) == 4:
    #     ax.text(x=0.755, y=1.0, s=f'{perc_assigned}%', transform=ax.transAxes + sublabel_loc,
    #         fontsize=8, fontweight='bold',verticalalignment='top', color = 'k',
    #         bbox=dict(facecolor=patfreq_col, edgecolor='none', pad=1.5),zorder=3)
    # else:
    #     ax.text(x=0.795, y=1.0, s=f'{perc_assigned}%', transform=ax.transAxes + sublabel_loc,
    #         fontsize=8, fontweight='bold',verticalalignment='top', color = 'k',
    #         bbox=dict(facecolor=patfreq_col, edgecolor='none', pad=1.5),zorder=3)

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
#fig.subplots_adjust(left=0.065,right=0.985,bottom=0.08, top=0.95,hspace=0.05, wspace=0.05) #bottom colorbar
fig.subplots_adjust(left=0.065,right=0.985,bottom=0.082, top=0.985,hspace=0.05, wspace=0.05) #bottom colorbar


save_dir='I:\\Emma\\FIROWatersheds\\Figures\\NodeHistograms'
os.chdir(save_dir)
plt.savefig(f'{percentile}_{numpatterns}_NodeProportionbyNode.png',dpi=300)
plt.show()