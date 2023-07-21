# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 13:49:49 2023

@author: emrus2

Precipitation histograms for each node

UPDATED 7/11/2023
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
#%% IMPORT AVERAGE PRECIP DATA
avgprecip = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_{numpatterns}NodeAssignAvergagePrecip.npy')
avgprecip_rounded = np.round(a=avgprecip,decimals=1)

medprecip = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_{numpatterns}NodeAssignMedianPrecip.npy')
medprecip_rounded = np.round(a=medprecip,decimals=1)

allprecip = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_{numpatterns}NodeAssignAllPrecip.npy',allow_pickle=True)
#%% SUBPLOTS OF NODES

#months = ("O", "N", "D", "J", "F", "M")
#color = 'tomato'
colors = ('tomato','indianred','gold','lightgreen','mediumseagreen','cornflowerblue','royalblue','plum','darkorchid','magenta',)
binlist = [np.arange(50,240,10)]


#x = np.arange(len(months))  # the label locations
width = 0.8  # the width of the bars
multiplier = 0

fig = plt.figure(figsize=(7.2,5))
#fig.suptitle('Node Precipitation',fontsize=13,fontweight="bold",y=0.9875)

if percentile == 95:
    ymin = 0
    ymax = 11
    yint = 2
elif percentile == 90:
    xmin = 45 
    xmax = 215
    xint = 25
    ymin = 0
    ymax = 17.5
    yint = 3
    
ylabels = np.arange(ymin,ymax-1,yint)
xlabels = np.arange(xmin+5,xmax,xint)

for i,node in enumerate(allprecip):
    #print(i,node)
    perc_assigned = round(pat_prop[i]*100,1)
    precipavg = avgprecip_rounded[i]
    precipmed = medprecip_rounded[i]
    #create colors
    precipmednorm = 0.8-(((precipmed-min(medprecip_rounded))/(max(medprecip_rounded)+7-min(medprecip_rounded))))
    precipmed_col = (precipmednorm,1,precipmednorm)
    precipavgnorm = 0.8-(((precipavg-min(avgprecip_rounded))/(max(avgprecip_rounded)+7-min(avgprecip_rounded))))
    precipavg_col = (precipavgnorm,precipavgnorm,1)
    print(precipavg_col)
    offset = width * multiplier
    ax = fig.add_subplot(3,3,i+1)
    sublabel_loc = mtransforms.ScaledTranslation(4/72, -4/72, fig.dpi_scale_trans)
    ax.text(x=0.0, y=1.0, s=i+1, transform=ax.transAxes + sublabel_loc,
        fontsize=10, fontweight='bold', verticalalignment='top',
        bbox=dict(facecolor='1', alpha = 0.8, edgecolor='none', pad=1.5),zorder=3)
    # ax.text(x=0.77, y=1.0, s=f'{perc_assigned}%', transform=ax.transAxes + sublabel_loc,
    #     fontsize=9, fontweight='bold', verticalalignment='top', color = 'red',
    #     bbox=dict(facecolor='1', edgecolor='none', pad=1.5),zorder=3)
    ax.text(x=0.825, y=1.0, s=precipavg, transform=ax.transAxes + sublabel_loc,
        fontsize=8, fontweight='bold', verticalalignment='top', color = 'k',
        bbox=dict(facecolor=precipavg_col, edgecolor='none', pad=1.5),zorder=3)
    ax.text(x=0.825, y=0.88, s=precipmed, transform=ax.transAxes + sublabel_loc,
        fontsize=8, fontweight='bold', verticalalignment='top', color = 'k',
        bbox=dict(facecolor=precipmed_col, edgecolor='none', pad=1.5),zorder=3)
    freqs = plt.hist(node,color=colors[i],label=i+1,align='mid',bins=binlist[0])
    plt.axvline(x=precipavg,color=(0,0,1))
    plt.axvline(x=precipmed,color=(0,1,0))
    
    #ax.annotate(precipavg,xy=(precipavg+3,14),color='blue',fontsize='9',fontweight='bold',zorder=3)
    #ax.annotate(precipmed,xy=(precipavg+3,12),color='magenta',fontsize='9',fontweight='bold',zorder=4)

    
    
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(xmin,xmax)
    if i == 0 or i == 3:
        ax.set_yticks(ylabels)
        ax.set_xticks(xlabels, [])
        if i == 3:
            ax.set_ylabel('Days',fontweight='bold')
    elif i == 7 or i == 8:
        ax.set_yticks(ylabels,[])
        ax.set_xticks(xlabels)
        if i == 7:
            ax.set_xlabel('Precipitation (mm)',fontweight='bold')
    elif i == 6:
        ax.set_yticks(ylabels)
        ax.set_xticks(xlabels)
    else:
        ax.set_yticks(ylabels,[])
        ax.set_xticks(xlabels, [])
    ax.tick_params(direction='in',which='both',axis='y')


#CUSTOMIZE SUBPLOT SPACING
#fig.subplots_adjust(left=0.065,right=0.985,bottom=0.08, top=0.95,hspace=0.05, wspace=0.05) #bottom colorbar
fig.subplots_adjust(left=0.065,right=0.985,bottom=0.082, top=0.985,hspace=0.05, wspace=0.05) #bottom colorbar


save_dir='I:\\Emma\\FIROWatersheds\\Figures\\NodeHistograms'
os.chdir(save_dir)
plt.savefig(f'{percentile}_{numpatterns}_PrecipbyNode.png',dpi=300)
plt.show()