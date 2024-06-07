# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 13:49:49 2023

@author: emrus2

Calculates the monthly frequency of each node, saves as:
    
    I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_{numpatterns}NodeAssignMonthlyHistogram.npy

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
metvar = 'IVT'
numpatterns = 9
percentile = 90

# change directory and import SOM data from .mat file
mat_dir='I:\\Emma\\FIROWatersheds\\Data\\SOMs\\SomOutput'
os.chdir(mat_dir)
soms = scipy.io.loadmat(f'{metvar}_{percentile}_{numpatterns}_sompatterns.mat')
pat_prop = np.squeeze(soms['pat_prop'])

#%% IMPORT EVENT DATA
data_dir='I:\\Emma\\FIROWatersheds\\Data\\'
os.chdir(data_dir)
# import extreme ndays and node assignment rray
nodeassign = np.load(f'{percentile}Percentile_ExtremeDaysand{numpatterns}NodeAssign.npy')
#%% IDENTIFY SUCCEEDING NODE DAYS

def somsucceedinghist(node):
#for node in np.arange(1):
    somsucceeding = []
    # loop through each extreme day
    for i, data in enumerate(nodeassign):
        assignmentnext = float(data[0])
        assignment = float(nodeassign[i-1,0])
        datenext = int(data[1])
        date = int(nodeassign[i-1,1])
        # if the days are succeeding
        if datenext == date + 1:
            # and if the first day is assigned to node
            if assignment == node + 1:
                somsucceeding.append(assignmentnext)
    #somsucceedingcount = [[x,somsucceeding.count(x)] for x in somsucceeding]
    return(somsucceeding) 
            
#%% SUBPLOTS OF NODES
colors = ('tomato','indianred','gold','lightgreen','mediumseagreen','cornflowerblue','royalblue','plum','darkorchid','magenta','tomato','indianred','gold','lightgreen','mediumseagreen','cornflowerblue','royalblue','plum','darkorchid','magenta')
binlist = np.arange(1,11,1)

fig = plt.figure(figsize=(7.2,5))
#fig.suptitle('Event Node Succession',fontsize=13,fontweight="bold",y=0.9875)

for node in np.arange(numpatterns):
    #run function to define succeeding patterns
    somsucceeding = somsucceedinghist(node)
    print(somsucceeding)
    ax = fig.add_subplot(3,3,node+1)
    sublabel_loc = mtransforms.ScaledTranslation(4/72, -4/72, fig.dpi_scale_trans)
    # add node labels
    ax.text(x=0.0, y=1.0, s=node+1, transform=ax.transAxes + sublabel_loc,
        fontsize=10, fontweight='bold', verticalalignment='top',
        bbox=dict(facecolor='1', alpha = 0, edgecolor='none', pad=1.5),zorder=3)
    # add percentage labels
    # ax.text(x=0.77, y=1.0, s=f'{perc_assigned}%', transform=ax.transAxes + sublabel_loc,
    #     fontsize=9, fontweight='bold', verticalalignment='top', color = 'red',
    #     bbox=dict(facecolor='1', alpha=0,edgecolor='none', pad=1.5),zorder=3)
    # loop through events
    freq, bins = np.histogram(somsucceeding,bins=binlist)
    freqprop = freq/len(somsucceeding)*100
    histog = ax.bar(np.arange(1,10,1), height=freqprop,color=colors[node],width=0.8)
                       #,bins=binlist,align='left',color=colors[node],rwidth=0.8)
    ax.plot(node+1,2.3,'*',color='k')
    # # adjust plotting       
    ax.set_ylim(0, 70)
    ax.set_xlim(0.4,9.6)
    xlabels = np.arange(1,numpatterns+1,1)
    ylabels = np.arange(0,70,15)
    if node == 0 or node == 3:
        ax.set_yticks(ylabels)
        ax.set_xticks(xlabels, [])
        if node == 3:
            ax.set_ylabel('Days (%)',fontweight='bold')
    elif node == 7 or node == 8:
        ax.set_yticks(ylabels,[])
        ax.set_xticks(xlabels)
        if node == 7:
            ax.set_xlabel('Node',fontweight='bold')
    elif node == 6:
        ax.set_yticks(ylabels)
        ax.set_xticks(xlabels)
    else:
        ax.set_yticks(ylabels,[])
        ax.set_xticks(xlabels, [])
    ax.tick_params(direction='in',which='both',axis='y')
    #ax.grid(visible=True,axis='y')
        
    # # plot accumulated precipitation
    # ax2 = ax.twinx()
    # ax2.set_zorder(1)
    # ax2.set_facecolor('none')
    # ax2.bar(np.arange(1,maxnodeduration+1),avgprecip_accum,color='slategrey',alpha=0.4,width=1)
    # ax2.set_ylim(0.5, 650)
    # ylabels = np.arange(0,670,100)
    # if node in np.arange(2,9,3):
    #     ax2.set_yticks(ylabels,minor='true')
    #     if node == 5:
    #         ax2.set_ylabel('Precipitation (mm)',fontweight='bold')
    # else:
    #     ax2.set_yticks(ylabels,[],minor='true')
    #     ax2.set_yticklabels([])
    # ax2.tick_params(direction='in',which='both',axis='y')

#CUSTOMIZE SUBPLOT SPACING
fig.subplots_adjust(left=0.065,right=0.985,bottom=0.082, top=0.985,hspace=0.05, wspace=0.05) #bottom colorbar


save_dir='I:\\Emma\\FIROWatersheds\\Figures\\NodeHistograms'
os.chdir(save_dir)
plt.savefig(f'{percentile}_{numpatterns}_SingleNodeSuccessionbyNode.png',dpi=300)
plt.show()