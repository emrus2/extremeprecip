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
alleventdates = np.load(f'{percentile}Percentile_{numpatterns}NodeClusteredEvents.npy',allow_pickle=True)
alleventnodes = np.load(f'{percentile}Percentile_{numpatterns}NodeClusteredEventsNodes.npy',allow_pickle=True)
alleventprecip = np.load(f'{percentile}Percentile_{numpatterns}NodeClusteredEventsPrecip.npy',allow_pickle=True)

#%% CALCULATE THE PATH FREQUENCY TO PLOT WIDTH BASED ON IT

# create list with node event paths
def pathfreq(node):
    nodepaths = []
    nodecounts = []
    for ent, k in enumerate(alleventnodes):
        if k[0] == node + 1:
            if len(k) > 1:
                nodepaths.append(k)
                
    # calculate longest duration event           
    longestevent = max(map(len,nodepaths))
    
    # expand all events to length of longest event
    for i, ev in enumerate(nodepaths):
        while len(ev) < longestevent:
            ev.append(0)
        #print(ev)
        
    for j in np.arange(longestevent):
        daynodes = []
        for i,ev in enumerate(nodepaths):
            daynodes.append(ev[j])
        nodecount = [[x,daynodes.count(x)] for x in daynodes]
        new_count = []
        for elem in nodecount:
            if elem not in new_count:
                if 0 not in elem:
                    new_count.append(elem)
        nodecounts.append(new_count)
    return(nodecounts)    
#%% SUBPLOTS OF NODES
colors = ('tomato','indianred','gold','lightgreen','mediumseagreen','cornflowerblue','royalblue','plum','darkorchid','magenta','tomato','indianred','gold','lightgreen','mediumseagreen','cornflowerblue','royalblue','plum','darkorchid','magenta')

fig = plt.figure(figsize=(7.2,5))
#fig.suptitle('Event Node Succession',fontsize=13,fontweight="bold",y=0.9875)

for node in np.arange(numpatterns):
    #print(node)
    counts = pathfreq(node)
    finaldestination = []
    nodeduration = []
    precipaccum = []
    avgprecipaccum = []
    #perc_assigned = round(pat_prop[node]*100,1)
    ax = fig.add_subplot(3,3,node+1)
    ax.set_zorder(2)
    ax.set_facecolor('none')
    sublabel_loc = mtransforms.ScaledTranslation(4/72, -4/72, fig.dpi_scale_trans)
    # define color of pattern frequency
    # add node labels
    ax.text(x=0.0, y=1.0, s=node+1, transform=ax.transAxes + sublabel_loc,
        fontsize=10, fontweight='bold', verticalalignment='top',
        bbox=dict(facecolor='1', alpha = 0, edgecolor='none', pad=1.5),zorder=3)
    # add percentage labels
    # ax.text(x=0.77, y=1.0, s=f'{perc_assigned}%', transform=ax.transAxes + sublabel_loc,
    #     fontsize=9, fontweight='bold', verticalalignment='top', color = 'red',
    #     bbox=dict(facecolor='1', alpha=0,edgecolor='none', pad=1.5),zorder=3)
    # loop through events
    for event,assign in enumerate(alleventnodes):
        assign = [a for a in assign if a != 0]
        # if starting node is the designated node pattern
        if assign[0] == node + 1:
            #print(assign)
            # define event dates and precip
            eventlist = alleventdates[event]
            eventprecip = alleventprecip[event]
            # count duration of event
            nodeduration.append(len(assign))
            # create x values
            xlist = np.arange(1,len(assign)+1)
            # calculate accumulated precipitation
            if len(assign) > 1:
                eventprecipaccum = [sum(eventprecip[:y]) for y in range(1, len(eventprecip) + 1)]
                precipaccum.append(eventprecipaccum)
            else:
                precipaccum.append(eventprecip)
            # calculate final destination of event
            endday = xlist[-1]
            endnode = assign[-1]
            finaldestination.append([endday,endnode])
            # define line size
            dayendcount = counts[xlist[-1]-1]
            for lst in dayendcount:
                if endnode == lst[0]:
                    linewidth = lst[1]
            #linewidth = 1
            # plot event
            ax.plot(np.arange(1,len(assign)+1),assign,color = colors[node],marker='o',markersize=5,linewidth = linewidth,zorder=2)
            
    # calculate longest lasting event in node
    maxnodeduration = max(nodeduration)
    
    # calculate mean accumulated precipitation
    avgprecip_accum = []
    # loop through event days
    for day in range(maxnodeduration):
        precip_forday = [] # list to store day precipitaiton
        for pr, precipevent in enumerate(precipaccum):
            if len(precipevent) >= day + 1:
                precip_day = precipevent[day]
                precip_forday.append(precip_day)
        avgprecip_forday = np.mean(precip_forday)
        avgprecip_accum.append(avgprecip_forday)
    
    # count the number of items in final destination                
    destinationcount = [[x,finaldestination.count(x)] for x in finaldestination]
    for destinations, count in destinationcount:
        ax.annotate(count,xy=destinations,fontsize='8',fontweight='bold',zorder=3)
    
    # calculate average event duration
    nodedurationaverage = np.mean(nodeduration)
    #plt.axvline(nodedurationaverage,color='magenta',alpha=0.8,zorder=1)
    #ax.annotate(f'{round(nodedurationaverage,2):.2f}',xy=(5.5,0.9),color='magenta',fontsize='9',fontweight='bold',zorder=3)
    nodedurnorm = 1.6 - nodedurationaverage
    nodedurnorm = 0.9-(((nodedurationaverage-1.13)/0.5))
    nodedur_col = (1,nodedurnorm,1)
    print(nodedur_col)
    ax.text(x=0.795, y=0.14, s=f'{round(nodedurationaverage,2):.2f}', transform=ax.transAxes + sublabel_loc,
        fontsize=8, fontweight='bold', verticalalignment='top', color = 'k',
        bbox=dict(facecolor=nodedur_col, edgecolor='none', pad=1.5),zorder=3)
    # plot number of events
    numevents = len(nodeduration)
    gbcolor = 1.1-(numevents/37)
    patfreq_col = (1,gbcolor,gbcolor)    
    ax.text(x=0.855, y=0.26, s=numevents, transform=ax.transAxes + sublabel_loc,
        fontsize=8, fontweight='bold',verticalalignment='top', color = 'k',
        bbox=dict(facecolor=patfreq_col, edgecolor='none', pad=1.5),zorder=3)
    # adjust plotting       
    ax.set_ylim(0.5, numpatterns+0.8)
    ax.set_xlim(0.5,6.5)
    ylabels = np.arange(1,numpatterns+1,1)
    xlabels = np.arange(1,7,1)
    if node in np.arange(6):
        ax.set_yticks(ylabels)
        ax.set_xticks(xlabels, [])
        if node == 3:
            ax.set_ylabel('Node',fontweight='bold')
    else:
        ax.set_yticks(ylabels)
        ax.set_xticks(xlabels)
        if node == 7:
            ax.set_xlabel('Duration (Days)',fontweight='bold')
    ax.tick_params(direction='in',which='both',axis='y')
    #ax.grid(visible=True,axis='y')
        
    # plot accumulated precipitation
    ax2 = ax.twinx()
    ax2.set_zorder(1)
    ax2.set_facecolor('none')
    ax2.bar(np.arange(1,maxnodeduration+1),avgprecip_accum,color='slategrey',alpha=0.4,width=1)
    ax2.set_ylim(0.5, 650)
    ylabels = np.arange(0,670,100)
    if node in np.arange(2,9,3):
        ax2.set_yticks(ylabels,minor='true')
        if node == 5:
            ax2.set_ylabel('Precipitation (mm)',fontweight='bold')
    else:
        ax2.set_yticks(ylabels,[],minor='true')
        ax2.set_yticklabels([])
    ax2.tick_params(direction='in',which='both',axis='y')

#CUSTOMIZE SUBPLOT SPACING
fig.subplots_adjust(left=0.051,right=0.93,bottom=0.085, top=0.985,hspace=0.05, wspace=0.09) #bottom colorbar


save_dir='I:\\Emma\\FIROWatersheds\\Figures\\NodeHistograms'
os.chdir(save_dir)
plt.savefig(f'{percentile}_{numpatterns}_NodeSuccessionbyNodeAccumPrecip_NumEvents.png',dpi=300)
plt.show()