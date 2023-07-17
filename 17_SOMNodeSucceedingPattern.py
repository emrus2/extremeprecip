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

#%% IMPORT AVERAGE PRECIP DATA
# avgprecip = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_{numpatterns}NodeAssignAvergagePrecip.npy')
# avgprecip_rounded = np.round(a=avgprecip,decimals=1)

# medprecip = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_{numpatterns}NodeAssignMedianPrecip.npy')
# medprecip_rounded = np.round(a=medprecip,decimals=1)

#%%

#extremedays = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_ExtremeDays.npy')

#assign_day = np.stack((assignment,extremedays),axis=1)
#np.save('I:\\Emma\\FIROWatersheds\\Data\\ExtremeDaysandNodeAssign.npy',assign_day)
#np.save(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_ExtremeDaysand{numpatterns}NodeAssign.npy',assign_day)
dateassign = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_ExtremeDaysand{numpatterns}NodeAssign.npy')
extremeprecip = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_ExtremeDaysPrecip.npy')        

extremeprecip = extremeprecip.reshape(len(extremeprecip),1)
dateassign = np.concatenate((dateassign,extremeprecip), axis=1)


n = 0
alleventdates = []
alleventnodes = []
alleventprecip = []
eventdates = []
eventnodes = []
eventprecip = []
for i, data in enumerate(dateassign):
    assignment = float(data[0])
    assignment_prev = float(dateassign[i-1,0])
    date = int(data[1])
    date_prev = int(dateassign[i-1,1])
    precip = float(data[2])
    precip_prev = float(dateassign[i-1,2])

    if date == date_prev + 1:
        n += 1
        eventdates.append(date_prev)
        eventdates.append(date)
        eventnodes.append(assignment_prev)
        eventnodes.append(assignment)
        eventprecip.append(precip_prev)
        eventprecip.append(precip)
        if len(eventdates) >= 2:
            eventdates.pop()
            eventnodes.pop()
            eventprecip.pop()
    else:
        eventdates.append(date_prev)
        eventnodes.append(assignment_prev)
        eventprecip.append(precip_prev)
        alleventdates.append(eventdates)
        alleventnodes.append(eventnodes)
        alleventprecip.append(eventprecip)
        n = 0
        eventdates = []
        eventnodes = []
        eventprecip = []

'''these have been saved'''
#np.save(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_{numpatterns}NodeClusteredEvents.npy',alleventdates)
#np.save(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_{numpatterns}NodeClusteredEventsNodes.npy',alleventnodes)
#np.save(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_{numpatterns}NodeClusteredEventsPrecip.npy',alleventprecip)


#%% SUBPLOTS OF NODES
colors = ('tomato','indianred','gold','lightgreen','mediumseagreen','cornflowerblue','royalblue','plum','darkorchid','magenta','tomato','indianred','gold','lightgreen','mediumseagreen','cornflowerblue','royalblue','plum','darkorchid','magenta')

fig = plt.figure(figsize=(7.9,5.4))
#fig.suptitle('Event Node Succession',fontsize=13,fontweight="bold",y=0.9875)

for s in np.arange(numpatterns):
    print(s)
    finaldestination = []
    nodeduration = []
    precipaccum = []
    avgprecipaccum = []
    perc_assigned = round(pat_prop[s]*100,1)
    # precipavg = avgprecip_rounded[s]
    # precipmed = medprecip_rounded[s]
    ax = fig.add_subplot(3,3,s+1)
    ax.set_zorder(2)
    ax.set_facecolor('none')
    sublabel_loc = mtransforms.ScaledTranslation(4/72, -4/72, fig.dpi_scale_trans)
    ax.text(x=0.0, y=1.0, s=s+1, transform=ax.transAxes + sublabel_loc,
        fontsize=9, fontweight='bold', verticalalignment='top',
        bbox=dict(facecolor='1', alpha = 0, edgecolor='none', pad=1.5),zorder=3)
    ax.text(x=0.77, y=1.0, s=f'{perc_assigned}%', transform=ax.transAxes + sublabel_loc,
        fontsize=9, fontweight='bold', verticalalignment='top', color = 'red',
        bbox=dict(facecolor='1', alpha=0,edgecolor='none', pad=1.5),zorder=3)
    # ax.text(x=0.4, y=1.0, s=precipavg, transform=ax.transAxes + sublabel_loc,
    #     fontsize=9, fontweight='bold', verticalalignment='top', color = 'blue',
    #     bbox=dict(facecolor='1', edgecolor='none', pad=1.5),zorder=3)
    # ax.text(x=0.4, y=0.9, s=precipmed, transform=ax.transAxes + sublabel_loc,
    #     fontsize=9, fontweight='bold', verticalalignment='top', color = 'g',
    #     bbox=dict(facecolor='1', edgecolor='none', pad=1.5),zorder=3)
    
    for i,l in enumerate(alleventnodes):
        eventlist = alleventdates[i]
        eventprecip = alleventprecip[i]
        #if starting node is the designated node pattern
        if l[0] == s + 1:
            nodeduration.append(len(l))
            xs = np.arange(1,len(l)+1)
            # calculate accumulated precipitation
            if len(l) > 1:
                eventprecipaccum = [sum(eventprecip[:y]) for y in range(1, len(eventprecip) + 1)]
                precipaccum.append(eventprecipaccum)
            else:
                precipaccum.append(eventprecip)
            
            ax.plot(np.arange(1,len(l)+1),l,color = colors[s],marker='o',markersize=5,zorder=2)
            finaldestination.append([xs[-1],l[-1]])
    maxnodeduration = max(nodeduration)
    
    #calculate mean accumulated precipitation
    avgprecip_accum = []
    for j in range(maxnodeduration):
        precip_forday = []
        for k, preciplist in enumerate(precipaccum):
            if len(preciplist) >= j + 1:
                precip_day = preciplist[j]
                precip_forday.append(precip_day)
        avgprecip_forday = np.mean(precip_forday)
        avgprecip_accum.append(avgprecip_forday)
                    
    destinationcount = [[x,finaldestination.count(x)] for x in finaldestination]
    nodedurationaverage = np.mean(nodeduration)
    #nodedurationmed = np.median(nodeduration)
    plt.axvline(nodedurationaverage,color='magenta',alpha=0.8,zorder=1)
    #plt.axvline(nodedurationmed,color='blueviolet',alpha=0.8,zorder=1)

    
    ax.annotate(f'{round(nodedurationaverage,2):.2f}',xy=(5.5,0.9),color='magenta',fontsize='9',fontweight='bold',zorder=3)
    #ax.annotate(f'{round(nodedurationmed,2):.2f}',xy=(5.5,1.9),color='blueviolet',fontsize='9',fontweight='bold',zorder=3)
    for destinations, count in destinationcount:
        ax.annotate(count,xy=destinations,fontsize='10',fontweight='bold',zorder=3)
            
    ax.set_ylim(0.5, numpatterns+0.8)
    ax.set_xlim(0.5,6.5)
    ylabels = np.arange(1,numpatterns+1,1)
    xlabels = np.arange(1,7,1)
    if s == 0 or s == 3:
        ax.set_yticks(ylabels)
        ax.set_xticks(xlabels, [])
        if s == 3:
            ax.set_ylabel('Node',fontweight='bold')
    elif s == 7 or s == 8:
        ax.set_yticks(ylabels,[])
        ax.set_xticks(xlabels)
        if s == 7:
            ax.set_xlabel('Duration (Days)',fontweight='bold')
    elif s == 6:
        #ax.set_ylabel('Node',fontweight='bold')
        ax.set_yticks(ylabels)
        ax.set_xticks(xlabels)
        #ax.set_xlabel('Duration (Days)',fontweight='bold')
    else:
        ax.set_yticks(ylabels,[])
        ax.set_xticks(xlabels, [])
        
    # plot accumulated precipitation
    ax2 = ax.twinx()
    ax2.set_zorder(1)
    ax2.set_facecolor('none')
    ax2.bar(np.arange(1,maxnodeduration+1),avgprecip_accum,color='slategrey',alpha=0.4,width=1)
    ax2.set_ylim(0.5, 650)
    ylabels = np.arange(0,670,150)
    if s == 2 or s == 5 or s == 8:
        ax2.set_yticks(ylabels)
        if s == 5:
            ax2.set_ylabel('Precipitation (mm)',fontweight='bold')
    else:
        ax2.set_yticks(ylabels,[])


#CUSTOMIZE SUBPLOT SPACING
fig.subplots_adjust(left=0.051,right=0.93,bottom=0.08, top=0.98,hspace=0.05, wspace=0.05) #bottom colorbar


save_dir='I:\\Emma\\FIROWatersheds\\Figures\\NodeHistograms'
os.chdir(save_dir)
plt.savefig(f'{percentile}_{numpatterns}_NodeSuccessionbyNodeAccumPrecip.png',dpi=300)
plt.show()