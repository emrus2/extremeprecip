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
#import scipy.io
#%% IMPORT SOM DATA
# change directory and import SOM data from .mat file
mat_dir='I:\\Emma\\FIROWatersheds\\Data\\SOMs\\SomOutput'
numpatterns = 9
percentile = 90

#os.chdir(mat_dir)
#soms = scipy.io.loadmat(f'IVT_{percentile}_{numpatterns}_sompatterns.mat')
#pats = np.squeeze(soms['pats'])
#asn_err = np.squeeze(soms['assignment'])
#assignment = asn_err[:,0]

#extremedays = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_ExtremeDays.npy')

#assign_day = np.stack((assignment,extremedays),axis=1)
#np.save('I:\\Emma\\FIROWatersheds\\Data\\ExtremeDaysandNodeAssign.npy',assign_day)
#np.save(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_ExtremeDaysand{numpatterns}NodeAssign.npy',assign_day)
dateassign = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_ExtremeDaysand{numpatterns}NodeAssign.npy')
n = 0
alleventdates = []
alleventnodes = []
eventdates = []
eventnodes = []
for i, data in enumerate(dateassign):
    assignment = float(data[0])
    assignment_prev = float(dateassign[i-1,0])
    date = int(data[1])
    date_prev = int(dateassign[i-1,1])
    
    if date == date_prev + 1:
        #print('succession')
        #print(eventdates)
        n += 1
        #print(len(eventdates))
        
            #print('Popping', eventdates)
        eventdates.append(date_prev)
        eventdates.append(date)
        eventnodes.append(assignment_prev)
        eventnodes.append(assignment)
        #print(eventdates)
        if len(eventdates) >= 2:
            eventdates.pop()
            eventnodes.pop()
        #print(eventdates)
    else:
        #print(eventdates)
        eventdates.append(date_prev)
        eventnodes.append(assignment_prev)
        alleventdates.append(eventdates)
        alleventnodes.append(eventnodes)
        print(eventdates)
        #print(alleventdates)
        n = 0
        eventdates = []
        eventnodes = []
        
#%%
        
# colors = ('tomato','indianred','gold','lightgreen','mediumseagreen','cornflowerblue','royalblue','plum','darkorchid','magenta',)

# fig, ax = plt.subplots(layout='constrained')

# for i,l in enumerate(alleventnodes):
#     plt.plot(np.arange(1,len(l)+1),l,color = colors[int(l[0]-1)],markersize=5,marker='o')
    
# ax.set_ylabel('Node',fontweight='bold')
# ax.set_xlabel('Event Day',fontweight='bold')
# plt.show()

#SUBPLOTS OF NODES
colors = ('tomato','indianred','gold','lightgreen','mediumseagreen','cornflowerblue','royalblue','plum','darkorchid','magenta','tomato','indianred','gold','lightgreen','mediumseagreen','cornflowerblue','royalblue','plum','darkorchid','magenta')
markers = ['>','x','d','s','.','p','>','x','d','s','.','p']

fig = plt.figure(figsize=(7.5,5.4))
fig.suptitle('Event Node Succession',fontsize=13,fontweight="bold",y=0.9875)

for s in np.arange(numpatterns):
    finaldestination = []
    nodeduration = []
    #perc_assigned = round(pat_prop[i]*100,1)
    #precipavg = precip_rounded[i]
    ax = fig.add_subplot(3,3,s+1)
    sublabel_loc = mtransforms.ScaledTranslation(4/72, -4/72, fig.dpi_scale_trans)
    ax.text(x=0.0, y=1.0, s=s+1, transform=ax.transAxes + sublabel_loc,
        fontsize=9, fontweight='bold', verticalalignment='top', 
        bbox=dict(facecolor='1', edgecolor='none', pad=1.5),zorder=3)
    # ax.text(x=0.77, y=1.0, s=f'{perc_assigned}%', transform=ax.transAxes + sublabel_loc,
    #     fontsize=9, fontweight='bold', verticalalignment='top', color = 'red',
    #     bbox=dict(facecolor='1', edgecolor='none', pad=1.5),zorder=3)
    # ax.text(x=0.45, y=1.0, s=precipavg, transform=ax.transAxes + sublabel_loc,
    #     fontsize=9, fontweight='bold', verticalalignment='top', color = 'blue',
    #     bbox=dict(facecolor='1', edgecolor='none', pad=1.5),zorder=3)
    #markerint = 0
    for i,l in enumerate(alleventnodes):
        eventlist = alleventdates[i]
        if l[0] == s + 1:
            nodeduration.append(len(l))
            xs = np.arange(1,len(l)+1)
            plt.plot(np.arange(1,len(l)+1),l,color = colors[s],marker='o',markersize=5,zorder=2)
            
            #markerint += 1
            #print(xs[-1],l[-1])
            finaldestination.append([xs[-1],l[-1]])
    destinationcount = [[x,finaldestination.count(x)] for x in finaldestination]
    nodedurationaverage = np.mean(nodeduration)
    plt.axvline(nodedurationaverage,color='slategrey',alpha=0.8,zorder=1)
    ax.annotate(f'{round(nodedurationaverage,2):.2f}',xy=(5.3,0.9),color='slategrey',fontsize='10',fontweight='bold',zorder=3)
    for destinations, count in destinationcount:
        ax.annotate(count,xy=destinations,fontsize='10',fontweight='bold',zorder=3)
            
    ax.set_ylim(0.5, numpatterns+0.8)
    ax.set_xlim(0.4,6.3)
    ylabels = np.arange(1,numpatterns+1,1)
    xlabels = np.arange(1,7,1)
    if s == 0 or s == 3:
        ax.set_yticks(ylabels)
        ax.set_ylabel('Node',fontweight='bold')
        ax.set_xticks(xlabels, [])
    elif s == 7 or s == 8:
        ax.set_yticks(ylabels,[])
        ax.set_xticks(xlabels)
        ax.set_xlabel('Duration (Days)',fontweight='bold')
    elif s == 6:
        ax.set_ylabel('Node',fontweight='bold')
        ax.set_yticks(ylabels)
        ax.set_xticks(xlabels)
        ax.set_xlabel('Duration (Days)',fontweight='bold')
    else:
        ax.set_yticks(ylabels,[])
        ax.set_xticks(xlabels, [])

#CUSTOMIZE SUBPLOT SPACING
fig.subplots_adjust(left=0.065,right=0.985,bottom=0.08, top=0.95,hspace=0.05, wspace=0.05) #bottom colorbar


save_dir='I:\\Emma\\FIROWatersheds\\Figures\\NodeHistograms'
os.chdir(save_dir)
plt.savefig(f'{percentile}_{numpatterns}_AllEventNodeSuccessionbyNode.png',dpi=300)
plt.show()

#%%
n = 0
alleventdates = []
alleventnodes = []
eventdates = []
eventnodes = []
for i, data in enumerate(dateassign):
    assignment = float(data[0])
    assignment_prev = float(dateassign[i-1,0])
    date = int(data[1])
    date_prev = int(dateassign[i-1,1])
    #print(i, assignment, date)
    
    if date == date_prev + 1:
        n += 1
        if len(eventdates) >= 2:
            eventdates.pop()
            eventnodes.pop()
        eventdates.append(date_prev)
        eventdates.append(date)
        eventnodes.append(assignment_prev)
        eventnodes.append(assignment)
        
    else:
        if len(eventdates) > 0:
            alleventdates.append(eventdates)
            alleventnodes.append(eventnodes)
            print(eventdates)
            print(eventnodes)
        n = 0
        eventdates = []
        eventnodes = []
        
# colors = ('tomato','indianred','gold','lightgreen','mediumseagreen','cornflowerblue','royalblue','plum','darkorchid','magenta',)

# fig, ax = plt.subplots(layout='constrained')

# for i,l in enumerate(alleventnodes):
#     plt.plot(np.arange(1,len(l)+1),l,color = colors[int(l[0]-1)],markersize=5,marker='o')
    
# ax.set_ylabel('Node',fontweight='bold')
# ax.set_xlabel('Event Day',fontweight='bold')
# plt.show()
#%%
colors = ('tomato','indianred','gold','lightgreen','mediumseagreen','cornflowerblue','royalblue','plum','darkorchid','magenta','tomato','indianred','gold','lightgreen','mediumseagreen','cornflowerblue','royalblue','plum','darkorchid','magenta')
markers = ['>','x','d','s','.','p','>','x','d','s','.','p']

fig = plt.figure(figsize=(7.5,5.4))
fig.suptitle('Event Node Succession',fontsize=13,fontweight="bold",y=0.9875)

for s in np.arange(numpatterns):
    finaldestination = []
    nodeduration = []
    #perc_assigned = round(pat_prop[i]*100,1)
    #precipavg = precip_rounded[i]
    ax = fig.add_subplot(3,3,s+1)
    sublabel_loc = mtransforms.ScaledTranslation(4/72, -4/72, fig.dpi_scale_trans)
    ax.text(x=0.0, y=1.0, s=s+1, transform=ax.transAxes + sublabel_loc,
        fontsize=9, fontweight='bold', verticalalignment='top', 
        bbox=dict(facecolor='1', edgecolor='none', pad=1.5),zorder=3)
    # ax.text(x=0.77, y=1.0, s=f'{perc_assigned}%', transform=ax.transAxes + sublabel_loc,
    #     fontsize=9, fontweight='bold', verticalalignment='top', color = 'red',
    #     bbox=dict(facecolor='1', edgecolor='none', pad=1.5),zorder=3)
    # ax.text(x=0.45, y=1.0, s=precipavg, transform=ax.transAxes + sublabel_loc,
    #     fontsize=9, fontweight='bold', verticalalignment='top', color = 'blue',
    #     bbox=dict(facecolor='1', edgecolor='none', pad=1.5),zorder=3)
    #markerint = 0
    for i,l in enumerate(alleventnodes):
        eventlist = alleventdates[i]
        if l[0] == s + 1:
            nodeduration.append(len(l))
            print(nodeduration)
            xs = np.arange(1,len(l)+1)
            plt.plot(np.arange(1,len(l)+1),l,color = colors[s],marker='o',markersize=5,zorder=2)
            
            #markerint += 1
            #print(xs[-1],l[-1])
            finaldestination.append([xs[-1],l[-1]])
    destinationcount = [[x,finaldestination.count(x)] for x in finaldestination]
    print(destinationcount)
    nodedurationaverage = np.mean(nodeduration)
    plt.axvline(nodedurationaverage,color='slategrey',alpha=0.8,zorder=1)
    ax.annotate(round(nodedurationaverage,2),xy=(nodedurationaverage,1.15),color='slategrey',fontsize='10',fontweight='bold',zorder=3)
    for destinations, count in destinationcount:
        ax.annotate(count,xy=destinations,fontsize='10',fontweight='bold',zorder=3)
            
    print(finaldestination)
    ax.set_ylim(0.5, numpatterns+0.75)
    ax.set_xlim(0.4,6.3)
    ylabels = np.arange(1,numpatterns+1,1)
    xlabels = np.arange(1,7,1)
    if s == 0 or s == 3:
        ax.set_yticks(ylabels)
        ax.set_ylabel('Node',fontweight='bold')
        ax.set_xticks(xlabels, [])
    elif s == 7 or s == 8:
        ax.set_yticks(ylabels,[])
        ax.set_xticks(xlabels)
        ax.set_xlabel('Duration (Days)',fontweight='bold')
    elif s == 6:
        ax.set_ylabel('Node',fontweight='bold')
        ax.set_yticks(ylabels)
        ax.set_xticks(xlabels)
        ax.set_xlabel('Duration (Days)',fontweight='bold')
    else:
        ax.set_yticks(ylabels,[])
        ax.set_xticks(xlabels, [])

#CUSTOMIZE SUBPLOT SPACING
fig.subplots_adjust(left=0.065,right=0.985,bottom=0.08, top=0.95,hspace=0.05, wspace=0.05) #bottom colorbar


save_dir='I:\\Emma\\FIROWatersheds\\Figures\\NodeHistograms'
os.chdir(save_dir)
plt.savefig(f'{percentile}_{numpatterns}_EventNodeSuccessionbyNode.png',dpi=300)
plt.show()