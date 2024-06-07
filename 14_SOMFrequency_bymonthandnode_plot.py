# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 13:49:49 2023

@author: emrus2

Uses:
    I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_{numpatterns}NodeAssignMonthlyHistogram.npy
    
to plot the monthly node frequency (for cummulative nodes)

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
clusters = 3
lettered = True

# change directory and import SOM data from .mat file
mat_dir='I:\\Emma\\FIROWatersheds\\Data\\SOMs\\SomOutput'
os.chdir(mat_dir)
soms = scipy.io.loadmat(f'{metvar}_{percentile}_{numpatterns}sompatterns_{clusters}d.mat')
pat_prop = np.squeeze(soms['pat_prop'])
pat_freq = np.squeeze(soms['pat_freq'])
#%% IMPORT HISTOGRAM DATA
data = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_{numpatterns}NodeAssignMonthlyHistogram_{clusters}d.npy')
#data[data == 0] = 'nan'
totaldata = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_{numpatterns}MonthlyTotals_{clusters}d.npy')

#convert data to a percentage of monthly total
datapercent = data / totaldata * 100
#%% CLUSTERED BAR CHART - 12 NODE

## SET PLOT A) PARAMETERS
#define month labels and colors
months = ("Oct", "Nov", "Dec", "Jan", "Feb", "Mar")
colors = ('tomato','cornflowerblue','lightgreen','darkorchid','gold','lightblue','plum','mediumseagreen','indianred','royalblue','grey',(.9,0,.9))
ymax = 110 # max y 
yint = 20 # y interval 
x = np.arange(len(months))  # the label locations
width = 1  # the width of the bars
# define number of columns and rows for subplot
numcols = 3
numrows = 6
# set heights for each row
rowheights = [2.5, 0.2, 1, 1, 1, 1]
# create figure
fig, axs = plt.subplots(numrows, numcols, figsize=(7, 11), height_ratios=rowheights)
# remove all axes so they aren't seen under plots
for i in range(numrows):
    for j in range(numcols):
        axs[i,j].axis('off')
        
## PLOT A) DATA
ax = fig.add_subplot(numrows,numcols,(1,3)) #add subplot
bottom = np.zeros(len(months)) # set bottom 'counters' for stacked histogram
fig.suptitle('a)',fontsize=11,fontweight="bold",y=0.99,x=0.03) # add title
# plot stacked data
for i, node in enumerate(datapercent):
    rects = ax.bar(x, node, label=i+1, align='center', bottom=bottom, color=colors[i])
    bottom += node
# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Percentage of Days in Month',fontweight='bold',labelpad=0)
ax.yaxis.set_label_coords(-0.045,0.5)
ax.set_xlabel('Month',fontweight='bold',labelpad=0)
ax.set_xticks(x, months)
ax.legend(loc='upper center',ncols=numpatterns,columnspacing=0.7,handletextpad=0.25,fontsize=10)
ax.set_ylim(0, ymax)
ax.set_yticks(np.arange(0,ymax,yint))
ax.set_xlim(-0.5,5.5)
ax.tick_params(direction='in',which='both',axis='y')


## SET PLOT B) PARAMETERS
months_abb = ("O", "N", "D", "J", "F", "M") # x labels
width = 0.8  # the width of the bars
if clusters == 1:
    ymax = 50
    yint = 10
elif clusters == 5:
    ymax = 71
    yint = 15
elif clusters == 3:
    ymax = 70
    yint = 15

## PLOT B) DATA
for i,node in enumerate(data):
    # perc_assigned = str(round(pat_prop[i]*100,1))
    node_rel = (node / pat_freq[i]) * 100
    round_node_rel = [ round(elem, 0) for elem in node_rel]
    # add subplots
    ax = fig.add_subplot(numrows,numcols,i+7)
    sublabel_loc = mtransforms.ScaledTranslation(4/72, -4/72, fig.dpi_scale_trans)
    ax.text(x=0.0, y=1.0, s=i+1, transform=ax.transAxes + sublabel_loc,
        fontsize=10, fontweight='bold', verticalalignment='top', 
        bbox=dict(facecolor='1', edgecolor='none', pad=1.5),zorder=3)
    freqs = plt.bar(x,height=round_node_rel,width=width,color=colors[i],label=i+1,align='center')
    ax.bar_label(freqs, padding=3)
    ax.set_ylim(0, ymax)
    ax.set_xlim(-0.5,5.5)
    ylabels = np.arange(0,ymax,yint)
    if i == 0:
        ax.set_title('b)',fontsize=11,fontweight="bold",y=0.865,x=-0.17)
    if i in range(0,7,3):
        ax.set_yticks(ylabels)
        ax.set_xticks(x, [])
        if i == 3:
            ax.set_ylabel('Percentage of Days in Node',fontweight='bold')
            ax.yaxis.set_label_coords(-0.135,0)
    elif i == numpatterns-2 or i == numpatterns-1:
        ax.set_yticks(ylabels,[])
        ax.set_xticks(x, months_abb)
        if i == 10:
            ax.set_xlabel('Month',fontweight='bold')
    elif i == numpatterns-3:
        ax.set_yticks(ylabels)
        ax.set_xticks(x, months_abb)
    else:
        ax.set_yticks(ylabels,[])
        ax.set_xticks(x, [])
    ax.tick_params(direction='in',which='both',axis='y')

#CUSTOMIZE SUBPLOT SPACING
fig.subplots_adjust(left=0.07,right=0.99,bottom=0.04, top=0.99,hspace=0.05, wspace=0.05) #bottom colorbar
save_dir='I:\\Emma\\FIROWatersheds\\Figures\\NodeHistograms'
os.chdir(save_dir)
# plt.savefig(f'{percentile}_{numpatterns}_NodeFreqMonthlyandIndiv_{clusters}d.png',dpi=300)
plt.show()