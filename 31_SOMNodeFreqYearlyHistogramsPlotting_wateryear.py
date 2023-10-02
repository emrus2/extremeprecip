# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 13:49:49 2023

@author: emrus2

Plots the annual composition of node patterns


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
#import scipy.io

#%% IMPORT HISTOGRAM DATA
numpatterns = 9
percentile = 90
data = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}_{numpatterns}NodeAssignYearlyHistogram_wateryear.npy')
#data[data == 0] = 'nan'

totaldays = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}_{numpatterns}NodeAssignTotalYearly_wateryear.npy')

#data1 = np.transpose(data)
#data_prop = data1 / totaldays[:,None]

#data_prop = np.transpose(data_prop)
#data_prop[np.isnan(data_prop)] = 0
    
#%% CLUSTERED BAR CHART

months = np.arange(1979,2022,1)
if numpatterns == 12:
    colors = ('orangered','indianred','gold','lightgreen','mediumseagreen','olive','cornflowerblue','royalblue','plum','darkorchid','magenta','slategrey')
elif numpatterns == 9:
    colors = ('tomato','indianred','gold','lightgreen','mediumseagreen','cornflowerblue','royalblue','plum','darkorchid','magenta',)

x = np.arange(len(months))  # the label locations
#width = 0.3  # the width of the bars


fig, ax = plt.subplots(layout='constrained')
bottom = np.zeros(len(months))
   
for i, node in enumerate(data):
    #print(i,node)
    rects = ax.bar(x, node, label=i+1, align='center', bottom=bottom, color=colors[i])
    #ax.bar_label(rects, padding=3)
    bottom += node
    
# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Days',fontweight='bold')
ax.set_xlabel('Year',fontweight='bold')
#ax.set_title('Node Frequency',fontweight='bold')
reducedyr = [str(month)[2:4] for month in months]
ticklabels = [f"{str(year)[2:4]}-{str(year+1)[2:4]}" for year in months]
ax.set_xticks(x, labels=ticklabels,rotation=70,fontsize=7.5)
ax.set_yticks(range(0,21,2))
#ax.legend(loc='upper center', ncols=3)
ax.legend(loc='upper center',ncols=9,columnspacing=1.4,handletextpad=0.4,fontsize=9.7)
ax.set_ylim(0,20.5)
ax.set_xbound(-0.75,len(months)-0.3)

box = ax.get_position()
#[left,bottom,width,height]
#ax.set_position([0.09, 0.13,0.78, 0.8])

save_dir='I:\\Emma\\FIROWatersheds\\Figures\\NodeHistograms'
os.chdir(save_dir)
plt.savefig(f'{percentile}_{numpatterns}_NodeFreqYearlyNum_wateryear.png',dpi=300,bbox_inches='tight',pad_inches=0.08)
plt.show()