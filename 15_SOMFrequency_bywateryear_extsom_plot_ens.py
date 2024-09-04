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
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
#from matplotlib.colors import ListedColormap
#import matplotlib.cm as cm
#import matplotlib.transforms as mtransforms
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
#from mpl_toolkits.basemap import Basemap #installed using 
    #conda install -c anaconda basemap
#import scipy.io

#%% IMPORT HISTOGRAM DATA

numpatterns = 12
percentile = 90
clusters = 3

data = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}_{numpatterns}NodeAssignYearlyHistogram_wateryear_{clusters}d.npy')
#data[data == 0] = 'nan'

totaldays = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}_{numpatterns}NodeAssignTotalYearly_wateryear_{clusters}d.npy')

#data1 = np.transpose(data)
#data_prop = data1 / totaldays[:,None]

#data_prop = np.transpose(data_prop)
#data_prop[np.isnan(data_prop)] = 0

#%% IMPORT ENSO TEXT DATA
#enso = pd.read_csv('C:/Users/emrus2/Downloads/meiv2.data.txt',header=None, sep='\s+',index_col=0)
ensoall = pd.read_csv('I:\\Emma\\FIROWatersheds\\data\\ONI_EnsoIndex.txt',header=0,sep='\s+')
enso = ensoall.loc[1980:2022]
#print(enso)
#plt.plot(enso)

# create datetime sequence of dates
years = np.array(enso.index)
datestr = [datetime.strftime(datetime(year,n,1),"%Y%m%d") for year in years for n in range(1,np.shape(enso)[1]+1)]
dates = np.array(datestr)
x = np.arange(len(dates))
enso = enso.to_numpy()
enso = enso.flatten()

ensodata = np.column_stack((dates,enso))
#plt.plot(ensodata[:,0],ensodata[:,1])
# create datetime series
#print(ensodata[:,0])

print(type(years))
fig, ax = plt.subplots(layout='constrained')

plt.fill_between(dates,enso,color=[0.8,0,0],where=enso>0.5)
plt.fill_between(dates, -3,3,where=enso>0.5,color='red',alpha=0.2)
plt.fill_between(dates, -3,3,where=enso<-0.5,color='blue',alpha=0.2)
plt.fill_between(dates,enso,color=[0.6,0.6,1],where=enso<0)
plt.fill_between(dates,enso,color=[0,0,0.9],where=enso<-0.5)
ax.set_ylim(-3,3)
#ax.set_ylabel('Index',fontweight='bold')
#ax.set_xlabel('Year',fontweight='bold')
#ax.set_xticks(range(x[0],x[-1],12),rotation=85)
ax.set_xticks([])
#ax.legend(loc='upper center', ncols=3)
# if numpatterns == 9:
#     ax.legend(loc='upper center',ncols=9,columnspacing=1.4,handletextpad=0.4,fontsize=9.7)
# else:
#     ax.legend(loc='upper center',ncols=12,columnspacing=0.5,handletextpad=0.2,fontsize=9)
#ax.set_xbound(-0.75,len(months)-0.3)
#plt.show()

#plt.fill_between(x,enso,color=[1,0.4,0.4],where=enso>0)
# plt.fill_between(range(len(ensodata)),enso,color=[0.8,0,0],where=enso>0.5)
# plt.fill_between(range(len(enso)), -3,3,where=enso>0.5,color='red',alpha=0.2)
# plt.fill_between(range(len(enso)), -3,3,where=enso<-0.5,color='blue',alpha=0.2)
# plt.fill_between(range(len(enso)),enso,color=[0.6,0.6,1],where=enso<0)
# plt.fill_between(range(len(enso)),enso,color=[0,0,0.9],where=enso<-0.5)
#save_dir='I:\\Emma\\FIROWatersheds\\figures\\NodeHistograms'
#os.chdir(save_dir)
#plt.savefig('enso.png',dpi=300,bbox_inches='tight')
#plt.show()
#enso = np.genfromtxt('C:/Users/emrus2/Downloads/meiv2.data.txt',skip_header=1)
#print(enso.loc[44])
    
#%% CLUSTERED BAR CHART

months = np.arange(1979,2022,1)
if numpatterns == 12:
    colors = ('','','gold','','mediumseagreen','olive','cornflowerblue','royalblue')
    colors = ('dodgerblue','gold','darkorchid','indianred','teal','orange','lightgreen','plum','grey','tomato','olivedrab','mediumslateblue')
    colors = ('tomato','dodgerblue','lightgreen','darkorchid','indianred','blue','gold','mediumslateblue',(.9,0,.9),'slategrey','mediumseagreen','plum','purple')
    colors = ('tomato','cornflowerblue','lightgreen','darkorchid','gold','lightblue','plum','mediumseagreen','indianred','royalblue','grey',(.9,0,.9))
elif numpatterns == 9:
    colors = ('tomato','indianred','gold','lightgreen','cornflowerblue','mediumseagreen','royalblue','plum','darkorchid','magenta',)

x = np.arange(len(months))  # the label locations
#width = 0.3  # the width of the bars


#fig, ax = plt.subplots(layout='constrained')
ax2 = ax.twinx()
bottom = np.zeros(len(months))
   
for i, node in enumerate(data):
    #print(i,node)
    rects = ax2.bar(x, node, label=i+1, align='center', bottom=bottom, color=colors[i])
    #ax.bar_label(rects, padding=3)
    bottom += node
    
# Add some text for labels, title and custom x-axis tick labels, etc.
ax2.set_ylabel('Days',fontweight='bold')
ax2.set_xlabel('Year',fontweight='bold')
#ax.set_title('Node Frequency',fontweight='bold')
reducedyr = [str(month)[2:4] for month in months]
ticklabels = [f"{str(year)[2:4]}-{str(year+1)[2:4]}" for year in months]
ax2.set_xticks(x, labels=ticklabels,rotation=85,fontsize=7.5)
ax2.set_yticks(range(0,20,2))
#ax.legend(loc='upper center', ncols=3)
if numpatterns == 9:
    ax2.legend(loc='upper center',ncols=9,columnspacing=1.4,handletextpad=0.4,fontsize=9.7)
else:
    ax2.legend(loc='upper center',ncols=12,columnspacing=0.5,handletextpad=0.2,fontsize=9)
ax2.set_ylim(0,20)
ax2.set_xbound(-0.75,len(months)-0.3)

box = ax2.get_position()
#[left,bottom,width,height]
#ax.set_position([0.09, 0.13,0.78, 0.8])

save_dir='I:\\Emma\\FIROWatersheds\\Figures\\NodeHistograms'
os.chdir(save_dir)
#plt.savefig(f'{percentile}_{numpatterns}_NodeFreqYearlyNum_wateryear_{clusters}d_test.png',dpi=300,bbox_inches='tight',pad_inches=0.08)
plt.show()