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
import seaborn as sns
#from matplotlib.colors import ListedColormap
#import matplotlib.cm as cm
import matplotlib.transforms as mtransforms
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
#from mpl_toolkits.basemap import Basemap #installed using 
    #conda install -c anaconda basemap
import scipy.io
import xarray as xr
from scipy.stats import gaussian_kde
#%% IMPORT SOM DATA

# define metvar
#metvars = ['Z500', 'SLP', '850T', '300W', 'IVT','Z500Anom','300W']
metvar = 'IVT'
numpatterns = 9
percentile = 90
watershed = "UpperYuba"

# change directory and import SOM data from .mat file
mat_dir='I:\\Emma\\FIROWatersheds\\Data'
os.chdir(mat_dir)
soms = scipy.io.loadmat(f'SOMs\\SomOutput\\{metvar}_{percentile}_{numpatterns}_sompatterns.mat')
pat_prop = np.squeeze(soms['pat_prop'])
pat_freq = np.squeeze(soms['pat_freq'])

#%% IMPORT ALL PRECIPITATION DATA
#define filepath and file of daily mean total precipitation
folderpath = f'Gridmet\\{watershed}\\DailyMeans'
filename = f'{watershed}dailymeansTOTAL_1980_2021.nc'
filepath = os.path.join(folderpath, filename)

#open filepath with xarray
NCdatatotal = xr.open_dataarray(filepath)
print(NCdatatotal)

#reduce total data to wet season data (October-March)
winter = [1, 2, 3, 10, 11, 12]
NCdatawinter =  NCdatatotal.where(NCdatatotal['day.month'].isin(winter), drop=True)
print(NCdatawinter)


#%% IMPORT SOM PRECIP DATA
avgprecip = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_{numpatterns}NodeAssignAvergagePrecip.npy')
avgprecip_rounded = np.round(a=avgprecip,decimals=1)

medprecip = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_{numpatterns}NodeAssignMedianPrecip.npy')
medprecip_rounded = np.round(a=medprecip,decimals=1)

allprecip = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_{numpatterns}NodeAssignAllPrecip.npy',allow_pickle=True)

#%% CALCULATE MAX AND MIN PRECIP VALUES FOR EACH NODE
maxprecip = [max(s) for s in allprecip]
minprecip = [min(s) for s in allprecip]
#%% SUBPLOTS OF NODES
# colors = ('tomato','indianred','gold','lightgreen','mediumseagreen','cornflowerblue','royalblue','plum','darkorchid','magenta',)
# zorders = (-1,-4,-2,-7,-6,-5,-8,-9,-3)

fig, ax = plt.subplots(layout='constrained')
plot = sns.histplot(NCdatawinter, kde=True, bins=75, zorder=1, \
                 color = 'slategrey',stat='frequency',alpha=0.7)

extrlim = 51.53282318
ymax = 21
xmax = 208
# ax.axvline(x=extrlim,color='red',linewidth=2,alpha=0.5)
# plt.fill_betweenx(y=plot, x1=extrlim, x2=xmax,color='red')
for line in ax.lines:
    x, y = line.get_xydata().T
    ax.fill_between(x, 0, y, color='royalblue', where=x<extrlim,alpha=0.5,label = '<90th Percentile')
    ax.fill_between(x, 0, y, color='red', where=x>extrlim,alpha=0.5,label='>90th Percentile')
    # ax.fill_betweenx(y, x, extrlim, color='red')


# for i,node in enumerate(allprecip):
#     line = plt.axvline(x=minprecip[i],ymax=0.9,color=colors[i],linewidth=1.5,alpha=0.8,zorder=zorders[i],label=i+1)
#     plt.axvline(x=maxprecip[i],ymax=0.9,color=colors[i],linewidth=1.5,alpha=0.8,zorder=zorders[i])
#     plt.fill_betweenx(y=np.arange(0,700), x1=minprecip[i], x2=maxprecip[i],color=colors[i],zorder=zorders[i],alpha=0.5)

#CUSTOMIZE SUBPLOT SPACING
ax.set_ylabel('Frequency',fontweight='bold')
ax.set_xlabel('Precipitation (mm)',fontweight='bold')
ax.legend(loc='upper center',ncols=9,columnspacing=1.4,handletextpad=0.4,fontsize=9.7)
# ax.set_ylim(0,ymax)
ax.set_xlim(-1,xmax)


save_dir='I:\\Emma\\FIROWatersheds\\Figures\\NodeHistograms'
os.chdir(save_dir)
plt.savefig(f'{percentile}_{numpatterns}_PrecipDistribution2.png',dpi=300)
plt.show()