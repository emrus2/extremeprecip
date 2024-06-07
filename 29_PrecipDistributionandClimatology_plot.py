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
import netCDF4 as nc
from datetime import datetime, timedelta
from matplotlib.legend_handler import HandlerTuple
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
    fillblue = ax.fill_between(x, 0, y, color='b', where=x<=extrlim,alpha=0.5,label = '<90th Percentile')
    fillred = ax.fill_between(x, 0, y, color='red', where=x>extrlim,alpha=0.5,label='>90th Percentile')
    # ax.fill_betweenx(y, x, extrlim, color='red')

#CUSTOMIZE SUBPLOT SPACING
ax.set_ylabel('Days',fontweight='bold')
ax.set_xlabel('Precipitation (mm)',fontweight='bold')
ax.legend(loc='lower right',bbox_to_anchor=(0.5, 0.01, 0.5, 0.5),ncols=9,columnspacing=1.4,handletextpad=0.4,fontsize=9.7)

# ax.set_ylim(0,ymax)
ax.set_xlim(-1,xmax)

#%%add plot within plot
inset_ax = fig.add_axes([0.305, .374, .68, .61],frame_on=True)
bbox = inset_ax.get_tightbbox(fig.canvas.get_renderer()) 
x0, y0, width, height = bbox.transformed(fig.transFigure.inverted()).bounds
# slightly increase the very tight bounds:
xpad = 0.013 * width
ypad = 0.03 * height
fig.add_artist(plt.Rectangle((x0-xpad-0.021, y0-ypad-0.03), width+2*xpad, height+2*ypad, edgecolor='0.1', linewidth=0.7, fill=False))


# define file paths for precipitation days
nc_dir = f'I:\\Emma\\FIROWatersheds\\Data\\Gridmet\\{watershed}'
path_totals = os.path.join(nc_dir,f'DailyMeans\\{watershed}DailymeansWINTERTOTAL_1980_2021.nc')
path_extremes = os.path.join(nc_dir,f'Extremes\\{watershed}DailymeanWINTEREXTREMES{percentile}_1980_2021.nc')

# import total precipitation data
totals = nc.Dataset(path_totals,mode='r')
print(totals.dimensions)
print(totals.variables)
precip = totals.variables['precipitation_amount'][:]
days = totals.variables['day'][:]
dates = [datetime(1979,1,1)+n*timedelta(days=1) for n in days]

# import extreme precipitation data 
extremes = nc.Dataset(path_extremes,mode='r')
ex_precip = extremes.variables['precipitation_amount'][:]

# plot data
markred = inset_ax.plot(dates,ex_precip,marker='.',linewidth=0,color='r',label=f'> {percentile}th Percentile',zorder=5)
lineblue = inset_ax.plot(dates,precip,linewidth=0.1,color='0.1')
markblue = inset_ax.plot(dates,precip,marker='.',linewidth=0,color='b',label=f'< {percentile}th Percentile')


if percentile == 95:
    extrlim = 63.17607498
elif percentile == 90:
    extrlim = 51.53282318
inset_ax.axhline(y=extrlim,color='red',linewidth=2,alpha=0.5)
inset_ax.set_xlabel('Year',fontweight='bold',fontsize=9)
inset_ax.set_ylabel('Precipitation (mm)',fontweight='bold',fontsize=9)
inset_ax.tick_params(direction='in',which='both',axis='y',labelsize=8)
inset_ax.tick_params(axis='x',labelsize=8)

#plt.title('Winter Season Daily Mean Precipitation',fontweight='bold')
# inset_ax.legend(loc='upper center',ncols=2,columnspacing=1.4,handletextpad=0.4,fontsize=9.7)
inset_ax.set_xbound(34E2,192E2)
inset_ax.set_ybound(0,215)

#%% CREATE CUSTOM LEGEND
import matplotlib.patches as mpatches

class BlueObject:
    pass
class RedObject:
    pass

class BlueObjectHandler:
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        rect = mpatches.Rectangle([x0, y0], width, height, facecolor='b',
                                   edgecolor='b',linewidth=1,alpha = 0.5,
                                   transform=handlebox.get_transform())
        circ = mpatches.Circle([x0 + 28, y0 + height/2],radius=2.5, facecolor='b',
                                   transform=handlebox.get_transform())
        handlebox.add_artist(rect)
        handlebox.add_artist(circ)
        return rect, circ
    
class RedObjectHandler:
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        rect = mpatches.Rectangle([x0, y0], width, height, facecolor='r',
                                   edgecolor='r',linewidth=1,alpha = 0.5,
                                   transform=handlebox.get_transform())
        circ = mpatches.Circle([x0 + 28, y0 + height/2],radius=2.5, facecolor='r',
                                   transform=handlebox.get_transform())
        handlebox.add_artist(rect)
        handlebox.add_artist(circ)
        return rect, circ


ax.legend([BlueObject(),RedObject()], [f'< {percentile}th Percentile',f'> {percentile}th Percentile'],
          handler_map={BlueObject: BlueObjectHandler(),RedObject: RedObjectHandler()},loc='lower right',
          handletextpad=1.5)

save_dir='I:\\Emma\\FIROWatersheds\\Figures\\NodeHistograms'
os.chdir(save_dir)
plt.savefig(f'{percentile}_{numpatterns}_PrecipDistributionandClim.png',dpi=300)
plt.show()