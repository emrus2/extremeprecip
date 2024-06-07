# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 14:50:11 2023

@author: emrus2
"""
#%% IMPORT PACKAGES
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import os

#%% IMPORT GRIDMET DAILY AVERAGE DATA
ds = xr.open_dataarray('I:/Emma/FIROWatersheds/Data/Gridmet/UpperYuba/DailyMeans/UpperYubadailymeansTOTAL_1980_2021.nc')
print(ds)
# calculate average monthly precipitation for dataset
ds_month = ds.groupby('day.month',squeeze=False).mean()
print(ds_month)

#%% calculate total precip
allprecip = sum(ds_month)
print(allprecip)
#%% calculate wet-season precip
winter = [1,2,3,10,11,12]
ds_winter = ds_month.where(ds_month.month.isin(winter),drop=True)
wintersum = sum(ds_winter)
print(wintersum)
print(wintersum/allprecip)
#%% IMPORT EXTREME PRECIPITATION DAYS DATA
extremedays = np.load('I:/Emma/FIROWatersheds/Data/90Percentile_ExtremeDays.npy')

# define months for each day
extrememonths = [date[4:6] for date in extremedays]
# convert month strings to integers
monthsint =  [int(y) for y in extrememonths]
# count the number of month occurrences
monthcounts = np.array([[month,monthsint.count(month)] for month in set(monthsint)])
#%%

# define month labeling
months = ("J", "F", "M",'A','M','J','J','A','S',"O", "N", "D")

# set plot parameters
barwidth = 0.7
offset = 0.2
axcolor = 'cornflowerblue'
ax2color = 'blue'


fig, ax = plt.subplots(layout='constrained')
ax.bar(ds_month.month,ds_month,color=axcolor,width=barwidth)
# ax.plot(ds_month.month,ds_month,color='cornflowerblue')

# define axis 1 parameters
ax.set_zorder(1)
ax.set_facecolor('none')
# labels
ax.set_ylabel('Average Precipitation (mm)',fontweight='bold',color=axcolor)
ax.set_xlabel('Month',fontweight='bold')
# x and y limits
ymax = float(max(ds_month)) + 1
xmin, xmax = (0.4,12.8)
ax.set_ylim(0,ymax)
ax.set_xlim(xmin,xmax)
# ticks
ax.set_yticks(range(0,21,5),range(0,21,5),color=axcolor)
ax.set_xticks(ds_month.month,months)
ax.tick_params(direction='in',which='both',axis='y',color=axcolor)
# axis color
ax.spines['left'].set_color(axcolor)
ax.spines['right'].set_color(ax2color)
# twin x axis for second axis
ax2 = ax.twinx()
ax2.bar(monthcounts[:,0]+offset,monthcounts[:,1],color=ax2color,width=barwidth)
# ax2.plot(monthcounts[:,0]+0.1,monthcounts[:,1],color='blue',marker='o',linewidth=0)

# define axis 2 parameters
ax2.set_zorder(0)
ax2.set_facecolor('none')
# labels
ax2.set_ylabel('Number of Extreme Days',fontweight='bold',color=ax2color)
# ticks
ax2.tick_params(direction='in',which='both',axis='y',color=ax2color)
ax2.set_yticks(range(0,51,10),range(0,51,10),color=ax2color)

# save file
os.chdir('I:/Emma/FIROWatersheds/Figures/PrecipDistributions')
plt.savefig('AverageMonthlyPrecip_ExtremeDays.png',dpi=300)
plt.show()

