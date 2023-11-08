# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 12:24:51 2023

@author: emrus2
"""

""" 
We now have tehe following saved files:
    UpperYubadailymeanTOTAL.nc
    UpperYubadailymeanEXTREMES.nc
    UpperYubadailymeanEXTREMES_1980_2021.nc
    
Let's create some preliminary plots
"""

import os
#import xarray as xr
import netCDF4 as nc
import matplotlib.pyplot as plt
#import numpy as np
from datetime import datetime, timedelta
#from cftime import num2date, date2num

watershed_name = "UpperYuba"
nc_dir = f'I:\\Emma\\FIROWatersheds\\Data\\Gridmet\\{watershed_name}'

percentile = 90


path_totals = os.path.join(nc_dir,f'DailyMeans\\{watershed_name}DailymeansWINTERTOTAL_1980_2021.nc')
path_extremes = os.path.join(nc_dir,f'Extremes\\{watershed_name}DailymeanWINTEREXTREMES{percentile}_1980_2021.nc')

#totals = xr.open_dataarray(path_totals)
totals = nc.Dataset(path_totals,mode='r')

print(totals.dimensions)
print(totals.variables)

precip = totals.variables['precipitation_amount'][:]
days = totals.variables['day'][:]
dates = [datetime(1979,1,1)+n*timedelta(days=1) for n in days]

# print(totals)
# print(list(totals.data_vars))
# print(list(totals.coords))
# extreme_thresh = totals.quantile(0.95)
# print(extreme_thresh) #63.1760

#totals = xr.open_dataarray(path_totals)
extremes = nc.Dataset(path_extremes,mode='r')
ex_precip = extremes.variables['precipitation_amount'][:]

#%%
#line1 = plt.plot(dates,extremes,marker='.',linewidth=0.1,color='red',label='> 95th Percentile',zorder=5)
#line2 = plt.plot(dates,precip,marker='.',linewidth=0.1,color='blue',label ='< 95th Percentile')
fig, ax = plt.subplots(layout='constrained')

line3 = plt.plot(dates,precip,marker='.',linewidth=0,color='blue',label=f'< {percentile}th Percentile')
line1 = plt.plot(dates,ex_precip,marker='.',linewidth=0,color='red',label=f'> {percentile}th Percentile',zorder=5)
line2 = plt.plot(dates,precip,marker='.',linewidth=0.1,color='blue')


if percentile == 95:
    extrlim = 63.17607498
elif percentile == 90:
    extrlim = 51.53282318
ax.axhline(y=extrlim,color='red',linewidth=2,alpha=0.5)
ax.set_xlabel('Year',fontweight='bold')
ax.set_ylabel('Precipitation (mm)',fontweight='bold')
ax.tick_params(direction='in',which='both',axis='y')
#plt.title('Winter Season Daily Mean Precipitation',fontweight='bold')
ax.legend(loc='upper center',ncols=2,columnspacing=1.4,handletextpad=0.4,fontsize=9.7)
ax.set_xbound(3348,19300)

# save_dir=f'I:\\Emma\\FIROWatersheds\\Figures\\{watershed_name}'
# os.chdir(save_dir)
# plt.savefig(f'DailyMeanPrecip_{percentile}Extremes_nolegend.png',bbox_inches='tight',pad_inches=0.1,dpi=300)
plt.show()
#%%
# plt.hist(totals,bins=50,color='blue')
# plt.hist(extremes,bins=50,color="red")
# plt.xlabel('Number of Days')
# plt.ylabel('Precipitation (mm)')
# plt.title('Histogram of Daily Mean Precipitation in Upper Yuba Watershed (1979-2022)')