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
import seaborn as sns
import xarray as xr
#from cftime import num2date, date2num

watershed_name = "UpperYuba"
nc_dir = f'I:\\Emma\\FIROWatersheds\\Data\\Gridmet\\{watershed_name}'

percentile = 90

path_totals = os.path.join(nc_dir,f'DailyMeans\\{watershed_name}DailymeansWINTERTOTAL_1980_2021.nc')
path_extremes = os.path.join(nc_dir,f'Extremes\\{watershed_name}DailymeanWINTEREXTREMES{percentile}_1980_2021.nc')

#totals = xr.open_dataarray(path_totals)
totals = nc.Dataset(path_totals,mode='r')
NCdatawinter = xr.open_dataarray(path_totals)

print(totals.dimensions)
print(totals.variables)

precip = totals.variables['precipitation_amount'][:]
days = totals.variables['day'][:]
dates = [datetime(1979,1,1)+n*timedelta(days=1) for n in days]

#totals = xr.open_dataarray(path_totals)
extremes = nc.Dataset(path_extremes,mode='r')
ex_precip = extremes.variables['precipitation_amount'][:]

#%%
#line1 = plt.plot(dates,extremes,marker='.',linewidth=0.1,color='red',label='> 95th Percentile',zorder=5)
#line2 = plt.plot(dates,precip,marker='.',linewidth=0.1,color='blue',label ='< 95th Percentile')
fig, ax = plt.subplots(layout='constrained')

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
# ax.legend(loc='upper left')
ax.set_xbound(3348,19300)


ax2 = ax.twiny()
plot = sns.histplot(data=ex_precip, kde=True, bins=75, zorder=1, \
                  color = 'slategrey',stat='percent',alpha=0.7, \
                      ax=ax2)

extrlim = 51.53282318
ymax = 21
xmax = 208
# ax.axvline(x=extrlim,color='red',linewidth=2,alpha=0.5)
# plt.fill_betweenx(y=plot, x1=extrlim, x2=xmax,color='red')
for line in ax2.lines:
    x, y = line.get_xydata().T
    ax2.fill_between(x, 0, y, color='royalblue', where=x<extrlim,alpha=0.5,label = '<90th Percentile')
    ax2.fill_between(x, 0, y, color='red', where=x>extrlim,alpha=0.5,label='>90th Percentile')
    # ax.fill_betweenx(y, x, extrlim, color='red')

# #CUSTOMIZE SUBPLOT SPACING
ax2.set_ylabel('Percent of Days (%)',fontweight='bold')
ax2.set_xlabel('Precipitation (mm)',fontweight='bold')
# ax2.legend(loc='upper center',ncols=9,columnspacing=1.4,handletextpad=0.4,fontsize=9.7)
# ax2.set_ylim(0,ymax)
# ax2.set_xlim(-1,xmax)

save_dir=f'I:\\Emma\\FIROWatersheds\\Figures\\{watershed_name}'
os.chdir(save_dir)
plt.savefig(f'DailyMeanPrecip_{percentile}Extremes_nolegend2.png',bbox_inches='tight',pad_inches=0.1,dpi=300)
plt.show()
#%%
# plt.hist(totals,bins=50,color='blue')
# plt.hist(extremes,bins=50,color="red")
# plt.xlabel('Number of Days')
# plt.ylabel('Precipitation (mm)')
# plt.title('Histogram of Daily Mean Precipitation in Upper Yuba Watershed (1979-2022)')