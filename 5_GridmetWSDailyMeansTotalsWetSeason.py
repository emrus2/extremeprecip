# -*- coding: utf-8 -*-
"""
Created on Wed May 10 15:48:01 2023

@author: emrus2

Calculated extreme precipitation days in Upper Yuba watershed from 
file of daily means. Computation is limited to winter months. 

Saved files:
    'UpperYubaDailymeansWINTERTOTAL{percentile}_1980_2021'
    'UpperYubaDailymeanWINTEREXTREMES{percentile}_1980_2021'

"""
#import required packages
import os
import xarray as xr

#define filepath and file of daily mean total precipitation
watershed = "UpperYuba"
folderpath = f'I:\\Emma\\FIROWatersheds\\Data\\Gridmet\\{watershed}\\DailyMeans'
filename = f'{watershed}dailymeansTOTAL_1980_2021.nc'
filepath = os.path.join(folderpath, filename)

#open filepath with xarray
NCdatatotal = xr.open_dataarray(filepath)
print(NCdatatotal)

#reduce total data to wet season data (October-March)
winter = [1, 2, 3, 10, 11, 12]
NCdatawinter =  NCdatatotal.where(NCdatatotal['day.month'].isin(winter), drop=True)
print(NCdatawinter)

percentile = 90
#calculate 90th percentile threshold
percent_90 = NCdatawinter.quantile(percentile/100)
print(percent_90)

#Calculate extreme days as those exceeding 95th percentile
extremes_winter = NCdatawinter.where(NCdatawinter>percent_90)
print(extremes_winter)
NCdatawinter.plot.line(marker='o',lw=0,markersize=2)
extremes_winter.plot.line(marker='o',lw=0,markersize=2)

#save values to nc files
# newfolder = f'I:\\Emma\\FIROWatersheds\\Data\\Gridmet\\{watershed}\\Extremes'
# NCdatawinter.to_netcdf(os.path.join(folderpath,f'{watershed}DailymeansWINTERTOTAL_1980_2021.nc'))
# extremes_winter.to_netcdf(os.path.join(newfolder,f'{watershed}DailymeanWINTEREXTREMES{percentile}_1980_2021.nc'))