# -*- coding: utf-8 -*-
# file_name.py
"""
Created on Fri December 9 15:27:00 2022

@author: emrus2

Creates file of GRIDMET extreme precipitation days

Saved files:
    # 'MERRA2_Z500_Yuba_Extremes_Daily_1980-2021_WINTERDIST.nc'
    # 'MERRA2_SLP_Yuba_Extremes_Daily_1980-2021_WINTERDIST.nc'
    # 'MERRA2_300W_Yuba_Extremes_Daily_1980-2021_WINTERDIST.nc'
    # 'MERRA2_850T_Yuba_Extremes_Daily_1980-2021_WINTERDIST.nc'
    # 'MERRA2_IVT_Yuba_Extremes_Daily_1980-2021_WINTERDIST.nc'
    # 'MERRA2_Z850_Yuba_Extremes_Daily_1980-2021_WINTERDIST.nc'

UPDATED 7/17/2023

"""
#%% IMPORT EXTREME DATES
#import desired modules
import numpy as np
import netCDF4 as nc
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import glob
import os
import xarray as xr
#%% IMPORT EXTREME DAYS
percentile = 90
extremedays = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_ExtremeDays.npy')

#%% ACCUMULATE GRIDMET DATA FOR EXTREME DAYS TOGETHER

# extremedaysprecip = np.zeros((1,585,1386))

# for year in range(1980,2022):
#     print(year)
#     gridmetfile =(f'I:\\GRIDMET\\pr\\pr_{year}.nc')
#     gridfile = nc.Dataset(gridmetfile,mode='r')
#     # print(gridfile)
#     # print(gridfile.dimensions)
#     # print(gridfile.variables)
#     precip = gridfile.variables['precipitation_amount'][:,:,:]
#     days = gridfile.variables['day'][:]
#     dates = [datetime.strftime(datetime(1900,1,1)+timedelta(days=float(days[n])),"%Y%m%d") for n in range(days.shape[0])]
#     gridfile.close()
    
#     for i,datestr in enumerate(dates):
#         if datestr in extremedays:
#             print(i,datestr)
#             precipday = (precip[i,:,:])
#             extremedaysprecip = np.concatenate((extremedaysprecip,np.expand_dims(precipday,axis=0)))
            
#%%
folderpath = 'I:\\GRIDMET\\pr'
files = glob.glob(os.path.join(folderpath,"*.nc"))
ds = xr.open_mfdataset(files, combine='nested', concat_dim='day')
print(ds)        
days = ds['day']
precip = ds['precipitation_amount']

#create datestring array for easy assignment
dates = [datetime.strftime(datetime(1979,1,1)+n*timedelta(days=1),"%Y%m%d") for n in range(days.shape[0])]

extremedaysprecip = np.zeros((1,585,1386))
for i,datestr in enumerate(dates):
    if datestr in extremedays:
        precipday = precip[i,:,:]
        preciparray = precipday.values
        extremedaysprecip = np.concatenate((extremedaysprecip,np.expand_dims(preciparray,axis=0)))
        
extremedaysprecip = extremedaysprecip[1:,:,:]

np.save('I:\\Emma\\FIROWatersheds\\Data\\GRIDMET_pr_Yuba_Extremes90_Daily_1980-2021_WINTERDIST',extremedaysprecip)
