# -*- coding: utf-8 -*-
# file_name.py
"""
Created on Fri December 9 15:27:00 2022

@author: emrus2

Compares MERRA-2 file names to extreme precipitation dates
Creates file of MERRA-2 extreme precipitation days

Saved files:
    'MERRA2_Z500_Yuba_Extremes_Daily_1980-2021_WINTERDIST.nc'
    'MERRA2_SLP_Yuba_Extremes_Daily_1980-2021_WINTERDIST.nc'
    'MERRA2_300W_Yuba_Extremes_Daily_1980-2021_WINTERDIST.nc'
    'MERRA2_850T_Yuba_Extremes_Daily_1980-2021_WINTERDIST.nc'
    'MERRA2_IVT_Yuba_Extremes_Daily_1980-2021_WINTERDIST.nc'
    'MERRA2_Z850_Yuba_Extremes_Daily_1980-2021_WINTERDIST.nc'

UPDATED 4/30/2023

#still need to get a few files foor Z850
"""
#%% IMPORT EXTREME DATES
#import desired modules
import os
import numpy as np
import xarray as xr
import glob
import netCDF4 as nc
#import numpy as np
from datetime import datetime, timedelta

#import extreme values from GRIDMET precip data
watershed = "UpperYuba"
percentile = 90
gridmetpath = f'I:\\Emma\\FIROWatersheds\\Data\\Gridmet\\{watershed}\\Extremes'
gridmetfile = f"{watershed}DailymeanWINTEREXTREMES{percentile}_1980_2021.nc" 
    # 7655 days of precipitation, with missing data for days not extreme
    # time in days since 1979-01-01
gridmet = os.path.join(gridmetpath,gridmetfile)
extremes = nc.Dataset(gridmet,mode='r')
print(extremes) #7655 in length
print(extremes.dimensions)
print(extremes.variables)
precipdays = extremes.variables['day'][:]
precip = extremes.variables['precipitation_amount'][:] #time x height x lat x lon
extremes.close()

#convert datetime to matching string of MERRA2 files
datestr = [datetime.strftime(datetime(1979,1,1)+timedelta(days=float(precipdays[n])),"%Y%m%d") for n in range(precip.shape[0])]

#create list of datetimes where precip is extreme
extremedates = []
precipvalues = []
#extremes are all where precip is not NaN in dataset
for n in range(precip.shape[0]):
    if precip[n]:
        extremedates.append(datestr[n])
        precipvalues.append(precip[n])
        
#np.save(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_ExtremeDays.npy',np.array(extremedates))
#np.save(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_ExtremeDaysPrecip.npy',np.array(precipvalues))

        #should be equal to 130 for 95th, 260 for 90th
#remove extreme days in summer (keep OCT-MAR)
#extremestr = [date for date in extremedates if int(date[5]) not in range(4,10)]
    # has already been done in the gridmet calculations ^^^
    #should be equal to 260

#%% DEFINE EXTREME MERRA2 FILES
#define files of interest
metvars = ['IVT']
metpath = {'Z500':'500_hPa_Geopotential_Height_3hourly','SLP':'Sea_level_pressure_3hourly','IVT':'IVT_daily', \
               '300W':'East_and_North_wind_components_at_300_hPa','850T':'Temperature_at_850_hPa_3hourly', \
                   'Z850':'850_hPa_Geopotential_Height_3hourly','850W':'East_and_North_wind_components_at_850_hPa'}

for metvar in metvars:
    folderpath = f'I:\\MERRA2\\Daily_and_Subdaily\\{metpath[metvar]}'
    os.chdir(folderpath)
    extremefiles = [file for file in os.listdir(folderpath) if any(date in file for date in extremedates)]
    # sort extremefiles in case it is out of order
    extremefiles.sort() # works for IVT

            
    #%% CREATE EXTREMES DATASET
    #create dataset of extreme files
    if metvar == 'IVT':
        ds = xr.open_mfdataset(extremefiles, combine='nested', concat_dim='time')
        ds = ds.drop_vars(['UFLXQL','VFLXQL'])
    else:
        ds = xr.open_mfdataset(extremefiles, combine='nested',concat_dim='time')
        #calculate daily average merra data
        ds = ds.groupby('time.date').mean() #this works great!!
        ds.squeeze() #should result in 162(260) dates
        ds['date'] = ds['date'].astype('datetime64')
        
    print(ds)
    
    #calculate composite mean for merra
    #total_means = np.squeeze(ds.mean(dim='time'))
    #print(total_means) #should be 1 time dimension
    
    #%% SAVE EXTREMES DATASET
    #save as file
    #define new files directory and names
    os.chdir(f'I:\\Emma\\FIROWatersheds\\Data\\DailyMERRA2\\{metvar}')
    savefile = f'MERRA2_{metvar}_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST_updated.nc'
    #savefile2 = f'{metvar}_extremes_composite.nc'
    ds.to_netcdf(savefile)
    #total_means.to_netcdf(savefile2)