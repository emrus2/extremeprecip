# -*- coding: utf-8 -*-
# file_name.py
"""
Created on Fri December 9 15:27:00 2022

@author: emrus2

Compares MERRA-2 file names to extreme precipitation dates
Creates file of MERRA-2 extreme precipitation days and 
    the three days prior to each extreme day

Saved files:
    'MERRA2_IVT_Yuba_Extremes_Daily_1980-2021_WINTERDIST_5d.nc'


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

#%% ADD EXTREME DATES 3 DAYS PRIOR

# define function to return day previous
def dateprior(date,dayprev):
    dayspriordt = datetime.strptime(date,'%Y%m%d') - timedelta(days=dayprev)
    daysprior = datetime.strftime(dayspriordt,'%Y%m%d')
    return daysprior

extremedatesprior = []

daysprior = 4
# iterate through each string in list
for date in extremedates:
    for i in range(daysprior,0,-1):
        dayprior = dateprior(date,i)
        extremedatesprior.append(dayprior)
    extremedatesprior.append(date)
    
# remove any repeating values
#extremedatesreduced = sorted(list(set(extremedatesprior)))
    # should be 801 days
#%% DEFINE EXTREME MERRA2 FILES
#define files of interest
metvars = ['Z500','850T','Z850','850W']
metpath = {'Z500':'500_hPa_Geopotential_Height_3hourly','SLP':'Sea_level_pressure_3hourly','IVT':'IVT_daily', \
               '300W':'East_and_North_wind_components_at_300_hPa','850T':'Temperature_at_850_hPa_3hourly', \
                   'Z850':'850_hPa_Geopotential_Height_3hourly','850W':'East_and_North_wind_components_at_850_hPa'}

for metvar in metvars:
    folderpath = f'I:\\MERRA2\\Daily_and_Subdaily\\{metpath[metvar]}'
    os.chdir(folderpath)
    extremefiles = [file for file in os.listdir(folderpath) if any(date in file for date in extremedatesprior)]
    # sort extremefiles in case it is out of order
    #print(extremefiles)
    extremefiles.sort() # works for IVT
    extremefiles = [string for string in extremefiles if '.nc.1' not in string]
        # should be 801 days
            
    #%% CREATE EXTREMES DATASET
    ''' can I concatenate every three days? (horizontally) '''
    #create dataset of extreme files
    if metvar == 'IVT':
        ds = xr.open_mfdataset(extremefiles, combine='nested', concat_dim='time')
        ds = ds.drop_vars(['UFLXQL','VFLXQL'])        
        # for some reason 19800106 has hourly data, so calculate average
        
    else:
        ds = xr.open_mfdataset(extremefiles, combine='nested',concat_dim='time')
    #calculate daily average merra data
    ds = ds.groupby('time.date').mean() #this works great!!
    ds.squeeze() #should result in 162(260) dates
    ds['date'] = ds['date'].astype('datetime64')
        
    print(ds)
        # time dimension should be 801 days
    #calculate composite mean for merra
    #total_means = np.squeeze(ds.mean(dim='time'))
    #print(total_means) #should be 1 time dimension
    
    #%% SAVE EXTREMES DATASET
    #save as file
    #define new files directory and names
    os.chdir(f'I:\\Emma\\FIROWatersheds\\Data\\DailyMERRA2\\{metvar}')
    savefile = f'MERRA2_{metvar}_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST_{daysprior+1}d.nc'
    #savefile2 = f'{metvar}_extremes_composite.nc'
    ds.to_netcdf(savefile)
    #total_means.to_netcdf(savefile2)
    
    #%% SAVE LIST OF SUCCEEDING EXTREME DAYS
    # np.save(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_ExtremeDays_{daysprior+1}d.npy',np.array(extremedatesprior))

#%% CALCULATE AVERAGE OF FOLDER
# import os
# import xarray as xr
# metvar = 'IVT'
# daysprior = 2
# percentile = 90
# os.chdir(f'I:\\Emma\\FIROWatersheds\\Data\\DailyMERRA2\\{metvar}')
# savefile = f'MERRA2_{metvar}_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST_{daysprior+1}d.nc'
# ds = xr.open_mfdataset(savefile, combine='nested',concat_dim='time')
# print(ds.date)
