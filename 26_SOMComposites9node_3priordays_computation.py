# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 13:49:49 2023

@author: emrus2

Here we will bin the days assigned to each IVT SOM together and then plot the
composite patterns of other variables, like SLP, Z500Anom, etc. 

For a 9-node SOM

UPDATED 6/12/2023
"""
#%% IMPORT MODULES
import numpy as np
import os
import xarray as xr
from datetime import datetime, timedelta
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
#%% IMPORT SOM DATA
# change directory and import SOM assignment and date array
mat_dir='I:\\Emma\\FIROWatersheds\\Data'
os.chdir(mat_dir)
extremeassign = np.load('90Percentile_ExtremeDaysand9NodeAssign.npy')

# define file paths for meteorological variables
metpath = {'Z500':'500_hPa_Geopotential_Height_3hourly','SLP':'Sea_level_pressure_3hourly','IVT':'IVT_daily', \
           '300W':'East_and_North_wind_components_at_300_hPa','850T':'Temperature_at_850_hPa_3hourly', \
               'Z850':'850_hPa_Geopotential_Height_3hourly','850W':'East_and_North_wind_components_at_850_hPa'}

#%% IDENTIFY DATESTRINGS OF INTEREST
# define variable of interest
metvar = '850W'
# loop through som patterns 1-9
#for sompattern in range(1,10):
for sompattern in range(1,10):
    # loop through previous days 0-3
    for dayprev in range(0,4):
    #for dayprev in range(0,4):
        # reduce extreme days to those assigned to node
        extremenodedays = [str(j) for i,j in extremeassign if float(i) == float(sompattern)]
        # create list of days prior datestrings
        dayspriordt = [datetime.strptime(j,'%Y%m%d') - timedelta(days=dayprev) for j in extremenodedays]
        daysprior = [datetime.strftime(j,'%Y%m%d') for j in dayspriordt]
        
        #%% IMPORT MERRA-2 DATA
        folderpath = f'I:\\MERRA2\\Daily_and_Subdaily\\{metpath[metvar]}'
        os.chdir(folderpath)
        # identify files of desired node-assigned days
        files = [file for file in os.listdir(folderpath) if any(date in file for date in daysprior)]
        
        #%% CREATE DATASET FROM FILES
        if metvar == 'IVT':
            ds = xr.open_mfdataset(files, combine='nested', concat_dim='time')
            print(ds)
            ds = ds.drop_vars(['UFLXQL','VFLXQL'])
            print(ds.variables)
            # calculate day average U,V components
            means = np.squeeze(ds.mean(dim='time'))
            print(means) #should be 1 time dimension
            # calculate average IVT magnitude
            means = means.assign(IVTmag = np.sqrt(means.UFLXQV**2 + means.VFLXQV**2))
            means = means.drop_vars(['UFLXQV','VFLXQV'])
            print(means)            
        elif 'W' in metvar:
            ds = xr.open_mfdataset(files, combine='nested', concat_dim='time')
            print(ds.variables)
            # calculate day average U,V components
            means = np.squeeze(ds.mean(dim='time'))
            print(means) #should be 1 time dimension
            # calculate average IVT magnitude
            means = means.assign(Windmag = np.sqrt(means.U**2 + means.V**2))
            means = means.drop_vars(['U','V'])
            print(means)            
        else:
            ds = xr.open_mfdataset(files, combine='nested',concat_dim='time')
            #calculate daily average merra data
            ds = ds.groupby('time.date').mean() #this works great!!
            ds.squeeze() #should result in 162(260) dates
            ds['date'] = ds['date'].astype('datetime64')
            # calculate day composites
            means = np.squeeze(ds.mean(dim='date'))
            print(means) #should be 1 time dimension
            
        #%% SAVE PREVIOUS DAY COMPOSITE AS DATASET
        os.chdir(f'I:\\Emma\\FIROWatersheds\\Data\\SOMPreviousDaysComposites\\SOM{sompattern}')
        savefile = f'Node{sompattern}_{dayprev}dayprior_composite{metvar}.nc'
        means.to_netcdf(savefile)
