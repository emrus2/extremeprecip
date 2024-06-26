# -*- coding: utf-8 -*-
# file_name.py
"""
Created on Fri December 9 15:27:00 2022

@author: emrus2

Creates file of GRIDMET extreme precipitation days, 260 days

Saved files:
    'GRIDMET_pr_Yuba_Extremes90_Daily_1980-2021_WINTERDIST.nc'

UPDATED 7/17/2023

"""
#%% IMPORT EXTREME DATES
#import desired modules
import numpy as np
#import netCDF4 as nc
#from datetime import datetime, timedelta
#import matplotlib.pyplot as plt
import glob
import os
import xarray as xr
import pandas as pd
#%% IMPORT EXTREME DAYS
percentile = 90
extremedays = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_ExtremeDays.npy')

#convert list of strings to datetime dataset
extremedays_dt = pd.to_datetime(extremedays)
#%%
#metvars = ['precip','mintemp','maxtemp','10mwinddir','10mwindspeed']
metvar = '10mwindspeed'
metabbr = {'precip':'pr','mintemp':'tmmn','maxtemp':'tmmx','10mwinddir':'th','10mwindspeed':'vs'}
metfolder = {'precip':'pr','mintemp':'temp','maxtemp':'temp','10mwinddir':'th_(10m_wind_direction)','10mwindspeed':'vs_(10m_wind_speed)'}


folderpath = f'I:\\GRIDMET\\{metfolder[metvar]}'    
files = glob.glob(os.path.join(folderpath,f"{metabbr[metvar]}*.nc")) #should be length 44
ds = xr.open_mfdataset(files, combine='nested', concat_dim='day') #16071 days
print(ds)      
ds_extreme = ds.sel(day = extremedays_dt) #should be 260 days

#%% SAVE EXTREMES DATASET
#save as file
#define new files directory and names
os.chdir('I:\\Emma\\FIROWatersheds\\Data\\Gridmet')
savefile = f'GRIDMET_{metabbr[metvar]}_Yuba_Extremes90_Daily_1980-2021_WINTERDIST.nc'
ds_extreme.to_netcdf(savefile)
