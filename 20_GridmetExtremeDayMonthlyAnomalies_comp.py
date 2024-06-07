# -*- coding: utf-8 -*-
# file_name.py
"""
Created on Fri December 9 15:27:00 2022

@author: emrus2

Description: This script calculates the mean geopotential height for the OCT-MAR
climatological period for geopotential height

Final result is two .nc files, one for mean and another for standard deviation of values

UPDATED 07/11/2023

STATUS: WORKS, TAKES A LONG TIME FOR CALCULATION AND FILE SAVING
"""

import os
import xarray as xr
import glob
import numpy as np


#define files of interest, in this case 1980-2021 years
folderpath = 'I:\\GRIDMET\\pr\\'
files = glob.glob(os.path.join(folderpath,"*.nc"))
files.pop(0) #remove 1979
files.pop(-1) #remove 2022
#files will be 42 years, 15341 days


#concatenate all files along time dimension and define as dataset
ds = xr.open_mfdataset(files,combine='nested',concat_dim="day")

#open extreme dataset
extremefile = 'I:\\Emma\\FIROWatersheds\\Data\\Gridmet\\GRIDMET_pr_Yuba_Extremes90_Daily_1980-2021_WINTERDIST.nc'
extremes = xr.open_dataset(extremefile)

data1 = ds.groupby('day.month').mean()
smalldata = np.array(data1.precipitation_amount[0])
data2 = extremes.groupby('day.month')

anoms = extremes.groupby('day.month') / ds.groupby('day.month').mean()


#define new files directory and names
os.chdir('I:\\Emma\\FIROWatersheds\\Data\\Gridmet')
newfile = 'GRIDMET_pr_Yuba_Extremes90_DailyAnomProp_1980-2021_WINTERDIST.nc'
anoms.to_netcdf(newfile)