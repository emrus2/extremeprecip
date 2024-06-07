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

# open file of average wet season precip
meanfile = 'I:\\Emma\\FIROWatersheds\\Data\\WinterClimatologies\\Gridmet_WetSeason_Climatology_1980-2021.nc'
wsmean = xr.open_dataarray(meanfile) # 260 days
print(wsmean)
#open extreme dataset
extremefile = 'I:\\Emma\\FIROWatersheds\\Data\\Gridmet\\GRIDMET_pr_Yuba_Extremes90_Daily_1980-2021_WINTERDIST.nc'
extremes = xr.open_dataarray(extremefile) # 260 days
print(extremes)

#divide data to get anomalies
anoms = extremes / wsmean

#define new files directory and names
os.chdir('I:\\Emma\\FIROWatersheds\\Data\\Gridmet')
newfile = 'GRIDMET_pr_Yuba_Extremes90_DailyAnomProp_1980-2021_WINTERDIST_wetseason.nc'
anoms.to_netcdf(newfile)
