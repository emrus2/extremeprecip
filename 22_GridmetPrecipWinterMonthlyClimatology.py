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
#import desired modules
import os
import xarray as xr
import glob
import numpy as np



print("**************************************************")
print("Calculating Climatology...")
print("**************************************************")

#define files of interest, in this case 1980-2021 years
folderpath = 'I:\\GRIDMET\\pr\\'
files = glob.glob(os.path.join(folderpath,"*.nc"))
files.pop(0) #remove 1979
files.pop(-1) #remove 2022
#files will be 42 years, 15341 days

winterlist = [10,11,12,1,2,3]

#concatenate all files along time dimension and define as dataset
ds = xr.open_mfdataset(files,combine='nested',concat_dim="day")
ds_month = ds.groupby('day.month',squeeze=False).mean()
print(ds_month.precipitation_amount)

#define new files directory and names
os.chdir('I:\\Emma\\FIROWatersheds\\Data\\WinterClimatologies')
newfile = 'Gridmet_Monthly_Climatology_1980-2021.nc'
ds_month.to_netcdf(newfile)

print("**************************************************")
print("********************Finished!!!*******************")
print("**************************************************")