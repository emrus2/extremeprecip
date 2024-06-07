# -*- coding: utf-8 -*-
# file_name.py
"""
Created on Fri December 9 15:27:00 2022

@author: emrus2

Description: This script calculates the wet season average precipitation
from GRIDMET


UPDATED 07/11/2023

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
files.pop(-1) #remove 2023
files.pop(-1) # remove 2022
#files will be 42 years, 15341 days

winterlist = [10,11,12,1,2,3]

#concatenate all files along time dimension and define as dataset
ds = xr.open_mfdataset(files,combine='nested',concat_dim="day")

# group by month
ds['months'] = ds['day'].dt.month
# reduce dataset to wet season only
ds_wetseason = ds.sel(day=ds['day.month'].isin(winterlist))
    # wet season should be reduced to 7655 days
print(ds_wetseason)
# calculate wet season mean
ds_ws_mean = np.squeeze(ds_wetseason.mean(dim='day'))

#define new files directory and names
os.chdir('I:\\Emma\\FIROWatersheds\\Data\\WinterClimatologies')
newfile = 'Gridmet_WetSeason_Climatology_1980-2021.nc'
ds_ws_mean.to_netcdf(newfile)

print("**************************************************")
print("********************Finished!!!*******************")
print("**************************************************")