# -*- coding: utf-8 -*-
# file_name.py
"""
Created on Fri December 9 15:27:00 2022

@author: emrus2

Description: This script calculates the mean and standard deviation climatologies
for the OCT-MAR  period for desired meteorological variable

Final result is two .nc files, one for mean and another for standard deviation of values

UPDATED 07/11/2023

STATUS: WORKS, TAKES A LONG TIME FOR CALCULATION AND FILE SAVING

Saved Files:
    'I:\\Emma\\FIROWatersheds\\Data\\WinterClimatologies'
        f'{metvar}_Climatology_1980-2021_OCT-MAR.nc'
        f'{metvar}_Climatology_Std_1980-2021_OCT-MAR.nc'
    
"""
#import desired modules
import os
import xarray as xr
import glob
import numpy as np

metvar = 'Z300'
filepath = {'Z300':'300_hPa_Geopotential_Height_3hourly','Z500':'500_hPa_Geopotential_Height_3hourly','SLP':'Sea_level_pressure_3hourly','850T':'Temperature_at_850_hPa_3hourly'}

print("**************************************************")
print("Calculating Climatology...")
print("**************************************************")

#define files of interest, in this case 1980-2021 years
folderpath = f'I:\\MERRA2\\Daily_and_Subdaily\\{filepath[metvar]}'
files = glob.glob(os.path.join(folderpath,"MERRA2*.inst3_3d_asm_Np.1????[0-3]*.SUB.nc")) + glob.glob(os.path.join(folderpath,'MERRA2*.inst3_3d_asm_Np.20[0-1]??[0-3]*.SUB.nc')) + glob.glob(os.path.join(folderpath,'MERRA2*.inst3_3d_asm_Np.202[0-1]?[0-3]*.SUB.nc'))
#files will be 7655

#concatenate all files along time dimension and define as dataset
ds = xr.open_mfdataset(files,combine='nested',concat_dim="time")
#calculate the mean of all files across time dimension
ds_mean = np.squeeze(ds.mean(dim='time'))
ds_std = np.squeeze(ds.std(dim='time'))

#define new files directory and names
os.chdir('I:\\Emma\\FIROWatersheds\\Data\\WinterClimatologies')
newfile_mean = f'{metvar}_Climatology_1980-2021_OCT-MAR.nc'
newfile_std = f'{metvar}_Climatology_Std_1980-2021_OCT-MAR.nc'
ds_mean.to_netcdf(newfile_mean)
print('Mean Done')
ds_std.to_netcdf(newfile_std)
print('Std Done')

print("**************************************************")
print("********************Finished!!!*******************")
print("**************************************************")