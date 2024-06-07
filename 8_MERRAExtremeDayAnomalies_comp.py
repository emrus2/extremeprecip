# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 13:49:49 2023

@author: emrus2

Calculated MERRA-2 anomalies for meteorological variables of interest

saves:
    f'I:\\Emma\\FIROWatersheds\\Data\\DailyMERRA2\\{metvar}'
    'MERRA2_{metvar}Anom_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST.npy

UPDATED 6/27/2023
"""
#%% IMPORT MODULES
import netCDF4 as nc #if this doesn't work, try to reinstall in anaconda prompt using
    #conda install netcdf4
import numpy as np
import os
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
#%% IMPORT MERRA2

metvar = 'Z300'
merravar = {'Z300':'H','Z500':'H','SLP':'SLP','850T':'T'}
percentile = 90
clusters = 5

folderpath = f'I:\\Emma\\FIROWatersheds\\Data\\DailyMERRA2\\{metvar}'
filename = f'MERRA2_{metvar}_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST_{clusters}d.nc'
filepath = os.path.join(folderpath,filename)

#COLLECT VARIABLE DATA FROM MERRA2 FILE
#open the netcdf file in read mode
gridfile = nc.Dataset(filepath,mode='r')
gridlat = gridfile.variables['lat'][:]
gridlon = gridfile.variables['lon'][:]
merra = np.squeeze(gridfile.variables[merravar[metvar]][:])
dates = gridfile.variables['date'][:]
gridfile.close()
  
#calculate Z500 anomaly
gridfile_mean = nc.Dataset(f'I:\\Emma\\FIROWatersheds\\Data\\WinterClimatologies\\{metvar}_Climatology_1980-2021_OCT-MAR.nc',mode='r')
climmean = gridfile_mean.variables[merravar[metvar]][:]
gridfile_mean.close()
gridfile_std = nc.Dataset(f'I:\\Emma\\FIROWatersheds\\Data\\WinterClimatologies\\{metvar}_Climatology_Std_1980-2021_OCT-MAR.nc',mode='r')
climstd = gridfile_mean.variables[merravar[metvar]][:]
gridfile_std.close()
merraanom = merra - climmean
merra = merraanom / climstd

merra = np.array(merra)
os.chdir(f'I:\\Emma\\FIROWatersheds\\Data\\DailyMERRA2\\{metvar}')
np.save(f'MERRA2_{metvar}Anom_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST_{clusters}d',merra)