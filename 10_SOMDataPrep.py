# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 13:49:49 2023

@author: emrus2

Prepares MERRA2 Extreme precipitation data for use in SOM Toolbox in MATLAB

saves {metvar}_{percentile}_data.mat for matlab use in SOM Toolbox

UPDATED 6/27/2023
"""
#%% IMPORT MODULES
import netCDF4 as nc #if this doesn't work, try to reinstall in anaconda prompt using
    #conda install netcdf4
import numpy as np
import os
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
import scipy.io
#%% IMPORT MERRA2

#define MERRA2 data location
#metvars = ['SLP', '300W', 'Z500Anom','SLPAnom']
percentile = 90
metvar = 'IVT'
#for metvar in metvars:
merravar = {'Z500':'H','SLP':'SLP','850T':'T','Z850':'H','Z500Anom':'H'}

if metvar == 'Z500Anom':
    folderpath = 'I:\\Emma\\FIROWatersheds\\Data\\DailyMERRA2\\Z500'
    filename = f'MERRA2_Z500_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST.nc'
elif metvar == 'SLPAnom':
    folderpath = 'I:\\Emma\\FIROWatersheds\\Data\\DailyMERRA2\\SLP'
    filename = f'MERRA2_SLP_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST.nc'
else:
    folderpath = f'I:\\Emma\\FIROWatersheds\\Data\\DailyMERRA2\\{metvar}'
    filename = f'MERRA2_{metvar}_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST.nc'
filepath = os.path.join(folderpath,filename)

#COLLECT VARIABLE DATA FROM MERRA2 FILE
#open the netcdf file in read mode
gridfile = nc.Dataset(filepath,mode='r')
gridlat = gridfile.variables['lat'][:]
gridlon = gridfile.variables['lon'][:]
if metvar == 'IVT':
    Uvapor = gridfile.variables['UFLXQV'][:]
    Vvapor = gridfile.variables['VFLXQV'][:]
    merra = np.sqrt(Uvapor**2 + Vvapor**2)
elif metvar == '300W':
    Uwind = gridfile.variables['U'][:]
    Vwind = gridfile.variables['V'][:]
    merra = np.sqrt(Uwind**2 + Vwind**2)
elif metvar == 'Z500Anom':
    merra = np.load(os.path.join(folderpath,f'MERRA2_Z500Anom_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST.npy'))
elif metvar == 'SLPAnom':
    merra = np.load(os.path.join(folderpath,f'MERRA2_SLPAnom_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST.npy'))
else:
    merra = gridfile.variables[merravar[metvar]][:]
gridfile.close()

merra = np.squeeze(merra)

# convert to hPa (mb)
#if metvar == 'SLP':
#    merra = merra/100

#%% REDUCE TRAINING DATA TO REGION OF INTEREST
#reduce data to desired restrictions
#define lat lon restrictions
latmin, latmax = (15.5,65.5)
lonmin, lonmax = (-170.25,-105.75)

#reduce lat
latlims = np.logical_and(gridlat > latmin, gridlat < latmax)
latind = np.where(latlims)[0]
gridlatreduced = gridlat[latind]
#reduce lon
lonlims = np.logical_and(gridlon > lonmin, gridlon < lonmax)
lonind = np.where(lonlims)[0]
gridlonreduced = gridlon[lonind]
#reduce pressure
merrareduced = merra[:,latind,:]
merrareduced = merrareduced[:,:,lonind]

#%% PREPARE TRAINING DATA
#define data as numpy array (for minisom package)
data = np.asarray(merrareduced) #(130x99x103)

#normalize by temporal standard deviation
#can wait to do this step, seems only important for multivariate SOMs
''' Loikith 2017 uses temporal std, while Aragon 2020
and Taylor 2023 use the spatial std to normalize '''
# for index, day in enumerate(data):
#     print(np.mean(day),np.std(day))
#     day -= np.mean(day)
#     day /= np.std(day)
        
#weight data based on area (square root of the cosine of latitude)
latweights = np.sqrt(np.cos((np.radians(gridlatreduced))))
for i in range(data.shape[1]): #lats
    data[:,i,:] *= latweights[i]
    
#flatten data from 3D to 2D (260x10197)
data = data.reshape(len(data),-1)

#%% SAVE DATA ARRAY TO OPEN IN MATLAB
mat_dir='I:\\Emma\\FIROWatersheds\\Data\\SOMs\\SomInput'
os.chdir(mat_dir)
datadict = {'lat': gridlatreduced, 'lon': gridlonreduced, 'merra': data}
scipy.io.savemat(f'{metvar}_{percentile}_data.mat', mdict = datadict)