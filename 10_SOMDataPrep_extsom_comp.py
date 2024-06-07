# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 13:49:49 2023

@author: emrus2

Prepares MERRA2 Extreme precipitation data for use in SOM Toolbox in MATLAB
Saves list of Extreme Days clustered by 5-days

Saved Files:
    'I:\\Emma\\FIROWatersheds\\Data\\SOMs\\SomInput'
        f'{metvar}_{percentile}_data_{clusters}d.mat'
    'I:\\Emma\\FIROWatersheds\\Data\\'    
        f'ExtremeEvents_{clusters}d.npy'
        
UPDATED 6/27/2023
"""
#%% IMPORT MODULES
import netCDF4 as nc #if this doesn't work, try to reinstall in anaconda prompt using
    #conda install netcdf4
import numpy as np
import os
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
import scipy.io
from datetime import datetime, timedelta
#%% IMPORT MERRA2
os.chdir('I:\\Emma\\FIROWatersheds\\Data')
#define MERRA2 data location
percentile = 90
daysprior = 2
clusters = daysprior + 1
metvar = 'IVT'

folderpath = f'DailyMERRA2\\{metvar}'
filename = f'MERRA2_{metvar}_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST_{daysprior+1}d.nc'
filepath = os.path.join(folderpath,filename)

#%% IMPORT EXTREME DATA
# import list of extreme days
extremedates = np.load(f'90Percentile_ExtremeDays_{clusters}d.npy')
startdate = extremedates [0] # first day of merra2 data
#convert extremedates tp datetime
extremedates = [datetime.strptime(i,'%Y%m%d') for i in extremedates]

#%% CLUSTER DAYS TOGETHER
# cluster dates into events
extremeevents = [] # should be equal to number of extreme days (260)
for d in range(0,len(extremedates),clusters):
    # print(d)
    event = []
    for i in range(clusters):
        # print(i)
        event.append(extremedates[d+i])
    extremeevents.append(event)
#%%
#COLLECT VARIABLE DATA FROM MERRA2 FILE
#open the netcdf file in read mode
gridfile = nc.Dataset(filepath,mode='r')
print(gridfile.variables) # look at units of date
date = gridfile.variables['date'][:]
date = [datetime.strptime(startdate,'%Y%m%d') + timedelta(days=int(j)) for j in date] #3d - 632
gridlat = gridfile.variables['lat'][:]
gridlon = gridfile.variables['lon'][:]
Uvapor = gridfile.variables['UFLXQV'][:]
Vvapor = gridfile.variables['VFLXQV'][:]
merra = np.sqrt(Uvapor**2 + Vvapor**2)
gridfile.close()

# reduce merra to 3D array
merra = np.squeeze(merra)
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

#%% LATITUDINAL NORMALIZATION
#weight data based on area (square root of the cosine of latitude)
latweights = np.sqrt(np.cos((np.radians(gridlatreduced))))
for i in range(merrareduced.shape[1]): #lats
    merrareduced[:,i,:] *= latweights[i]
    
#%% FLATTEN MERRA INTO 2D ARRAY, THEN BY DAT
# flatten merra by combining lat,lon dimension
merraflat = merrareduced.reshape(len(merrareduced),-1)

flatlen = len(merraflat[1]) #should be lat x lon
# create array to store flattened data
alldata = np.zeros((len(extremeevents),flatlen*clusters))
# loop through extreme precipitation events
for j,event in enumerate(extremeevents):
    # create empty array to store event data
    merraevent = np.zeros((clusters,flatlen))
    for count,day in enumerate(event):
        # find index of merra data for the extreme date
        dayindx = date.index(day)
        # define merra data as that of index, append to storage array
        merraday = merraflat[dayindx,:] # each day
        merraevent[count,:] = merraday # whole event
    # horizontally concatenate merra data for each event
    data = np.hstack(merraevent)
    alldata[j,:] = data # should be numextremedays x clusters*lonsize*latsize

#%% SAVE DATA ARRAY TO OPEN IN MATLAB
mat_dir='I:\\Emma\\FIROWatersheds\\Data\\SOMs\\SomInput'
os.chdir(mat_dir)
datadict = {'lat': gridlatreduced, 'lon': gridlonreduced, 'merra': alldata}
# scipy.io.savemat(f'{metvar}_{percentile}_data_{clusters}d.mat', mdict = datadict)

# save list of extreme events
data_dir = 'I:\\Emma\\FIROWatersheds\\Data\\'
np.save(os.path.join(data_dir,f'ExtremeEvents_{clusters}d.npy'),extremeevents)
