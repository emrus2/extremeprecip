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
import xarray as xr
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
import scipy.io
from datetime import datetime, timedelta
#%% IMPORT MERRA2
os.chdir('I:\\Emma\\FIROWatersheds\\Data')
#define MERRA2 data location
#metvars = ['SLP', '300W', 'Z500Anom','SLPAnom']
percentile = 90
daysprior = 4
clusters = daysprior + 1
metvar = 'IVT'
#for metvar in metvars:
merravar = {'Z500':'H','SLP':'SLP','850T':'T','Z850':'H','Z500Anom':'H'}

folderpath = f'DailyMERRA2\\{metvar}'
filename = f'MERRA2_{metvar}_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST_{daysprior+1}d.nc'
filepath = os.path.join(folderpath,filename)
ds = xr.open_mfdataset(filepath, combine='nested',concat_dim='time')

#%%
#COLLECT VARIABLE DATA FROM MERRA2 FILE
#open the netcdf file in read mode
gridfile = nc.Dataset(filepath,mode='r')
date = gridfile.variables['date'][:]
date = [datetime.strptime('19800105','%Y%m%d') + timedelta(days=int(j)) for j in date]
gridlat = gridfile.variables['lat'][:]
gridlon = gridfile.variables['lon'][:]
Uvapor = gridfile.variables['UFLXQV'][:]
Vvapor = gridfile.variables['VFLXQV'][:]
merra = np.sqrt(Uvapor**2 + Vvapor**2)
gridfile.close()

# reduce merra to 3D array
merra = np.squeeze(merra)

#%%
# import list of extreme days
extremedates = np.load(f'90Percentile_ExtremeDays_{clusters}d.npy')
#convert extremedates tp datetime
extremedates = [datetime.strptime(i,'%Y%m%d') for i in extremedates]

# cluster dates into events
extremeevents = []
for d in range(0,len(extremedates),clusters):
    print(d)
    event = []
    for i in range(clusters):
        print(i)
        event.append(extremedates[d+i])
    extremeevents.append(event)
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

#%%
#weight data based on area (square root of the cosine of latitude)
latweights = np.sqrt(np.cos((np.radians(gridlatreduced))))
for i in range(merrareduced.shape[1]): #lats
    merrareduced[:,i,:] *= latweights[i]
    
#%% FLATTEN MERRA INTO 2D ARRAY, THEN BY DAT
merraflat = merrareduced.reshape(len(merrareduced),-1)

flatlen = len(merraflat[1])
# loop through extreme precipitation events
alldata = np.zeros((len(extremeevents),flatlen*clusters))
for j,event in enumerate(extremeevents):
    # create empty array to store event data
    merraevent = np.zeros((clusters,flatlen))
    for count,day in enumerate(event):
        # find associated index
        # if day in date: # can remove once the 19800107 is done
            # find index of merra data for the extreme date
            dayindx = date.index(day)
            # define merra data as that of index, append to storage array
            merraday = merraflat[dayindx,:]
            merraevent[count,:] = merraday
    # horizontally concatenate merra data for each event
    # merraevent = merraevent.reshape(len(merraevent),-1)
    # data = np.concatenate(merraevent,axis=0)
    data = np.hstack(merraevent)
    alldata[j,:] = data

#%% HORIZONTALLY CONCATENATE EXTREME EVENTS, THEN FLATTEN
# loop through extreme precipitation events
# alldata = np.zeros((len(extremeevents),len(gridlatreduced),len(gridlonreduced)*4))
# for j,event in enumerate(extremeevents):
#     # create empty array to store event data
#     merraevent = np.zeros((clusters,len(gridlatreduced),len(gridlonreduced)))
#     for count,day in enumerate(event):
#         # find associated index
#         if day in date: # can remove once the 19800107 is done
#             # find index of merra data for the extreme date
#             dayindx = date.index(day)
#             # define merra data as that of index, append to storage array
#             merraday = merrareduced[dayindx,:,:]
#             merraevent[count,:,:] = merraday
#     # horizontally concatenate merra data for each event
#     # merraevent = merraevent.reshape(len(merraevent),-1)
#     # data = np.concatenate(merraevent,axis=0)
#     data = np.hstack(merraevent)
#     alldata[j,:,:] = data

# #flatten data from 3D to 2D (260x10197)
# alldata = alldata.reshape(len(alldata),-1)

#%% SAVE DATA ARRAY TO OPEN IN MATLAB
mat_dir='I:\\Emma\\FIROWatersheds\\Data\\SOMs\\SomInput'
os.chdir(mat_dir)
datadict = {'lat': gridlatreduced, 'lon': gridlonreduced, 'merra': alldata}
# scipy.io.savemat(f'{metvar}_{percentile}_data_{clusters}d.mat', mdict = datadict)

# save list of extreme events
np.save(f'ExtremeEvents_{clusters}d.npy',extremeevents)
