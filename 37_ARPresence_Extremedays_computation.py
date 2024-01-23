# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 15:35:52 2024

@author: emrus2
"""
#%% IMPORT PACKAGES
import os
import netCDF4 as nc
import numpy as np
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
from mpl_toolkits.basemap import Basemap #installed using 
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

#%% DEFINE WATERSHED SHAPEFILE
watershed = 'UpperYuba'
ws_directory = f'I:\\Emma\\FIROWatersheds\\Data\\WatershedShapefiles\\California\\{watershed}\\'

#%% OPEN AR CATALOG DATA
folder = 'I:\\ARDetection\\CMIP6_Historical\\'
file = 'cmip6_IPSL-CM6A-LR_historical_r1i1p1f1_gr.ar_tag.Guan&Waliser_v2.6hr.195001010600-201501010000.nc4'
path = os.path.join(folder,file)

ds = nc.Dataset(path, mode='r')
time = ds.variables['time'][:]
arindicator = ds.variables['ar_binary_tag'][:]
gridlat = ds.variables['lat'][:]
gridlon = ds.variables['lon'][:]
ds.close()

# convert lon to -180 to 180
gridlon = np.array([arr-360 for i,arr in enumerate(gridlon) if arr > 180] + [arr for i,arr in enumerate(gridlon) if arr <= 180])

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

#%% convert time from hours to datetimes
dates = [datetime.strftime(datetime(1950,1,1)+timedelta(hours=time[n]),"%Y%m%d") for n in range(time.shape[0])]
#%% REDUCE CATALOG TO EXTREME DAYS
# import extreme dates of upper yuba
extremedates = np.load('I:/Emma/FIROWatersheds/Data/90Percentile_ExtremeDays.npy',allow_pickle=True)

# remove duplicate dates in date array
dates_catal = sorted(list(set(dates)))

# create list of extreme dates in catalog
catextremes = [string for string in dates_catal if string in extremedates]

# sum binary tags for whole days
days_summed = np.zeros((int((len(dates)/4)),len(gridlat),len(gridlon)))
counter = 0
for i in range(0,len(arindicator),4):
    print(i)
    daysum = np.sum(arindicator[i:i+3,:,:],axis=0)
    days_summed[counter,:,:] = daysum
    counter += 1

# reduce summed catalog to extreme days
counter = 0
ar_extremedays = np.zeros((len(catextremes),len(gridlat),len(gridlon)))
for i,arr in enumerate(days_summed):
    if dates_catal[i] in extremedates:
        ar_extremedays[counter,:,:] = arr
        counter += 1

# reduce ar catalog to lat and lon bounds        
artestreduced = ar_extremedays[:,latind,:]
artestreduced = artestreduced[:,:,lonind]

#%% IDENTIFY AR ACTIVITY
# reduce data to only watershed gridcells
arwatershed = ar_extremedays[:,102,:]
arwatershed = arwatershed[:,21:23]

# add up all gridcells for each day
arwatershedsum = np.sum(arwatershed,axis=1)
# if gridcells > 0, ar is present, which means True
arpresence = arwatershedsum > 0

# combine ar presence data with dates data
dates_arpresence = np.array(list(zip(catextremes,arpresence)))

# save as numpy array
savefolder = 'I:/Emma/FIROWatersheds/Data/'
np.save(os.path.join(savefolder,'ARPresence_90_Percentile_Days.npy'),dates_arpresence)

#%%
# # plot
# for i,j in enumerate(catextremes):
#     artest = artestreduced[i,:,:]
#     fig = plt.figure(figsize=(7.2,5.1))
#     plt.title(j)
#     #MAP DESIRED VARIABLE
#     #convert lat and lon into a 2D array
#     lon, lat = np.meshgrid(gridlonreduced,gridlatreduced)  
#     #define area threshold for basemap
#     area_thresh = 1E4
#     #create equidistant cylindrical projection basemap
#     map = Basemap(projection='cyl',llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,\
#               urcrnrlon=lonmax,resolution='c',area_thresh=area_thresh)
#     xi, yi = map(lon,lat)
#     border_c = '0.4'
#     border_w = 0.4
#     #create map features
#     map.drawcoastlines(color=border_c, linewidth=border_w,zorder=2)
#     map.drawstates(color=border_c, linewidth=border_w,zorder=2)
#     map.drawcountries(color=border_c, linewidth=border_w,zorder=2)
    
#     map.readshapefile(os.path.join(ws_directory,f'{watershed}'), watershed,zorder=6,color='w')
#     plt.scatter(120.9,39.5,color='w',marker='*',linewidths=0.7,zorder=4)


#     #create colormap of MERRA2 data
#     colorm = map.pcolor(xi,yi,artest,shading='auto')
#     #define border color and thickness
    
#     plt.show()