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
import shapefile 

#%% DEFINE WATERSHED SHAPEFILE
watershed = 'UpperYuba'
ws_directory = f'I:\\Emma\\FIROWatersheds\\Data\\WatershedShapefiles\\California\\{watershed}\\'

#%% 
# import extreme dates of upper yuba
extremedates = np.load('I:/Emma/FIROWatersheds/Data/90Percentile_ExtremeDays.npy',allow_pickle=True)

#%% OPEN AR CATALOG DATA
folder = 'I:\\CCAR-2022\\'
file = 'globalARcatalog_MERRA2_1980-2021_v3.0.nc'
path = os.path.join(folder,file)
ds = nc.Dataset(path, mode='r')
time = ds.variables['time'][:] # in seconds since 19800101, every 6 hours
gridlat = ds.variables['lat'][:]
gridlon = ds.variables['lon'][:] # in 0 to 360
# ds.close()

# convert time from hours to datetimes
dates = np.array([datetime.strftime(datetime(1900,1,1)+timedelta(minutes=time[n]),"%Y%m%d") for n in range(time.shape[0])])
    
#%%
for i in range(len(time)):
    # print(kid[i])
# for i in dated:
    latmin, latmax = (-90,90)
    lonmin, lonmax = (0,360)
    
    #reduce lat
    latlims = np.logical_and(gridlat > latmin, gridlat < latmax)
    latind = np.where(latlims)[0]
    gridlatreduced = gridlat[latind]
    #reduce lon
    lonlims = np.logical_and(gridlon > lonmin, gridlon < lonmax)
    lonind = np.where(lonlims)[0]
    gridlonreduced = gridlon[lonind]
    
    shape = np.squeeze(ds.variables['kidmap'][:,i,:,:,:])
    datareduced = shape[latind,:]
    datareduced = datareduced[:,lonind]
    
    shape2 = np.squeeze(ds.variables['shapemap'][:,i,:,:,:])
    datareduced2 = shape2[latind,:]
    datareduced2 = datareduced2[:,lonind]
    #plot
    fig = plt.figure(figsize=(7.2,5.1))
    # plt.title(f'{dates_red[i]},{arwatershedsum[i]}')
    #MAP DESIRED VARIABLE
    #convert lat and lon into a 2D array
    lon, lat = np.meshgrid(gridlonreduced,gridlatreduced)  
    #define area threshold for basemap
    area_thresh = 1E4
    #create equidistant cylindrical projection basemap
    map = Basemap(projection='cyl',llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,\
              urcrnrlon=lonmax,resolution='c',area_thresh=area_thresh)
    xi, yi = map(lon,lat)
    border_c = '0.4'
    border_w = 0.4
    #create map features
    map.drawcoastlines(color=border_c, linewidth=border_w,zorder=2)
    map.drawstates(color=border_c, linewidth=border_w,zorder=2)
    map.drawcountries(color=border_c, linewidth=border_w,zorder=2)
    
    parallels = np.arange(20.,71.,2.)
    meridians = np.arange(-160.,-109.,2.)
    # map.drawparallels(parallels, labels=[1,0,0,0])
    # map.drawmeridians(meridians, labels=[0,0,0,1])
    
    map.readshapefile(os.path.join(ws_directory,f'{watershed}'), watershed,zorder=6,color='k')
    plt.scatter(360-120.9,39.5,color='b',marker='*',linewidths=0.7,zorder=4)
    
    
    #create colormap of MERRA2 data
    colorm = map.pcolor(xi,yi,datareduced,shading='auto',cmap='tab20')
    # colorm = map.pcolor(xi,yi,datareduced2,shading='auto',cmap='tab20',alpha=0.5)
    #define border color and thickness
    plt.colorbar()
    plt.show()

