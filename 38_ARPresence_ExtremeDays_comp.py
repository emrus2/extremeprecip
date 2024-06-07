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
time = ds.variables['time'][:] 
print(ds.variables['time']) # in minutes since 19000101, every 6 hours (4 timesteps per day)
gridlat = ds.variables['lat'][:]
gridlon = ds.variables['lon'][:] # in 0 to 360

# convert time from hours to datetimes
dates = np.array([datetime.strftime(datetime(1900,1,1)+timedelta(minutes=time[n]),"%Y%m%d") for n in range(time.shape[0])])
#%% REDUCE CATALOG TO EXTREME DATES
# find extreme dates in AR catalog dataset
#reduce dates
datelims = np.isin(dates,extremedates) # find dates that are present in extreme dates list
dateind = np.where(datelims)[0] # find those indices
datereduced = dates[dateind] # reduce dates to indices of those presnt in extreme dates list
dates_red = sorted(list(set(datereduced))) # sort dates back to chronological order

#%% IMPORT SHAPE DATA
shape = np.squeeze(ds.variables['kidmap'][:,dateind,:,:,:]) #should be numextremedays x 4 (1040)
    # kidmap shows the shape of AR with values filled with the k-id if each AR
    # reduced to the extreme day indices, still 6-hourly data
    
#%% SUM AR PRESENCE ACROSS EACH DAY
# create an empty array
days_summed = np.zeros((int((len(shape)/4)),len(gridlat),len(gridlon))) #empty array
counter = 0
for i in range(0,len(shape),4): # loop through all AR timesteps
    day = shape[i:i+4,:,:] # snip to the day (4 timesteps)
    daysum = np.sum(day,axis=0) # calculate sum of AR frequency in the day
    days_summed[counter,:,:] = daysum # add to empty array
    counter += 1

#%% REDUCE AR CATALOG DATA TO WATERSHED REGION
# import watershed data
watershed_name= 'UpperYuba'
ws_directory = f'I:\\Emma\\FIROWatersheds\\Data\\WatershedShapefiles\\California\\{watershed_name}\\'
ws_filename = (f'{watershed_name}.shp')
ws_path = os.path.join(ws_directory,ws_filename)
shp = shapefile.Reader(ws_path, crs="epsg:4326")
feature = shp.shapes()

# Extract array of lat/lon coordinates
coords = np.squeeze(np.array([s.points for s in shp.shapes()]))
shp_lon = np.array(coords)[:,0] + 360 # convert to match AR catalago 0-360 lat coords
shp_lat = np.array(coords)[:,1]

# set gridcells of watershed region
latminws, latmaxws = (38.5,40.5)
lonminws, lonmaxws = (360-122.5,360-119.9)

# reduce lat
latlimsws = np.logical_and(gridlat > latminws, gridlat < latmaxws)
latindws = np.where(latlimsws)[0]
gridlatws = gridlat[latindws]
# reduce lon
lonlimsws = np.logical_and(gridlon > lonminws, gridlon < lonmaxws)
lonindws = np.where(lonlimsws)[0]
gridlonws = gridlon[lonindws]
# reduce AR data
arwatershed = days_summed[:,latindws,:]
arwatershed = arwatershed[:,:,lonindws]

#%% IDENTIFY AR PRESENCE IN WATERSHED REGION
# add up AR frequency across watershed gridcells
arwatershedsum = np.sum(arwatershed,axis=(1,2))
# if gridcells > 0, ar is present, which means True
arpresence = arwatershedsum > 0

# combine ar presence data with dates data
dates_arpresence = np.array(list(zip(dates_red,arpresence)))

# save as numpy array
savefolder = 'I:/Emma/FIROWatersheds/Data/'
np.save(os.path.join(savefolder,'ARPresence_90_Percentile_Days.npy'),dates_arpresence)

#%% PLOT AR FREQUENCY TO CONFIRM ALGORITHM
for i in range(len(days_summed)):
# for i in dated:
    latmin, latmax = (25,55)
    lonmin, lonmax = (360-150,360-110)
    
    #reduce lat
    latlims = np.logical_and(gridlat > latmin, gridlat < latmax)
    latind = np.where(latlims)[0]
    gridlatreduced = gridlat[latind]
    #reduce lon
    lonlims = np.logical_and(gridlon > lonmin, gridlon < lonmax)
    lonind = np.where(lonlims)[0]
    gridlonreduced = gridlon[lonind]
    
    datareduced = days_summed[i,:,:]
    datareduced = datareduced[latind,:]
    datareduced = datareduced[:,lonind]

    #plot
    fig = plt.figure(figsize=(7.2,5.1))
    plt.title(f'{dates_red[i]},{arwatershedsum[i]}')
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
    map.drawparallels(parallels, labels=[1,0,0,0])
    map.drawmeridians(meridians, labels=[0,0,0,1])
    
    map.readshapefile(os.path.join(ws_directory,f'{watershed}'), watershed,zorder=6,color='k')
    plt.scatter(360-120.9,39.5,color='b',marker='*',linewidths=0.7,zorder=4)
    
    
    #create colormap of MERRA2 data
    colorm = map.pcolor(xi,yi,datareduced,shading='auto',cmap='gnuplot2_r')
    #define border color and thickness
    plt.colorbar()
    plt.show()