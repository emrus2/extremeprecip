# -*- coding: utf-8 -*-
"""
Created on Sat Jan 14 17:18:01 2023

@author: emrus2
"""
#MAP WATERSHED VARIABLES
"""
maps nc data in region of watershed
"""
#%% IMPORTING PACKAGES
import os
#import pandas as pd
import netCDF4 as nc
import numpy as np
#import xarray as xr
import matplotlib.pyplot as plt
#import matplotlib.dates as mdates
#import matplotlib.transforms as mtransforms
#from matplotlib.lines import Line2D
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
from mpl_toolkits.basemap import Basemap

#%% UPLOADING TOPOGRAPHY DATA -only complete on initial HYSPLIT mapping, takes a long time
#define ETOPO data filepath
etopofile = 'ETOPO1_Bed_g_gmt4.grd'
etopodir = 'I:\\Emma\\LaborDayWildfires\\Data\\Topo_Data\\ETOPO1_Bedrock'
etopopath = os.path.join(etopodir,etopofile)

#COLLECT VARIABLE DATA FROM ETOPO FILE
gridfile = nc.Dataset(etopopath)
gridlon = gridfile.variables['x'][:]
gridlat = gridfile.variables['y'][:]
etopo = gridfile.variables['z'][:] #lat x lon
gridfile.close()

#define resolution of topography data
interval = 1 #the higher the interval, the lower the resolution
skip1 = slice(None,None,interval)
skip = (slice(None, None, interval), slice(None, None, interval))
gridlatreduced = gridlat[skip1]
gridlonreduced = gridlon[skip1]
etoporeduced = etopo[skip]

#%%
import rasterio
etopofile = 'gt30w140n40.tif'
etopodir = 'C:\\Users\\emrus2\\Downloads'
etopopath = os.path.join(etopodir,etopofile)

file = rasterio.open(etopopath)
dataset = file.read()
print(dataset.shape)

plt.imshow(dataset[0], cmap='Spectral')
plt.show()
#%% DEFINE WATERSHED
watershed = 'UpperYuba'
ws_directory = f'I:\\Emma\\FIROWatersheds\\Data\\WatershedShapefiles\\California\\{watershed}\\'

#%% DEFINE LAT LON BOUNDS FOR PLOTTING
latmin,latmax = (38,41)
lonmin, lonmax = (-123.5, -119.5)

#reduce the latitude, longitude, and etopo to Western US
latlims = np.logical_and(gridlatreduced > latmin, gridlatreduced < latmax)
latind = np.where(latlims)[0]
gridlatnew = gridlatreduced[latind]

lonlims = np.logical_and(gridlonreduced > lonmin, gridlonreduced < lonmax)
lonind = np.where(lonlims)[0]
gridlonnew = gridlonreduced[lonind]

etoponew = etoporeduced[latind,:]
etoponew = etoponew[:,lonind]

#convert lat and lon into a 2D array
lon, lat = np.meshgrid(gridlonnew,gridlatnew)

#%% PLOT FIGURE

fig = plt.figure(figsize=(7.2,5))

area_thresh = 1E4
map = Basemap(projection='cyl',llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,urcrnrlon=lonmax,resolution='l',area_thresh=area_thresh) #low resolution
xi, yi = map(lon,lat)

#create a colormap of topography data
vmin,vmax = (-2500,3750)
colorm = map.pcolor(xi,yi,etoponew,shading='auto',cmap='terrain',zorder=1,vmin=vmin,vmax=vmax)

#define border color and thickness
border_c = '0.1'
border_w = 0.5

#map underlying map, with  parallel labels on the left, and meridian labels on the bottom
map.drawlsmask(land_color='none',ocean_color='paleturquoise',zorder=2)
try:
    map.drawcoastlines(color=border_c,linewidth=border_w,zorder=3)
except:
    pass
map.drawstates(color=border_c, linewidth=border_w,zorder=3)
map.drawcountries(color=border_c, linewidth=border_w,zorder=3)
map.drawparallels(np.arange(37.,43.,1.), labels=[1,0,0,0], fontsize=10,color=border_c, linewidth=border_w)
map.drawmeridians(np.arange(-125.,-119.,1.), labels=[0,0,0,1], fontsize=10,color=border_c, linewidth=border_w)

#PLOT WATERSHED SHAPE
map.readshapefile(os.path.join(ws_directory,f'{watershed}'), watershed)

#create plot title
#plt.title('24 OCT 2021',fontweight='bold')

lowlim,highlim = (0,3750)
#create colorbar
cbar = map.colorbar(colorm, location='right', pad="5%",)
cbar.set_label('Elevation (m)')

#show map
save_dir = f'I:\\Emma\\FIROWatersheds\\Figures\\{watershed}'
save_file = f'{watershed}_elevation.png'
plt.savefig(os.path.join(save_dir,save_file),dpi=300)
plt.show()
    