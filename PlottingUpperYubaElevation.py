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
latmin, latmax = (30.75,43.5)
lonmin, lonmax = (-125.25,-113.5)


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

#%% PLOT MAIN FIGURE

fig, ax = plt.subplots(layout='constrained',figsize=(7.2,5))

area_thresh = 1E4
latminsm,latmaxsm = (38,41)
lonminsm, lonmaxsm = (-123.5, -119.5)
map = Basemap(projection='cyl',llcrnrlat=latminsm,urcrnrlat=latmaxsm,llcrnrlon=lonminsm,urcrnrlon=lonmaxsm,resolution='i',area_thresh=area_thresh) #low resolution
map.bluemarble()

xi, yi = map(lon,lat)

#create a colormap of topography data
vmin,vmax = (-2500,3750)
# colorm = map.pcolor(xi,yi,etoponew,shading='auto',cmap='terrain',zorder=1,vmin=vmin,vmax=vmax)

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
# map.drawparallels(np.arange(37.,43.,1.), labels=[1,0,0,0], fontsize=10,color=border_c, linewidth=border_w)
# map.drawmeridians(np.arange(-125.,-119.,1.), labels=[0,0,0,1], fontsize=10,color=border_c, linewidth=border_w)

#PLOT WATERSHED SHAPE
map.readshapefile(os.path.join(ws_directory,f'{watershed}'), watershed)
# plt.scatter(-120.9,39.5,color='r',marker='*',s=300,linewidth=1,edgecolors='k',zorder=4)
#create plot title
#plt.title('24 OCT 2021',fontweight='bold')

lowlim,highlim = (0,3750)
#create colorbar
# cbar = map.colorbar(colorm, location='right', pad="5%",)
# cbar.set_label('Elevation (m)')

#INSERT SMALLER AXIS OF ZOOMED OUT VERSION


# this is an inset axes over the main axes
right_inset_ax = fig.add_axes([-0.01, .58, .4, .4], facecolor='k')
map2 = Basemap(projection='cyl',llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,urcrnrlon=lonmax,resolution='l',area_thresh=area_thresh,ax=right_inset_ax) #low resolution

colorm = map2.pcolor(xi,yi,etoponew,shading='auto',cmap='terrain',zorder=1,vmin=vmin,vmax=vmax,ax=right_inset_ax)
#define border color and thickness
border_c = '0.1'
border_w = 0.5

#map underlying map, with  parallel labels on the left, and meridian labels on the bottom
map2.drawlsmask(land_color='none',ocean_color='paleturquoise',zorder=2)
try:
    map2.drawcoastlines(color=border_c,linewidth=border_w,zorder=3)
except:
    pass
map2.drawstates(color=border_c, linewidth=border_w,zorder=3)
map2.drawcountries(color=border_c, linewidth=border_w,zorder=3)
# map.drawparallels(np.arange(37.,43.,1.), labels=[1,0,0,0], fontsize=10,color=border_c, linewidth=border_w)
# map.drawmeridians(np.arange(-125.,-119.,1.), labels=[0,0,0,1], fontsize=10,color=border_c, linewidth=border_w)

#PLOT WATERSHED SHAPE
# map.readshapefile(os.path.join(ws_directory,f'{watershed}'), watershed)
plt.scatter(-120.9,39.5,color='r',marker='*',s=100,linewidth=1,edgecolors='k',zorder=4)


#show map
save_dir = 'I:\\Emma\\FIROWatersheds\\Figures\\'
save_file = f'{watershed}_elevation.png'
plt.savefig(os.path.join(save_dir,save_file),dpi=300)
plt.show()
    