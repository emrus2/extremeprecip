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

fig, ax = plt.subplots(layout='constrained',figsize=(4.8,5))

area_thresh = 1E4
latminsm,latmaxsm = (38,41)
lonminsm, lonmaxsm = (-123.5, -119.5)
map = Basemap(projection='cyl',llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,urcrnrlon=lonmax,resolution='l',area_thresh=area_thresh) #low resolution
# map.bluemarble()
xi, yi = map(lon,lat)

#create a colormap of topography data
vmin,vmax = (-1000,4015)

colorm = map.pcolor(xi,yi,etoponew,shading='auto',cmap='gray',zorder=1,vmin=vmin,vmax=vmax)
# etoponew[etoponew <= -15] = 1E5
# colorm = map.pcolor(xi,yi,etoponew,shading='auto',cmap='gray',zorder=2,vmin=vmin,vmax=vmax)
#define border color and thickness
border_c = 'k'
border_w = 3
#map underlying map, with  parallel labels on the left, and meridian labels on the bottom
water = map.drawlsmask(land_color='none',ocean_color='white',zorder=2)
try:
    map.drawcoastlines(color='0.3',linewidth=border_w+3,zorder=3)
    map.drawcoastlines(color=border_c,linewidth=border_w,zorder=3)
except:
    pass
map.drawstates(color=border_c, linewidth=border_w,zorder=3)
map.drawcountries(color=border_c, linewidth=border_w,zorder=3)

#PLOT WATERSHED SHAPE
# map.readshapefile(os.path.join(ws_directory,f'{watershed}'), watershed,color='r',linewidth=2)
plt.scatter(-120.9,39.5,color='r',marker='*',s=600,linewidth=1,edgecolors='k',zorder=4)


#show map
save_dir = 'I:\\Emma\\FIROWatersheds\\Figures\\'
save_file = f'{watershed}_elevation.png'
plt.savefig(os.path.join(save_dir,save_file),dpi=300,bbox_inches='tight',pad_inches=0.05)
plt.show()
    