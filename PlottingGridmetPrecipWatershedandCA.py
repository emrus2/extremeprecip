# -*- coding: utf-8 -*-
"""
Created on Sat Jan 14 17:18:01 2023

@author: emrus2
"""
#MAP WATERSHED VARIABLES
"""
maps nc data in region of watershed
"""

#IMPORT MODULES
import numpy as np
import netCDF4 as nc
import os
import xarray as xr
import matplotlib.pyplot as plt
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
from mpl_toolkits.basemap import Basemap

#DEFINE WATERSHED
watershed = 'UpperYuba'
ws_directory = f'I:\\Emma\\FIROWatersheds\\Data\\WatershedShapefiles\\California\\{watershed}\\'

#IMPORT NETCDF DATA
year = '1980'
#define NC location
#nc_file = f'GRIDMET_PRECIP_{watershed}_{year}.nc'
nc_file = f'pr_{year}.nc'
#nc_file = f'GRIDMET_PRECIP_{year}.nc'
#nc_dir = f'I:\\Emma\\FIROWatersheds\\Data\\Gridmet\\{watershed}'
nc_dir = 'I:\\GRIDMET\\pr'
filepath = os.path.join(nc_dir,nc_file)
#open the netcdf file in read mode
gridfile = nc.Dataset(filepath,mode='r')
print(gridfile)
gridlat = gridfile.variables['lat'][:]
gridlon = gridfile.variables['lon'][:]
gridday = gridfile.variables['day'][:]
precip = gridfile.variables['precipitation_amount'][:,:,:] #time x height x lat x lon

#REDUCE VARIABLES TO DESIRED AREA
#convert height to a 3D array
precipsm = np.squeeze(precip)
#define lat lon restrictions
latmin = 38
latmax = 42
lonmin = -124
lonmax = -119
#reduce lat
latlims = np.logical_and(gridlat > latmin, gridlat < latmax)
latind = np.where(latlims)[0]
latreduced = gridlat[latind]
#reduce lon
lonlims = np.logical_and(gridlon > lonmin, gridlon < lonmax)
lonind = np.where(lonlims)[0]
lonreduced = gridlon[lonind]
#reduce pressure
precipreduced = precip[:,latind,:]
precipreduced = precipreduced[:,:,lonind]

#convert lat and lon into a 2D array
lon, lat = np.meshgrid(lonreduced,latreduced)
#define area threshold for basemap
area_thresh = 1E4
#define map bounds
latmin, latmax = (min(latreduced),max(latreduced))
lonmin, lonmax = (min(lonreduced),max(lonreduced))

#%%

fig = plt.figure(figsize=(7,5))

i = 337
precipsm = precipreduced[i,:,:]
print(np.amax(precipsm))

#MAP DESIRED VARIABLE
#create equidistant cylindrical projection basemap
map = Basemap(projection='cyl',llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,\
              urcrnrlon=lonmax,resolution='l',area_thresh=area_thresh)
xi, yi = map(lon,lat)
#define z-axis bounds for colormap
lowlim = 0
highlim = 500
#create colormap of MERRA2 data
colorm = map.pcolor(xi,yi,precipsm,shading='auto',cmap='gist_ncar_r',vmin=lowlim,vmax=highlim)
#define border color and thickness
border_c = '0.4'
border_w = 0.4
#create map features
try:
    map.drawcoastlines(color=border_c,linewidth=border_w,zorder=3)
except:
    pass
map.drawstates(color=border_c, linewidth=border_w)
map.drawcountries(color=border_c, linewidth=border_w)
map.drawparallels(np.arange(37.,43.,1.), labels=[1,0,0,0], fontsize=10,color=border_c, linewidth=border_w)
map.drawmeridians(np.arange(-125.,-119.,1.), labels=[0,0,0,1], fontsize=10,color=border_c, linewidth=border_w)
#define contour color and thickness
#contour_c = '0.1'
#contour_w = 0.7
#create contour map
#contourm = map.contour(xi,yi,pressurenew)
map.readshapefile(os.path.join(ws_directory,f'{watershed}'), watershed)

#create plot title
plt.title('24 OCT 2021',fontweight='bold')

#create colorbar
cbar = map.colorbar(colorm, location='right', pad="5%",ticks=np.arange(lowlim,highlim+1,100))
cbar.set_label('mm')

#show map
save_dir = f'I:\\Emma\\FIROWatersheds\\Figures\\{watershed}'
save_file = f'{watershed}_precip_HIGHDAY_{year}.png'
#save_file = f'precip_HIGHDAY_{year}.png'
#plt.savefig(os.path.join(save_dir,save_file),dpi=300)
plt.show()
    