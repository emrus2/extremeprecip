# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 13:49:49 2023

@author: emrus2

Here we will bin the days assigned to each IVT SOM together and then plot the
composite precipitation pattern from GRIDMET

For a 9-node SOM

UPDATED 6/12/2023
"""
#%% IMPORT MODULES
import netCDF4 as nc #if this doesn't work, try to reinstall in anaconda prompt using
    #conda install netcdf4
import numpy as np
import os
import matplotlib.pyplot as plt
#from datetime import datetime, timedelta
#from matplotlib.colors import ListedColormap
#import matplotlib.cm as cm
import matplotlib.transforms as mtransforms
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
from mpl_toolkits.basemap import Basemap #installed using 
    #conda install -c anaconda basemap
#import scipy.io

#%% DEFINE WATERSHED SHAPEFILE
watershed = 'UpperYuba'
ws_directory = f'I:\\Emma\\FIROWatersheds\\Data\\WatershedShapefiles\\California\\{watershed}\\'

#%% IMPORT SOM AND GRIDMET DATA
# change directory and import SOM data from .mat file
numpatterns = 9
percentile = 90
pats_assignments = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_ExtremeDaysand{numpatterns}NodeAssign.npy')
assignment = pats_assignments[:,0]
assignment = [float(x) for x in assignment]

# define lat, lon region of data for plotting
latmin, latmax = (37.5,41.5)
lonmin, lonmax = (-124.5,-119.5)

filepath = 'I:\\GRIDMET\\th_(10m_wind_direction)\\th_1980.nc'
gridfile = nc.Dataset(filepath,mode='r')
print(gridfile)
griddata = np.squeeze(gridfile.variables['wind_from_direction'][:]) #degrees clockwise from north
gridlat = gridfile.variables['lat'][:]
gridlon = gridfile.variables['lon'][:]
days = gridfile.variables['day'][:]
gridfile.close()
filepath2 = 'I:\\GRIDMET\\vs_(10m_wind_speed)\\vs_1980.nc'
gridfile2 = nc.Dataset(filepath2,mode='r')
print(gridfile2)
print(gridfile2.variables)
griddata2 = np.squeeze(gridfile2.variables['wind_speed'][:])
gridfile2.close()

#%%

#REDUCE VARIABLES TO DESIRED AREA
#reduce lat
latlims = np.logical_and(gridlat > latmin, gridlat < latmax)
latind = np.where(latlims)[0]
latreduced = gridlat[latind]
#reduce lon
lonlims = np.logical_and(gridlon > lonmin, gridlon < lonmax)
lonind = np.where(lonlims)[0]
lonreduced = gridlon[lonind]
#reduce pressure
merrareduced = griddata[0,latind,:]
merrareduced = merrareduced[:,lonind]

merra2reduced = griddata2[0,latind,:]
merra2reduced = merra2reduced[:,lonind]

#convert lat and lon into a 2D array
lon, lat = np.meshgrid(lonreduced,latreduced)
#define area threshold for basemap
area_thresh = 1E4

#%% CALCULATE U AND V WIND COMPONENTS

merra_radians = np.radians(270 - merrareduced)
Uwind = merra2reduced*np.cos(merra_radians)
Vwind = merra2reduced*np.sin(merra_radians)
    
#%% PLOT NODES from MATLAB

#create subplot for mapping multiple timesteps
fig = plt.figure(figsize=(7.5,5.5))
#fig.suptitle(f'{plottitles[metvar]} Composites',fontsize=13,fontweight="bold",y=0.9875)

precipsm = merra2reduced

#MAP DESIRED VARIABLE
#create equidistant cylindrical projection basemap
map = Basemap(projection='cyl',llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,\
              urcrnrlon=lonmax,resolution='l',area_thresh=area_thresh)
xi, yi = map(lon,lat)
#define z-axis bounds for colormap
#lowlim = 0
#highlim = 20
#create colormap of MERRA2 data
colorm = map.pcolor(xi,yi,precipsm,shading='auto',cmap='gist_ncar_r')
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
map.drawparallels(np.arange(38.,42.,1.), labels=[1,0,0,0], fontsize=10,color=border_c, linewidth=border_w)
map.drawmeridians(np.arange(-124.,-119.,2.), labels=[0,0,0,1], fontsize=10,color=border_c, linewidth=border_w)
#define contour color and thickness
#contour_c = '0.1'
#contour_w = 0.7
#create contour map
interval = 5
size = 1500
skip = (slice(None, None, interval), slice(None, None, interval))
#vectorm = map.quiver(xi2[skip],yi2[skip],U_arrs[skip],V_arrs[skip],color='darkgreen')
vectorm = map.quiver(xi[skip],yi[skip],Uwind[skip],Vwind[skip],pivot='mid',color='k')
#contourm = map.contour(xi,yi,pressurenew)
map.readshapefile(os.path.join(ws_directory,f'{watershed}'), watershed)

#create plot title
#plt.title(dates[i],fontweight='bold')

#create colorbar
cbar = map.colorbar(colorm, location='right', pad="5%")
cbar.set_label('mm')

#show map
#save_dir = f'I:\\Emma\\FIROWatersheds\\Figures\\{watershed}'
#save_file = f'{watershed}_precip_HIGHDAY_{year}.png'
#save_file = f'precip_HIGHDAY_{year}.png'
#plt.savefig(os.path.join(save_dir,save_file),dpi=300)
plt.show()
