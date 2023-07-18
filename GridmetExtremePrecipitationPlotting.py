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
from datetime import datetime, timedelta
#from matplotlib.colors import ListedColormap
#import matplotlib.cm as cm
#import matplotlib.transforms as mtransforms
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
#%%
#IMPORT NETCDF DATA
#define NC location
filepath = ('I:\\Emma\\FIROWatersheds\\Data\\GRIDMET_pr_Yuba_Extremes90_Daily_1980-2021_WINTERDIST.nc')
#open the netcdf file in read mode
gridfile = nc.Dataset(filepath,mode='r')
print(gridfile)
gridlat = gridfile.variables['lat'][:]
gridlon = gridfile.variables['lon'][:]
precip = np.squeeze(gridfile.variables['precipitation_amount'][:])
days = gridfile.variables['day'][:]
dates = [datetime.strftime(datetime(1900,1,1)+timedelta(days=days[n]),"%Y%m%d") for n in range(days.shape[0])]

gridfile.close()

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
merrareduced = precip[:,latind,:]
merrareduced = merrareduced[:,:,lonind]

#convert lat and lon into a 2D array
lon, lat = np.meshgrid(lonreduced,latreduced)
#define area threshold for basemap
area_thresh = 1E4
#define map bounds
#latmin, latmax = (min(latreduced),max(latreduced))
#lonmin, lonmax = (min(lonreduced),max(lonreduced))
    
#%%
fig = plt.figure(figsize=(7,5))

for i in range(len(merrareduced)):
    precipsm = merrareduced[i,:,:]
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
    map.drawparallels(np.arange(38.,42.,1.), labels=[1,0,0,0], fontsize=10,color=border_c, linewidth=border_w)
    map.drawmeridians(np.arange(-124.,-119.,2.), labels=[0,0,0,1], fontsize=10,color=border_c, linewidth=border_w)
    #define contour color and thickness
    #contour_c = '0.1'
    #contour_w = 0.7
    #create contour map
    #contourm = map.contour(xi,yi,pressurenew)
    map.readshapefile(os.path.join(ws_directory,f'{watershed}'), watershed)
    
    #create plot title
    plt.title(dates[i],fontweight='bold')
    
    #create colorbar
    cbar = map.colorbar(colorm, location='right', pad="5%",ticks=np.arange(lowlim,highlim+1,100))
    cbar.set_label('mm')
    
    #show map
    #save_dir = f'I:\\Emma\\FIROWatersheds\\Figures\\{watershed}'
    #save_file = f'{watershed}_precip_HIGHDAY_{year}.png'
    #save_file = f'precip_HIGHDAY_{year}.png'
    #plt.savefig(os.path.join(save_dir,save_file),dpi=300)
    plt.show()