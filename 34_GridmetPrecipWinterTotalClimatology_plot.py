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
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.transforms as mtransforms
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
from mpl_toolkits.basemap import Basemap #installed using 
    #conda install -c anaconda basemap
import scipy.io

#%% DEFINE WATERSHED SHAPEFILE
watershed = 'UpperYuba'
ws_directory = f'I:\\Emma\\FIROWatersheds\\Data\\WatershedShapefiles\\California\\{watershed}\\'

#%% IMPORT SOM AND GRIDMET DATA
# change directory and import SOM data from .mat file
numpatterns = 9
percentile = 90
clusters = 5


# define lat, lon region of data for plotting
latmin, latmax = (37.5,41.5)
lonmin, lonmax = (-124.5,-119.5)
latmin, latmax = (37,42.25)
lonmin, lonmax = (-125.25,-118.75)

#%%
#IMPORT NETCDF DATA
#define NC location
# load in monthly gridmet data
os.chdir('I:\\Emma\\FIROWatersheds\\Data\\WinterClimatologies')
filepath = 'Gridmet_Monthly_Climatology_1980-2021.nc'
gridfile = nc.Dataset(filepath,mode='r')
print(gridfile)
gridlat = gridfile.variables['lat'][:]
gridlon = gridfile.variables['lon'][:]
months = gridfile.variables['month'][:]
precip = np.squeeze(gridfile.variables['precipitation_amount'][:])
# days = gridfile.variables['day'][:]
gridfile.close()
#dates = [datetime.strftime(datetime(1900,1,1)+timedelta(days=days[n]),"%Y%m%d") for n in range(days.shape[0])]

#Reduce to winter months
wintermonths = [1, 2, 3, 10, 11, 12]
#Reduce precip to winter months
winterprecip = np.array([pr for i,pr in enumerate(precip) if i+1 in wintermonths])
#Calculate winter climatological mean
meanwinterprecip = np.mean(winterprecip,axis=0)

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
merrareduced = meanwinterprecip[latind,:]
merrareduced = merrareduced[:,lonind]

#convert lat and lon into a 2D array
lon, lat = np.meshgrid(lonreduced,latreduced)
#define area threshold for basemap
area_thresh = 1E4
    
#%% CLUSTER ASSIGNED PATTERNS AND CALCULATE AVERAGE

# # create array to store 12 SOM composites
# som_composites = np.zeros((numpatterns,len(latreduced),len(lonreduced)))

# # loop through all 12 som patterns
# for som in range(numpatterns):
#     # create array to store assigned days data
#     som_merra = np.zeros((1,len(latreduced),len(lonreduced)))
#     # loop through all days
#     for day,arr in enumerate(merrareduced):
#         print(np.amax(arr))
#         # add data to som_merra if day is assigned to node
#         if assignment[day] == float(som + 1):
#             som_merra = np.ma.concatenate((som_merra,np.expand_dims(arr,axis=0)))
#             print(np.amax(som_merra))
#     # remove initial row of zeros
#     som_merra = som_merra[1:,:,:]
#     # calculate the mean of assigned days
#     som_mean = np.squeeze(np.mean(som_merra,axis=0))
#     # append to array of composites
#     som_composites[som] = som_mean
#%% DETERMINE MAX AND MIN VALIUES
zmax = 0
zmin = 1E8

for i, arr in enumerate(merrareduced):
    
    #determine zmax and zmin for all days
    highlim = np.nanmax(arr)
    lowlim = np.nanmin(arr)
    print(lowlim,highlim)
    if highlim > zmax:
        zmax = highlim
    if lowlim < zmin:
        zmin = lowlim

print(f'Lowest Value:{zmin} \nHighest Value:{zmax}')

#%% PLOT NODES from MATLAB

#create subplot for mapping multiple timesteps
fig = plt.figure(figsize=(7.2,6.9))
# fig.suptitle('Precipitation Composites',fontsize=13,fontweight="bold",y=0.9875)

lowlim = 0
highlim = 19
colors = ['white',"yellow",'greenyellow',"limegreen","lightseagreen",'royalblue','mediumblue','#7400E0','#B800E0','#E0ADB1'] #,'mediumorchid','#A600E0','pink']
colormap = LinearSegmentedColormap.from_list("mycmap", colors)
# titles = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

# for i, arr in enumerate(merrareduced):
#MAP DESIRED VARIABLE
#convert lat and lon into a 2D array
#define area threshold for basemap
area_thresh = 1E4
#create equidistant cylindrical projection basemap
map = Basemap(projection='cyl',llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,\
          urcrnrlon=lonmax,resolution='l',area_thresh=area_thresh)
xi, yi = map(lon,lat)
# ax = fig.add_subplot(4,3,i+1)
sublabel_loc = mtransforms.ScaledTranslation(4/72, -100/72, fig.dpi_scale_trans)
# ax.text(0.0, 1.0, titles[i], transform=ax.transAxes + sublabel_loc,
#     fontsize=9, fontweight='bold', verticalalignment='top', 
#     bbox=dict(facecolor='1', edgecolor='none', pad=1.5),zorder=3)
#create colormap of MERRA2 data
colorm = map.pcolor(xi,yi,merrareduced,shading='auto',cmap=colormap,vmin=lowlim,vmax=highlim,zorder=1)

#define border color and thickness
border_c = '0.4'
border_w = 0.4
#create map features
map.drawcoastlines(color=border_c, linewidth=border_w)
map.drawstates(color=border_c, linewidth=border_w)
map.drawcountries(color=border_c, linewidth=border_w)
gridlinefont = 8.5
parallels = np.arange(38.,42.,1.)
meridians = np.arange(-124.,-119.,2.)
if i in range(0,7,3):
    map.drawparallels(parallels, labels=[1,0,0,0], fontsize=gridlinefont,color=border_c,linewidth=border_w)
    map.drawmeridians(meridians,color=border_c,linewidth=border_w)
elif i == 10 or i == 11:
    map.drawparallels(parallels, color=border_c,linewidth=border_w)
    map.drawmeridians(meridians, labels=[0,0,0,1], fontsize=gridlinefont,color=border_c,linewidth=border_w)
elif i == 9:
    map.drawparallels(parallels, labels=[1,0,0,0], fontsize=gridlinefont,color=border_c,linewidth=border_w)
    map.drawmeridians(meridians, labels=[0,0,0,1], fontsize=gridlinefont,color=border_c,linewidth=border_w)
else:
    map.drawparallels(parallels, color=border_c,linewidth=border_w)
    map.drawmeridians(meridians,color=border_c,linewidth=border_w)
#define contour color and thickness
contour_c = '0.1'
contour_w = 0.7
#create contour map
#contourm = map.contour(xi,yi,arr,colors=contour_c,linewidths=contour_w,levels=np.arange(contourstart[metvar],highlims[metvar]+1,contourint[metvar]),zorder=2)
#plt.clabel(contourm,levels=contourm.levels[::2],fontsize=6,inline_spacing=1,colors='k',zorder=2,manual=False)
    
#add yuba shape
map.readshapefile(os.path.join(ws_directory,f'{watershed}'), watershed)
#plt.scatter(-120.9,39.5,color='tomato',edgecolors='r',marker='*',linewidths=0.8,zorder=4)
    
#CUSTOMIZE SUBPLOT SPACING
# fig.subplots_adjust(left=0.05,right=0.9,bottom=0.026, top=0.985,hspace=0.05, wspace=0.05) #bottom colorbar
#fig.add_axis([left,bottom, width,height])
cbar_ax = fig.add_axes([0.91,0.05,0.025,0.9]) #bottom colorbar
cbar = fig.colorbar(colorm, cax=cbar_ax,ticks=np.arange(0,highlim+1,4),orientation='vertical')
cbar.ax.tick_params(labelsize=8)
cbar.set_label('mm',fontsize=8.5,labelpad=0.5,fontweight='bold')

    
#SHOW MAP
save_dir='I:\\Emma\\FIROWatersheds\\Figures\\SOMs\\Composites'
os.chdir(save_dir)
# plt.savefig('GRIMET_pr_MonthlyAverages.png',dpi=300)
plt.show()