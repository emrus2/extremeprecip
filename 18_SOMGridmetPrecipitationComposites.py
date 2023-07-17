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
#from matplotlib.colors import ListedColormap
#import matplotlib.cm as cm
import matplotlib.transforms as mtransforms
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
from mpl_toolkits.basemap import Basemap #installed using 
    #conda install -c anaconda basemap
#import scipy.io
#%% IMPORT SOM AND GRIDMET DATA
# change directory and import SOM data from .mat file
numpatterns = 9
percentile = 90

pats_assignments = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_ExtremeDaysand{numpatterns}NodeAssign.npy')
assignment = pats_assignments[:,0]
extreme_griddedprecip = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\GRIDMET_pr_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST.npy')
extreme_griddedprecip = extremedaysprecip
#%%
#IMPORT NETCDF DATA
year = '1980'
#define NC location
nc_file = f'pr_{year}.nc'
nc_dir = 'I:\\GRIDMET\\pr'
filepath = os.path.join(nc_dir,nc_file)
#open the netcdf file in read mode
gridfile = nc.Dataset(filepath,mode='r')
gridlat = gridfile.variables['lat'][:]
gridlon = gridfile.variables['lon'][:]

#REDUCE VARIABLES TO DESIRED AREA
#convert height to a 3D array
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
merrareduced = extreme_griddedprecip[:,latind,:]
merrareduced = extreme_griddedprecip[:,:,lonind]

#convert lat and lon into a 2D array
lon, lat = np.meshgrid(lonreduced,latreduced)
#define area threshold for basemap
area_thresh = 1E4
#define map bounds
latmin, latmax = (min(latreduced),max(latreduced))
lonmin, lonmax = (min(lonreduced),max(lonreduced))

#%% CLUSTER ASSIGNED PATTERNS AND CALCULATE AVERAGE

# create array to store 12 SOM composites
som_composites = np.zeros((numpatterns,len(latreduced),len(lonreduced)))

# loop through all 12 som patterns
for som in range(numpatterns):
    # create array to store assigned days data
    som_merra = np.zeros((1,len(latreduced),len(lonreduced)))
    # loop through all days
    for day,arr in enumerate(merrareduced):
        # add data to som_merra if day is assigned to node
        if assignment[day] == som + 1:
            som_merra = np.concatenate((som_merra,np.expand_dims(arr,axis=0)))
    # remove initial row of zeros
    som_merra = som_merra[1:,:,:]
    # confirm correct number of days assigned to node
    #print(som+1,len(som_merra),pat_freq[som])
    # calculate the mean of assigned days
    som_mean = np.squeeze(np.mean(som_merra,axis=0))
    # append to array of composites
    som_composites[som] = som_mean

#%% DETERMINE MAX AND MIN VALIUES
zmax = 0
zmin = 1E8

for i, arr in enumerate(som_composites):
    
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
fig = plt.figure(figsize=(7.5,5.5))
#fig.suptitle(f'{plottitle[metvar]} Composites',fontsize=13,fontweight="bold",y=0.9875)

lowlim = 0
highlim = 500
colormap = 'gist_ncar_r'

for i, arr in enumerate(som_composites):
    #MAP DESIRED VARIABLE
    #convert lat and lon into a 2D array
    #define area threshold for basemap
    area_thresh = 1E4
    #create equidistant cylindrical projection basemap
    map = Basemap(projection='cyl',llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,\
              urcrnrlon=lonmax,resolution='l',area_thresh=area_thresh)
    xi, yi = map(lon,lat)
    ax = fig.add_subplot(3,3,i+1)
    sublabel_loc = mtransforms.ScaledTranslation(4/72, -4/72, fig.dpi_scale_trans)
    ax.text(0.0, 1.0, i+1, transform=ax.transAxes + sublabel_loc,
        fontsize=9, fontweight='bold', verticalalignment='top', 
        bbox=dict(facecolor='1', edgecolor='none', pad=1.5),zorder=3)
    #create colormap of MERRA2 data
    colorm = map.pcolor(xi,yi,arr,shading='auto',cmap=colormap,vmin=lowlim,vmax=highlim,zorder=1)
    
    #define border color and thickness
    border_c = '0.4'
    border_w = 0.4
    #create map features
    map.drawcoastlines(color=border_c, linewidth=border_w)
    map.drawstates(color=border_c, linewidth=border_w)
    map.drawcountries(color=border_c, linewidth=border_w)
    gridlinefont = 8.5
    parallels = np.arange(37.,43.,1.)
    meridians = np.arange(-125.,-119.,1.)
    if i == 0 or i == 3:
        map.drawparallels(parallels, labels=[1,0,0,0], fontsize=gridlinefont,color=border_c,linewidth=border_w)
        map.drawmeridians(meridians,color=border_c,linewidth=border_w)
    elif i == 7 or i == 8:
        map.drawparallels(parallels, color=border_c,linewidth=border_w)
        map.drawmeridians(meridians, labels=[0,0,0,1], fontsize=gridlinefont,color=border_c,linewidth=border_w)
    elif i == 6:
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
    #map.readshapefile(os.path.join(ws_directory,f'{watershed}'), watershed,linewidth=0.8,color='r')
    plt.scatter(-120.9,39.5,color='tomato',edgecolors='r',marker='*',linewidths=0.8,zorder=4)
    
#CUSTOMIZE SUBPLOT SPACING
fig.subplots_adjust(left=0.05,right=0.89,bottom=0.021, top=0.955,hspace=0.05, wspace=0.05) #bottom colorbar
#fig.add_axis([left,bottom, width,height])
cbar_ax = fig.add_axes([0.904,0.05,0.025,0.88]) #bottom colorbar
#cbar = fig.colorbar(colorm, cax=cbar_ax,ticks=np.arange(cbarstart[metvar],highlims[metvar]+1,cbarint[metvar]),orientation='vertical')
#cbar.ax.tick_params(labelsize=8)
#cbar.set_label(cbarlabs[metvar],fontsize=8.5,labelpad=0.5,fontweight='bold')

    
#SHOW MAP
save_dir='I:\\Emma\\FIROWatersheds\\Figures\\SOMs\\Composites'
os.chdir(save_dir)
#plt.savefig(f'{metvar}_{percentile}_{numpatterns}_SOM_composite.png',dpi=300)
plt.show()