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
bigarrows = True

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
#create directories for varying met data
#metvars = ['precip','mintemp','maxtemp','wind']
metvar = 'wind'

metabbr = {'precip':'pr','mintemp':'tmmn','maxtemp':'tmmx','wind':'vs'}
varnames = {'precip':'precipitation_amount','mintemp':'air_temperature','maxtemp':'air_temperature','wind':'wind_speed'}
plottitles = {'precip':'Precipitation','mintemp':'Minimum Air Temperature','maxtemp':'Maximum Air Temperature','wind':'10m Wind'}

lowlims = {'precip':0,'mintemp':-12,'maxtemp':-12,'wind':0}
highlims = {'precip':150,'mintemp':17,'maxtemp':17,'wind':15}

cbarstart = {'precip':0,'mintemp':-12,'maxtemp':-12,'wind':0}
cbarint = {'precip':25,'mintemp':4,'maxtemp':4,'wind':3}
    
colormap = {'precip':'gist_ncar_r','mintemp':'turbo','maxtemp':'turbo','wind':'hot_r'}
cbarlabs = {'precip':'mm','mintemp':u'\N{DEGREE SIGN}C','maxtemp':u'\N{DEGREE SIGN}C','wind':'m/s'}

#%%
#IMPORT NETCDF DATA
#define NC location
filepath = (f'I:\\Emma\\FIROWatersheds\\Data\\Gridmet\\GRIDMET_{metabbr[metvar]}_Yuba_Extremes90_Daily_1980-2021_WINTERDIST.nc')
#open the netcdf file in read mode
gridfile = nc.Dataset(filepath,mode='r')
print(gridfile)
gridlat = gridfile.variables['lat'][:]
gridlon = gridfile.variables['lon'][:]
griddata = np.squeeze(gridfile.variables[varnames[metvar]][:])
days = gridfile.variables['day'][:]
gridfile.close()
#dates = [datetime.strftime(datetime(1900,1,1)+timedelta(days=days[n]),"%Y%m%d") for n in range(days.shape[0])]

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
merrareduced = griddata[:,latind,:]
merrareduced = merrareduced[:,:,lonind]

#convert lat and lon into a 2D array
lon, lat = np.meshgrid(lonreduced,latreduced)
#define area threshold for basemap
area_thresh = 1E4

#%% CLUSTER ASSIGNED PATTERNS AND CALCULATE AVERAGE

# create array to store 12 SOM composites
som_composites = np.zeros((numpatterns,len(latreduced),len(lonreduced)))

# loop through all 12 som patterns
for som in range(numpatterns):
    # create array to store assigned days data
    som_merra = np.zeros((1,len(latreduced),len(lonreduced)))
    # loop through all days
    for day,arr in enumerate(merrareduced):
        #print(np.amax(arr))
        # add data to som_merra if day is assigned to node
        if assignment[day] == float(som + 1):
            som_merra = np.ma.concatenate((som_merra,np.expand_dims(arr,axis=0)))
            #print(np.amax(som_merra))
    # remove initial row of zeros
    som_merra = som_merra[1:,:,:]
    # calculate the mean of assigned days
    som_mean = np.squeeze(np.mean(som_merra,axis=0))
    # append to array of composites
    som_composites[som] = som_mean
    
#convert zeroes to NaNs
som_composites[som_composites == 0] = 'nan'

#convert K to C
if 'temp' in metvar:
    som_composites -= 273.15
    
#%% INCLUDE WIND DIRECTION ALSO
if metvar == 'wind':
    filepath2 = ('I:\\Emma\\FIROWatersheds\\Data\\Gridmet\\GRIDMET_th_Yuba_Extremes90_Daily_1980-2021_WINTERDIST.nc')
    gridfile2 = nc.Dataset(filepath2,mode='r')
    griddata2 = np.squeeze(gridfile2.variables['wind_from_direction'][:])
    gridfile2.close()
    winddir = griddata2[:,latind,:]
    winddir = winddir[:,:,lonind]
    
    # create array to store 12 SOM composites
    som_composites_winddir = np.zeros((numpatterns,len(latreduced),len(lonreduced)))

    # loop through all 12 som patterns
    for som in range(numpatterns):
        # create array to store assigned days data
        som_merra = np.zeros((1,len(latreduced),len(lonreduced)))
        # loop through all days
        for day,arr in enumerate(winddir):
            #print(np.amax(arr))
            # add data to som_merra if day is assigned to node
            if assignment[day] == float(som + 1):
                som_merra = np.ma.concatenate((som_merra,np.expand_dims(arr,axis=0)))
                print(np.amin(som_merra),np.amax(som_merra))
        # remove initial row of zeros
        som_merra = som_merra[1:,:,:]
        # calculate the mean of assigned days
        som_mean = np.squeeze(np.mean(som_merra,axis=0))
        # append to array of composites
        som_composites_winddir[som] = som_mean
    
    #convert zeroes to NaNs
    #som_composites_winddir[som_composites_winddir == 0] = 'nan'
    #print(np.nanmin(som_composites_winddir),np.nanmax(som_composites_winddir))
    
    winddir_radians = np.radians(270 - som_composites_winddir)
    Uwind = som_composites*np.cos(winddir_radians)
    Vwind = som_composites*np.sin(winddir_radians)

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
fig = plt.figure(figsize=(7.2,5))
#fig.suptitle(f'{plottitles[metvar]} Composites',fontsize=13,fontweight="bold",y=0.9875)

for i, arr in enumerate(som_composites):
    Uarr = Uwind[i,:,:]
    Varr = Vwind[i,:,:]
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
        bbox=dict(facecolor='1', edgecolor='none', pad=1.5),zorder=5)
    #create colormap of MERRA2 data
    colorm = map.pcolor(xi,yi,arr,shading='auto',cmap=colormap[metvar],vmin=lowlims[metvar],vmax=highlims[metvar],zorder=1)
    
    #define border color and thickness
    border_c = '0.4'
    border_w = 0.4
    #create map features
    map.drawcoastlines(color=border_c, linewidth=border_w,zorder=2)
    map.drawstates(color=border_c, linewidth=border_w,zorder=2)
    map.drawcountries(color=border_c, linewidth=border_w,zorder=2)
    gridlinefont = 8.5
    parallels = np.arange(38.,42.,1.)
    meridians = np.arange(-124.,-119.,2.)
    if i == 0 or i == 3:
        map.drawparallels(parallels, labels=[1,0,0,0], fontsize=gridlinefont,color=border_c,linewidth=border_w,zorder=2)
        map.drawmeridians(meridians,color=border_c,linewidth=border_w,zorder=2)
    elif i == 7 or i == 8:
        map.drawparallels(parallels, color=border_c,linewidth=border_w,zorder=2)
        map.drawmeridians(meridians, labels=[0,0,0,1], fontsize=gridlinefont,color=border_c,linewidth=border_w,zorder=2)
    elif i == 6:
        map.drawparallels(parallels, labels=[1,0,0,0], fontsize=gridlinefont,color=border_c,linewidth=border_w,zorder=2)
        map.drawmeridians(meridians, labels=[0,0,0,1], fontsize=gridlinefont,color=border_c,linewidth=border_w,zorder=2)
    else:
        map.drawparallels(parallels, color=border_c,linewidth=border_w,zorder=2)
        map.drawmeridians(meridians,color=border_c,linewidth=border_w,zorder=2)
    #define contour color and thickness
    #contour_c = '0.1'
    contour_c = 'b'
    contour_w = 0.7
    #create contour map
    if metvar == 'wind':
        if bigarrows == True:
            interval = 8
            size = 1000
            skip = (slice(None, None, interval), slice(None, None, interval))
            #vectorm = map.quiver(xi2[skip],yi2[skip],U_arrs[skip],V_arrs[skip],color='darkgreen')
            vectorm = map.quiver(xi[skip],yi[skip],Uarr[skip],Varr[skip],pivot='mid',color='b',zorder=4)
        else:
            interval = 5
            size = 1500
            skip = (slice(None, None, interval), slice(None, None, interval))
            #vectorm = map.quiver(xi2[skip],yi2[skip],U_arrs[skip],V_arrs[skip],color='darkgreen')
            vectorm = map.quiver(xi[skip],yi[skip],Uarr[skip],Varr[skip],pivot='mid',color='b',zorder=4)
    if metvar == 'mintemp':
        contourm = map.contour(xi,yi,arr,colors=contour_c,linewidths=contour_w,levels=np.arange(0,1,1),zorder=4)
    #plt.clabel(contourm,levels=contourm.levels[::2],fontsize=6,inline_spacing=1,colors='k',zorder=2,manual=False)
        
    #add yuba shape
    map.readshapefile(os.path.join(ws_directory,f'{watershed}'), watershed,linewidth=0.6,zorder=3)
    #plt.scatter(-120.9,39.5,color='tomato',edgecolors='r',marker='*',linewidths=0.8,zorder=4)
    
#CUSTOMIZE SUBPLOT SPACING
fig.subplots_adjust(left=0.05,right=0.9,bottom=0.026, top=0.985,hspace=0.05, wspace=0.05) #bottom colorbar
#fig.add_axis([left,bottom, width,height])
cbar_ax = fig.add_axes([0.91,0.05,0.025,0.9]) #bottom colorbar
cbar = fig.colorbar(colorm, cax=cbar_ax,ticks=np.arange(cbarstart[metvar],highlims[metvar]+1,cbarint[metvar]),orientation='vertical')
cbar.ax.tick_params(labelsize=8)
cbar.set_label(cbarlabs[metvar],fontsize=8.5,labelpad=0.5,fontweight='bold')

    
#SHOW MAP
save_dir='I:\\Emma\\FIROWatersheds\\Figures\\SOMs\\Composites'
os.chdir(save_dir)
if metvar == 'wind' and bigarrows == True:
    plt.savefig(f'GRIMET_{metvar}_SOM_composite_biggerarrows.png',dpi=300)
else:
    plt.savefig(f'GRIMET_{metvar}_SOM_composite',dpi=300)
plt.show()