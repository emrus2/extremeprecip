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
from palettable.colorbrewer.diverging import BrBG_10
from palettable.cmocean.sequential import Tempo_20
from palettable.cmocean.sequential import Algae_20
from palettable.colorbrewer.sequential import Greens_9
from palettable.colorbrewer.sequential import BuGn_9
from palettable.scientific.sequential import Davos_20_r

#%% DEFINE WATERSHED SHAPEFILE
watershed = 'UpperYuba'
ws_directory = f'I:\\Emma\\FIROWatersheds\\Data\\WatershedShapefiles\\California\\{watershed}\\'

# determine if you want to plot IVT data
IVT = True
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
filepath = ('I:\\Emma\\FIROWatersheds\\Data\\Gridmet\\GRIDMET_pr_Yuba_Extremes90_DailyAnomProp_1980-2021_WINTERDIST.nc')
#open the netcdf file in read mode
gridfile = nc.Dataset(filepath,mode='r')
print(gridfile)
gridlat = gridfile.variables['lat'][:]
gridlon = gridfile.variables['lon'][:]
precip = np.squeeze(gridfile.variables['precipitation_amount'][:])
days = gridfile.variables['day'][:]
gridfile.close()
#dates = [datetime.strftime(datetime(1900,1,1)+timedelta(days=days[n]),"%Y%m%d") for n in range(days.shape[0])]

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
    
#%% CLUSTER ASSIGNED PATTERNS AND CALCULATE AVERAGE

# create array to store 12 SOM composites
som_composites = np.zeros((numpatterns,len(latreduced),len(lonreduced)))

# loop through all 12 som patterns
for som in range(numpatterns):
    # create array to store assigned days data
    som_merra = np.zeros((1,len(latreduced),len(lonreduced)))
    # loop through all days
    for day,arr in enumerate(merrareduced):
        print(np.amax(arr))
        # add data to som_merra if day is assigned to node
        if assignment[day] == float(som + 1):
            som_merra = np.ma.concatenate((som_merra,np.expand_dims(arr,axis=0)))
            print(np.amax(som_merra))
    # remove initial row of zeros
    som_merra = som_merra[1:,:,:]
    # calculate the mean of assigned days
    som_mean = np.squeeze(np.mean(som_merra,axis=0))
    # append to array of composites
    som_composites[som] = som_mean
#%% DETERMINE MAX AND MIN VALIUES
zmax = 0
zmin = 1E8

som_composites[som_composites==0.] = np.nan

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

#%%
if IVT == True:
    # IMPORT MERRA2 DATA
    metvar = 'IVT'
    merrapath = f'I:\\Emma\\FIROWatersheds\\Data\\DailyMERRA2\\{metvar}'
    merraname = f'MERRA2_{metvar}_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST.nc'
    merrafile = os.path.join(merrapath,merraname)

    gridmerra = nc.Dataset(merrafile,mode='r')
    print(gridmerra)
    gridlatm = gridmerra.variables['lat'][:]
    gridlonm = gridmerra.variables['lon'][:]
    Uvapor = gridmerra.variables['UFLXQV'][:]
    Vvapor = gridmerra.variables['VFLXQV'][:]
    merra = np.sqrt(Uvapor**2 + Vvapor**2)
    gridmerra.close()

    #INCLUDE CALCULATION OF IVT VECTORS
    latlimsm = np.logical_and(gridlatm > latmin, gridlatm < latmax)
    latindm = np.where(latlimsm)[0]
    latreducedm = gridlatm[latindm]
    #reduce lon
    lonlimsm = np.logical_and(gridlonm > lonmin, gridlonm < lonmax)
    lonindm = np.where(lonlimsm)[0]
    lonreducedm = gridlonm[lonindm]
    #reduce pressure
    Uvaporreduced = Uvapor[:,latindm,:]
    Uvaporreduced = Uvaporreduced[:,:,lonindm]
    Vvaporreduced = Vvapor[:,latindm,:]
    Vvaporreduced = Vvaporreduced[:,:,lonindm]

    lonm, latm = np.meshgrid(lonreducedm,latreducedm)

    # create array to store 12 SOM composites
    U_composites = np.zeros((numpatterns,len(latreducedm),len(lonreducedm)))
    V_composites = np.zeros((numpatterns,len(latreducedm),len(lonreducedm)))
    # loop through all 12 som patterns
    for som in range(numpatterns):
        # create array to store assigned days data
        U_merra = np.zeros((1,len(latreducedm),len(lonreducedm)))
        V_merra = np.zeros((1,len(latreducedm),len(lonreducedm)))
        # loop through all days
        for day,arr in enumerate(merrareduced):
            U_array = Uvaporreduced[day,:,:]
            V_array = Vvaporreduced[day,:,:]
            # add data to som_merra if day is assigned to node
            if assignment[day] == float(som + 1):
                U_merra = np.concatenate((U_merra,np.expand_dims(U_array,axis=0)))
                V_merra = np.concatenate((V_merra,np.expand_dims(V_array,axis=0)))
        # remove initial row of zeros
        U_merra = U_merra[1:,:,:]
        V_merra = V_merra[1:,:,:]
        # confirm correct number of days assigned to node
        #print(som+1,len(U_merra),pat_freq[som])
        # calculate the mean of assigned days
        U_mean = np.squeeze(np.mean(U_merra,axis=0))
        V_mean = np.squeeze(np.mean(V_merra,axis=0))
        # append to array of composites
        U_composites[som] = U_mean
        V_composites[som] = V_mean

#%% PLOT NODES from MATLAB

#create subplot for mapping multiple timesteps
fig = plt.figure(figsize=(7.2,5.1))
fig.suptitle('b)',fontsize=11,fontweight="bold",y=0.997,x=0.07)

lowlim = 3.3
highlim = 18.4
#colormap = 'gist_ncar_r'
colormap = BrBG_10.mpl_colormap

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
        fontsize=10, fontweight='bold', verticalalignment='top', 
        bbox=dict(facecolor='1', edgecolor='none', pad=1.5),zorder=5)
    #create colormap of MERRA2 data
    colorm = map.pcolor(xi,yi,arr,shading='auto',cmap=colormap,vmin=lowlim,vmax=highlim,zorder=1)
    
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
        map.drawparallels(parallels, labels=[1,0,0,0], fontsize=gridlinefont,color=border_c,linewidth=border_w,zorder=3)
        map.drawmeridians(meridians,color=border_c,linewidth=border_w,zorder=3)
    elif i == 7 or i == 8:
        map.drawparallels(parallels, color=border_c,linewidth=border_w,zorder=3)
        map.drawmeridians(meridians, labels=[0,0,0,1], fontsize=gridlinefont,color=border_c,linewidth=border_w,zorder=3)
    elif i == 6:
        map.drawparallels(parallels, labels=[1,0,0,0], fontsize=gridlinefont,color=border_c,linewidth=border_w,zorder=3)
        map.drawmeridians(meridians, labels=[0,0,0,1], fontsize=gridlinefont,color=border_c,linewidth=border_w,zorder=3)
    else:
        map.drawparallels(parallels, color=border_c,linewidth=border_w,zorder=3)
        map.drawmeridians(meridians,color=border_c,linewidth=border_w,zorder=3)
    #define contour color and thickness
    contour_c = '0.1'
    contour_w = 0.7
    
    if IVT == True:
        #plot IVT vectors
        U_arrs = U_composites[i,:,:]
        V_arrs = V_composites[i,:,:]
        interval = 1
        size = 1500
        skip = (slice(None, None, interval), slice(None, None, interval))
        xim, yim = map(lonm,latm)
        #vectorm = map.quiver(xi2[skip],yi2[skip],U_arrs[skip],V_arrs[skip],color='darkgreen')
        vectorm = map.quiver(xim[skip],yim[skip],U_arrs[skip],V_arrs[skip],pivot='mid',zorder=4, \
                             scale=size, scale_units='inches',headlength=5,headwidth=3,color='b',width=0.007,alpha=0.7)

    #create contour map
    #contourm = map.contour(xi,yi,arr,colors=contour_c,linewidths=contour_w,levels=np.arange(contourstart[metvar],highlims[metvar]+1,contourint[metvar]),zorder=2)
    #plt.clabel(contourm,levels=contourm.levels[::2],fontsize=6,inline_spacing=1,colors='k',zorder=2,manual=False)
        
    #add yuba shape
    map.readshapefile(os.path.join(ws_directory,f'{watershed}'), watershed,zorder=3)
    #plt.scatter(-120.9,39.5,color='tomato',edgecolors='r',marker='*',linewidths=0.8,zorder=4)
    
#CUSTOMIZE SUBPLOT SPACING
fig.subplots_adjust(left=0.05,right=0.9,bottom=0.026, top=0.968,hspace=0.05, wspace=0.05) #bottom colorbar
#fig.add_axis([left,bottom, width,height])
cbar_ax = fig.add_axes([0.91,0.05,0.025,0.9]) #bottom colorbar
cbar = fig.colorbar(colorm, cax=cbar_ax,ticks=np.arange(4,highlim,2),orientation='vertical')
cbar.ax.tick_params(labelsize=8)
cbar.set_label('Proportion of Monthly Average',fontsize=8.5,labelpad=0.5,fontweight='bold')

    
#SHOW MAP
save_dir='I:\\Emma\\FIROWatersheds\\Figures\\SOMs\\Composites'
os.chdir(save_dir)
if IVT == True:
    plt.savefig('GRIMET_PrecipitiationMonthlyProportion_IVT_SOM_composite_lettered.png',dpi=300)
else:
    plt.savefig('GRIMET_PrecipitiationMonthlyProportion_SOM_composite_lettered.png',dpi=300)
plt.show()