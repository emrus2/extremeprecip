# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 13:49:49 2023

@author: emrus2

Plotting ACTUAL IVT during extreme precipitation events to compare to SOM patterns

UPDATED 6/12/2023
"""
#%% IMPORT MODULES
import netCDF4 as nc #if this doesn't work, try to reinstall in anaconda prompt using
    #conda install netcdf4
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.cm as cm
import matplotlib.transforms as mtransforms
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
from mpl_toolkits.basemap import Basemap #installed using 
    #conda install -c anaconda basemap
#import scipy.io
from datetime import datetime
#%% IMPORT EXTREME DAYS DATA
data_dir='I:\\Emma\\FIROWatersheds\\Data\\'
os.chdir(data_dir)
extremedays = np.load('ExtremePrecipitationDays_2022_23.npy',allow_pickle=True)


# define lat, lon region of data for plotting
latmin, latmax = (15.5,65.5)
lonmin, lonmax = (-170.25,-105.75)

percentile = 90
#%% IMPORT MERRA2 DATA
# define metvar
metvars = ['IVT']
folderpath = 'I:\\MERRA2\\Daily_and_Subdaily\\IVT_hourly_new\\'
#%%
fig = plt.figure(figsize=(8,4))

for i,arr in enumerate(extremedays):
    # arr = extremedays[0]
    date = arr[0]
    dt = datetime.strptime(date,'%Y%m%d')
    precip = arr[1]
    filename = f'MERRA2.tavg1_2d_int_Nx.{date}.SUB.nc'
    filepath = os.path.join(folderpath,filename)
    
    #COLLECT VARIABLE DATA FROM MERRA2 FILE
    gridfile = nc.Dataset(filepath,mode='r')
    print(gridfile)
    gridlat = gridfile.variables['lat'][:]
    gridlon = gridfile.variables['lon'][:]
    Uvapor = gridfile.variables['UFLXQV'][:]
    Vvapor = gridfile.variables['VFLXQV'][:]
    Uvaporavg = Uvapor.mean(axis=0)
    Vvaporavg = Vvapor.mean(axis=0)
    merra = np.sqrt(Uvaporavg**2 + Vvaporavg**2)
    gridfile.close()
    
    merra = np.squeeze(merra)
        
    
    #REDUCE VARIABLES TO DESIRED AREA
    #reduce lat
    latlims = np.logical_and(gridlat > latmin, gridlat < latmax)
    latind = np.where(latlims)[0]
    gridlatreduced = gridlat[latind]
    #reduce lon
    lonlims = np.logical_and(gridlon > lonmin, gridlon < lonmax)
    lonind = np.where(lonlims)[0]
    gridlonreduced = gridlon[lonind]
    #reduce pressure
    merrareduced = merra[latind,:]
    merrareduced = merrareduced[:,lonind]
    
    print(np.amin(merrareduced),np.amax(merrareduced))

    # DEFINE PLOTTING VARIABLES
    lowlims, highlims = (0,975)
    contourstart, contourint = (0,100)
    cbarstart, cbarint = (0,150)
    colormap = 'gnuplot2_r'
    cbarlabs = 'kg/m/s'
    plottitle = 'IVT'
    
    # PLOT NODES from MATLAB
    
    #create subplot for mapping multiple timesteps
    ax = fig.add_subplot(2,4,i+1)
    #MAP DESIRED VARIABLE
    #convert lat and lon into a 2D array
    lon, lat = np.meshgrid(gridlonreduced,gridlatreduced) 
    #define area threshold for basemap
    area_thresh = 1E4
    #create equidistant cylindrical projection basemap
    map = Basemap(projection='cyl',llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,\
              urcrnrlon=lonmax,resolution='l',area_thresh=area_thresh)
    xi, yi = map(lon,lat)
    ax.set_title('{:%d %b %Y}'.format(dt),pad=4,fontsize=12)
    # sublabel color
    sublabel_loc = mtransforms.ScaledTranslation(4/72, -4/72, fig.dpi_scale_trans)
    ax.text(0.63, 1.0, round(float(precip),2), transform=ax.transAxes + sublabel_loc,
        fontsize=9, fontweight='bold', verticalalignment='top', 
        bbox=dict(facecolor=(0,1,0), edgecolor='none', pad=1.5),zorder=3)
    #create colormap of MERRA2 data
    colorm = map.pcolor(xi,yi,merrareduced,shading='auto',cmap=colormap,vmin=lowlims,vmax=highlims,zorder=1)
    
    #define border color and thickness
    border_c = '0.4'
    border_w = 0.4
    #create map features
    map.drawcoastlines(color=border_c, linewidth=border_w)
    map.drawstates(color=border_c, linewidth=border_w)
    map.drawcountries(color=border_c, linewidth=border_w)
    gridlinefont = 8.5
    parallels = np.arange(20.,71.,20.)
    meridians = np.arange(-160.,-109.,20.)
    if i in np.arange(0,10,4):
        map.drawparallels(parallels, labels=[1,0,0,0], fontsize=gridlinefont,color=border_c,linewidth=border_w)
        map.drawmeridians(meridians, labels=[0,0,0,1], fontsize=gridlinefont,color=border_c,linewidth=border_w)
    else:
        map.drawparallels(parallels, color=border_c,linewidth=border_w)
        map.drawmeridians(meridians, labels=[0,0,0,1], fontsize=gridlinefont,color=border_c,linewidth=border_w)
    #define contour color and thickness
    contour_c = '0.1'
    contour_w = 0.7
    #create contour map
    contourm = map.contour(xi,yi,merrareduced,colors=contour_c,linewidths=contour_w,levels=np.arange(contourstart,highlims+1,contourint),zorder=2)
    plt.clabel(contourm,levels=contourm.levels[::2],fontsize=6,inline_spacing=1,colors='k',zorder=2,manual=False)
        
    #add yuba shape
    #map.readshapefile(os.path.join(ws_directory,f'{watershed}'), watershed,linewidth=0.8,color='r')
    plt.scatter(-120.9,39.5,color='w',marker='*',linewidths=0.7,zorder=4)

fig.subplots_adjust(left=0.05,right=0.89,bottom=0.01, top=0.97,hspace=0.05, wspace=0.05) #bottom colorbar

#CUSTOMIZE SUBPLOT SPACING
#fig.add_axis([left,bottom, width,height])
cbar_ax = fig.add_axes([0.9,0.05,0.025,0.88]) #bottom colorbar
cbar = fig.colorbar(colorm, cax=cbar_ax,ticks=np.arange(cbarstart,highlims+1,cbarint),orientation='vertical')
cbar.ax.tick_params(labelsize=8)
cbar.set_label(cbarlabs,fontsize=8.5,labelpad=0.5,fontweight='bold')


#SHOW MAP
save_dir='I:\\Emma\\FIROWatersheds\\Figures\\Events'
os.chdir(save_dir)
plt.savefig('202223_ExtremeDays_IVT.png',dpi=300,bbox_inches='tight')
plt.show()