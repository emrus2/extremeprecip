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
#import matplotlib.transforms as mtransforms
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
from mpl_toolkits.basemap import Basemap #installed using 
    #conda install -c anaconda basemap
#import scipy.io
from datetime import datetime
#%% IMPORT EXTREME DAYS DATA
# change directory and import SOM data from .mat file
#mat_dir='I:\\Emma\\FIROWatersheds\\Data\\SOMs\\SomOutput'
numpatterns = 9
percentile = 90

# define lat, lon region of data for plotting
latmin, latmax = (15.5,65.5)
lonmin, lonmax = (-170.25,-105.75)

# import extreme precipitation days and events
data_dir='I:\\Emma\\FIROWatersheds\\Data\\'
os.chdir(data_dir)
extremedays = np.load(f'{percentile}Percentile_ExtremeDays.npy',allow_pickle=True)
alleventdates = np.load(f'{percentile}Percentile_{numpatterns}NodeClusteredEvents.npy',allow_pickle=True)

# find 6-day extreme event
eventlength = 6
for i in alleventdates:
    if len(i) == eventlength:
        extremeevent = i

# find index of first day of extreme event 
eventstart, = np.where(extremedays == str(extremeevent[0]))
eventstart = eventstart[0]
#%% IMPORT MERRA2 DATA
# define metvar
metvars = ['SLP', '300W','Z500Anom','SLPAnom','Z850','850T','850TAnom']
metvars = ['IVT']
#metvar = '300W'
for metvar in metvars:
    #define composite location
    #metvar = input('Enter MERRA-2 Variable: Z500, SLP, 850T, 300W, or IVT \n')
    if metvar == 'Z500Anom':
        folderpath = 'I:\\Emma\\FIROWatersheds\\Data\\DailyMERRA2\\Z500'
        filename = f'MERRA2_Z500_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST.nc'
    elif metvar == 'SLPAnom':
        folderpath = 'I:\\Emma\\FIROWatersheds\\Data\\DailyMERRA2\\SLP'
        filename = f'MERRA2_SLP_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST.nc'
    elif metvar == '850TAnom':
        folderpath = 'I:\\Emma\\FIROWatersheds\\Data\\DailyMERRA2\\850T'
        filename = f'MERRA2_850T_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST.nc'
    else:
        folderpath = f'I:\\Emma\\FIROWatersheds\\Data\\DailyMERRA2\\{metvar}'
        filename = f'MERRA2_{metvar}_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST.nc'
    filepath = os.path.join(folderpath,filename)
    
    
    #COLLECT VARIABLE DATA FROM MERRA2 FILE
    merravar = {'Z500':'H','SLP':'SLP','850T':'T','Z850':'H'}
    #open the netcdf file in read mode
    gridfile = nc.Dataset(filepath,mode='r')
    print(gridfile)
    gridlat = gridfile.variables['lat'][:]
    gridlon = gridfile.variables['lon'][:]
    if metvar == 'IVT':
        Uvapor = gridfile.variables['UFLXQV'][:]
        Vvapor = gridfile.variables['VFLXQV'][:]
        merra = np.sqrt(Uvapor**2 + Vvapor**2)
    elif metvar == '300W':
        Uwind = gridfile.variables['U'][:]
        Vwind = gridfile.variables['V'][:]
        merra = np.sqrt(Uwind**2 + Vwind**2)
    elif 'Anom' in metvar:
        merra = np.load(os.path.join(folderpath,f'MERRA2_{metvar}_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST.npy'))
    else:
        merra = gridfile.variables[merravar[metvar]][:]
    gridfile.close()
    
    merra = np.squeeze(merra)
    
    if metvar == 'SLP':
        merra = merra/100
        
    #%% REDUCE LAT AND LON TO DESIRED AREA
    
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
    merrareduced = merra[:,latind,:]
    merrareduced = merrareduced[:,:,lonind]
    
    #print(np.amin(merrareduced),np.amax(merrareduced))
    
    #%% DETERMINE MAX AND MIN VALUES
    zmax = 0
    zmin = 1E8
    
    for n in np.arange(eventstart,eventstart+eventlength,1):
        
        arr = merrareduced[n,:,:]
        
        #determine zmax and zmin for all days
        highlim = np.nanmax(arr)
        lowlim = np.nanmin(arr)
        print(lowlim,highlim)
        if highlim > zmax:
            zmax = highlim
        if lowlim < zmin:
            zmin = lowlim
    
    print(f'Lowest Value:{zmin} \nHighest Value:{zmax}')
    
    #%% CREATE ANOMALY MAP 
    #GENERATE CUSTOM COLORMAP
    def center_colormap(lowlim, highlim, center=0):
        dv = max(-lowlim, highlim) * 2
        N = int(256 * dv / (highlim-lowlim))
        bwr = cm.get_cmap('seismic', N)
        newcolors = bwr(np.linspace(0, 1, N))
        beg = int((dv / 2 + lowlim)*N / dv)
        end = N - int((dv / 2 - highlim)*N / dv)
        newmap = ListedColormap(newcolors[beg:end])
        return newmap
    
    #%% DEFINE PLOTTING VARIABLES
    if metvar == 'Z500Anom':
        lowanom, highanom = (-2.0, 1.1)
        newmap = center_colormap(lowanom, highanom, center=0)
    elif metvar == 'SLPAnom':
        lowanom, highanom = (-2.4, 0.9)
        newmap = center_colormap(lowanom, highanom, center=0)
    else:
        lowanom, highanom = (-1.3, 1.2)
        newmap = center_colormap(lowanom, highanom, center=0)
    lowlims = {'Z500':2850,'SLP':985,'IVT':0,'300W':0,'850T':252,'Z500Anom':lowanom,'Z850':1187,'SLPAnom':lowanom,'850TAnom':lowanom}
    highlims = {'Z500':5700,'SLP':1022,'IVT':1212,'300W':56,'850T':293,'Z500Anom':highanom,'Z850':1548,'SLPAnom':highanom,'850TAnom':highanom}
    
    contourstart = {'Z500':3000,'SLP':990,'IVT':0,'300W':5,'850T':250,'Z500Anom':-1.75,'Z850':1190,'SLPAnom':-2.25,'850TAnom':-1.2}
    contourint = {'Z500':200,'SLP':4,'IVT':100,'300W':5,'850T':2.5,'Z500Anom':0.25,'Z850':30,'SLPAnom':0.25,'850TAnom':0.15}
    
    cbarstart = {'Z500':3000,'SLP':990,'IVT':0,'300W':0,'850T':250,'Z500Anom':-2.0,'Z850':1200,'SLPAnom':-2.4,'850TAnom':-1.2}
    cbarint = {'Z500':500,'SLP':5,'IVT':150,'300W':10,'850T':5,'Z500Anom':0.5,'Z850':50,'SLPAnom':0.4,'850TAnom':0.3}
    
    colormap = {'Z500':'jet','SLP':'rainbow','IVT':'gnuplot2_r','300W':'hot_r','850T':'turbo','Z500Anom':newmap,'Z850':'turbo','SLPAnom':newmap,'850TAnom':newmap}
    cbarlabs = {'Z500':'m','SLP':'hPa','IVT':'kg $\mathregular{m^{-1}}$ $\mathregular{s^{-1}}$','300W':'m/s','850T':'K','Z500Anom':r'$\mathbf{\sigma}$','Z850':'m','SLPAnom':r'$\mathbf{\sigma}$','850TAnom':r'$\mathbf{\sigma}$'}
    plottitle = {'Z500':'Z500','SLP':'SLP','IVT':'IVT','300W':'300 hPa Wind','850T':'850 hPa Temperature','Z500Anom':'Z500 Anomaly','Z850':'Z850','SLPAnom':'SLP Anomaly','850TAnom':'850 hPa Temperature Anomaly'}
    #%% PLOT NODES from MATLAB
    
    #create subplot for mapping multiple timesteps
    fig = plt.figure(figsize=(7.2,4))
    i = 0
    for n in np.arange(eventstart,eventstart+eventlength,1):
        #MAP DESIRED VARIABLE
        # define date of plot
        datetitle = datetime.strptime(str(extremeevent[i]),'%Y%m%d')
        # reduce merra to desired day
        arr = merrareduced[n,:,:]
        #convert lat and lon into a 2D array
        lon, lat = np.meshgrid(gridlonreduced,gridlatreduced) 
        #define area threshold for basemap
        area_thresh = 1E4
        #create equidistant cylindrical projection basemap
        map = Basemap(projection='cyl',llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,\
                  urcrnrlon=lonmax,resolution='l',area_thresh=area_thresh)
        xi, yi = map(lon,lat)
        ax = fig.add_subplot(2,3,i+1)
        ax.set_title('{:%d %b}'.format(datetitle),pad=4,fontsize=12)
        # sublabel_loc = mtransforms.ScaledTranslation(4/72, -4/72, fig.dpi_scale_trans)
        # ax.text(0.0, 1.0, extremeevent[i], transform=ax.transAxes + sublabel_loc,
        #     fontsize=9, fontweight='bold', verticalalignment='top', 
        #     bbox=dict(facecolor='1', edgecolor='none', pad=1.5),zorder=3)
        #create colormap of MERRA2 data
        colorm = map.pcolor(xi,yi,arr,shading='auto',cmap=colormap[metvar],vmin=lowlims[metvar],vmax=highlims[metvar],zorder=1)
        
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
        if i in np.arange(0,10,3):
            map.drawparallels(parallels, labels=[1,0,0,0], fontsize=gridlinefont,color=border_c,linewidth=border_w)
            map.drawmeridians(meridians, labels=[0,0,0,1], fontsize=gridlinefont,color=border_c,linewidth=border_w)
        else:
            map.drawparallels(parallels, color=border_c,linewidth=border_w)
            map.drawmeridians(meridians, labels=[0,0,0,1], fontsize=gridlinefont,color=border_c,linewidth=border_w)
        #define contour color and thickness
        contour_c = '0.1'
        contour_w = 0.7
        #create contour map
        contourm = map.contour(xi,yi,arr,colors=contour_c,linewidths=contour_w,levels=np.arange(contourstart[metvar],highlims[metvar]+1,contourint[metvar]),zorder=2)
        plt.clabel(contourm,levels=contourm.levels[::2],fontsize=6,inline_spacing=1,colors='k',zorder=2,manual=False)
            
        #add yuba shape
        #map.readshapefile(os.path.join(ws_directory,f'{watershed}'), watershed,linewidth=0.8,color='r')
        plt.scatter(-120.9,39.5,color='w',marker='*',linewidths=0.7,zorder=4)
        
        i += 1
        
    #CUSTOMIZE SUBPLOT SPACING
    fig.subplots_adjust(left=0.05,right=0.89,bottom=0.021, top=0.955,hspace=0.05, wspace=0.05) #bottom colorbar
    #fig.add_axis([left,bottom, width,height])
    cbar_ax = fig.add_axes([0.9,0.05,0.025,0.88]) #bottom colorbar
    cbar = fig.colorbar(colorm, cax=cbar_ax,ticks=np.arange(cbarstart[metvar],highlims[metvar]+1,cbarint[metvar]),orientation='vertical')
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label(cbarlabs[metvar],fontsize=8.5,labelpad=0.5,fontweight='bold')
    
        
    #SHOW MAP
    save_dir='I:\\Emma\\FIROWatersheds\\Figures\\Events'
    os.chdir(save_dir)
    plt.savefig(f'{eventlength}day_event_{metvar}.png',dpi=300)
    plt.show()