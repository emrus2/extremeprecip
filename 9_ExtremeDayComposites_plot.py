# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 14:13:40 2022

@author: emrus2

Plots composite mean maps for MERRA2 extreme precip days


UPDATED ON 7/11/2023

"""

#IMPORT MODULES
import netCDF4 as nc #if this doesn't work, try to reinstall in anaconda prompt using
    #conda install netcdf4
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.cm as cm
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
from mpl_toolkits.basemap import Basemap #installed using 
    #conda install -c anaconda basemap
    
watershed = 'UpperYuba'
ws_directory = f'I:\\Emma\\FIROWatersheds\\Data\\WatershedShapefiles\\California\\{watershed}\\'
    
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

#%% IMPORT MERRA2

percentile = 90

#define composite location
#metvar = input('Enter MERRA-2 Variable: Z500, SLP, 850T, 300W, Z500Anom, SLPAnom or IVT \n')
metvars = ['IVT', 'SLP', 'SLPAnom', 'Z500Anom', '300W','850T']
metvars = ['300W']
#metvar = '850T'
for metvar in metvars:
    #define composite location
    if metvar == 'Z500Anom':
        folderpath = 'I:\\Emma\\FIROWatersheds\\Data\\DailyMERRA2\\Z500'
        filename = f'MERRA2_Z500_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST.nc'
    elif metvar == 'SLPAnom':
        folderpath = 'I:\\Emma\\FIROWatersheds\\Data\\DailyMERRA2\\SLP'
        filename = f'MERRA2_SLP_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST.nc'
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
    elif metvar == 'Z500Anom':
        merra = np.load(os.path.join(folderpath,f'MERRA2_Z500Anom_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST.npy'))
    elif metvar == 'SLPAnom':
        merra = np.load(os.path.join(folderpath,f'MERRA2_SLPAnom_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST.npy'))
    else:
        merra = gridfile.variables[merravar[metvar]][:]
    gridfile.close()
    
    merra = np.squeeze(merra)
    
    if metvar == 'SLP':
        merra = merra/100
    
    #%% REDUCE MERRA2
    
    #calculate mean merra variables
    merramean = np.squeeze(np.mean(merra,axis=0))
        
    #REDUCE VARIABLES TO DESIRED AREA
    #define lat lon restrictions
    latmin, latmax = (15.5,65.5)
    lonmin, lonmax = (-170.25,-105.75)
    #reduce lat
    latlims = np.logical_and(gridlat > latmin, gridlat < latmax)
    latind = np.where(latlims)[0]
    gridlatreduced = gridlat[latind]
    #reduce lon
    lonlims = np.logical_and(gridlon > lonmin, gridlon < lonmax)
    lonind = np.where(lonlims)[0]
    gridlonreduced = gridlon[lonind]
    #reduce pressure
    merrareduced = merramean[latind,:]
    merrareduced = merrareduced[:,lonind]
    
    print(np.amin(merrareduced),np.amax(merrareduced))
    
    #%%
    #Add Z850 contours
    if metvar == '850T':
        #include Z850 contours
        folderpath = 'I:\\Emma\\FIROWatersheds\\Data\\DailyMERRA2\\Z850'
        filename = f'MERRA2_Z850_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST.nc'
        filepath = os.path.join(folderpath,filename)
        gridfile = nc.Dataset(filepath,mode='r')
        print(gridfile)
        merra_contours = gridfile.variables['H'][:]
        gridfile.close()
        #calculate mean
        merramean_contours = np.squeeze(np.mean(merra_contours,axis=0))
        #reduce area
        merrareduced_contours = merramean_contours[latind,:]
        merrareduced_contours = merrareduced_contours[:,lonind]
        
        print(np.amin(merrareduced_contours),np.amax(merrareduced_contours))
        #1255.35 and 1633.44
    
    #%% PLOT VARIABLE
    if percentile == 95:
        #create anomaly colormap for Z500
        #lowanom, highanom = (-240, 40)
        if metvar == 'Z500Anom':
            lowanom, highanom = (-1.6, 0.5)
            newmap = center_colormap(lowanom, highanom, center=0)
        else:
            lowanom, highanom = (-2.0, 0.6)
            newmap = center_colormap(lowanom, highanom, center=0)
        #create subplot for mapping multiple timesteps
        #for regular Z500: low:5110, high:5851, start:5200, int:70, labs:m, int: 150
        #for regular Z500: low:-240, high:40, start:-200, int:30, labs:m, int: 50
    
        colormap = {'Z500Anom':newmap,'SLP':'rainbow','IVT':'gnuplot2_r','300W':'hot_r','850T':'turbo','SLPAnom':newmap}
        lowlims = {'Z500Anom':lowanom,'SLP':994,'IVT':0,'300W':0,'850T':250,'SLPAnom':lowanom}
        highlims = {'Z500Anom':highanom,'SLP':1020,'IVT':530,'300W':48,'850T':294,'SLPAnom':highanom}
        
        contourstart = {'Z500Anom':-1.4,'SLP':996,'IVT':0,'300W':10,'850T':240,'SLPAnom':-1.8}
        contourint = {'Z500Anom':0.25,'SLP':3,'IVT':50,'300W':4,'850T':3,'SLPAnom':0.25}
        
        cbarstart = {'Z500Anom':-1.5,'SLP':995,'IVT':0,'300W':0,'850T':240,'SLPAnom':-2.0}
        cbarint = {'Z500Anom':0.5,'SLP':5,'IVT':100,'300W':10,'850T':10,'SLPAnom':0.5}
        
    elif percentile == 90:  
        #create anomaly colormap for Z500
        #lowanom, highanom = (-240, 40)
        if metvar == 'Z500Anom':
            lowanom, highanom = (-1.6, 0.4)
            newmap = center_colormap(lowanom, highanom, center=0)
        else:
            lowanom, highanom = (-1.9, 0.5)
            newmap = center_colormap(lowanom, highanom, center=0)
        #create subplot for mapping multiple timesteps
        #for regular Z500: low:5110, high:5851, start:5200, int:70, labs:m, int: 150
        #for regular Z500: low:-240, high:40, start:-200, int:30, labs:m, int: 50
    
        colormap = {'Z500Anom':newmap,'SLP':'rainbow','IVT':'gnuplot2_r','300W':'hot_r','850T':'turbo','SLPAnom':newmap}
        lowlims = {'Z500Anom':lowanom,'SLP':995,'IVT':0,'300W':0,'850T':249,'SLPAnom':lowanom}
        highlims = {'Z500Anom':highanom,'SLP':1020,'IVT':452,'300W':45,'850T':293,'SLPAnom':highanom}
        
        contourstart = {'Z500Anom':-1.5,'SLP':997,'IVT':0,'300W':10,'850T':250,'SLPAnom':-1.8}
        contourint = {'Z500Anom':0.25,'SLP':3,'IVT':50,'300W':4,'850T':3,'SLPAnom':0.25}
        
        cbarstart = {'Z500Anom':-1.6,'SLP':995,'IVT':0,'300W':0,'850T':240,'SLPAnom':-1.6}
        cbarint = {'Z500Anom':0.4,'SLP':5,'IVT':75,'300W':10,'850T':10,'SLPAnom':0.4}     
        
    cbarlabs = {'Z500Anom':r'$\mathbf{\sigma}$','SLP':'hPa','IVT':'kg $\mathregular{m^{-1}}$ $\mathregular{s^{-1}}$','300W':'m/s','850T':'K','SLPAnom':r'$\mathbf{\sigma}$'}
    plottitle = {'Z500Anom':'Z500 Anomaly','SLP':'SLP','IVT':'IVT','300W':'300 hPa Wind','850T':'850 hPa Temperature','SLPAnom':'SLP Anomaly'}
    
#%%
    #MAP DESIRED VARIABLE
    fig = plt.figure(figsize=(5.7,4))
    #convert lat and lon into a 2D array
    lon, lat = np.meshgrid(gridlonreduced,gridlatreduced) 
    #define area threshold for basemap
    area_thresh = 1E4
    #create equidistant cylindrical projection basemap
    map = Basemap(projection='cyl',llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,\
              urcrnrlon=lonmax,resolution='l',area_thresh=area_thresh)
    xi, yi = map(lon,lat)
    #create colormap of MERRA2 data
    colorm = map.pcolor(xi,yi,merrareduced,shading='auto',cmap=colormap[metvar],vmin=lowlims[metvar],vmax=highlims[metvar],zorder=1)
    
    #define border color and thickness
    border_c = '0.4'
    border_w = 0.4
    #create map features
    map.drawcoastlines(color=border_c, linewidth=border_w)
    map.drawstates(color=border_c, linewidth=border_w)
    map.drawcountries(color=border_c, linewidth=border_w)
    gridlinefont = 9
    map.drawparallels(np.arange(20.,71.,20.), labels=[1,0,0,0], fontsize=gridlinefont,color=border_c, linewidth=border_w)
    map.drawmeridians(np.arange(-160.,-109.,20.), labels=[0,0,0,1], fontsize=gridlinefont,color=border_c, linewidth=border_w)
    #define contour color and thickness
    contour_c = '0.1'
    contour_w = 0.7
    #create contour map
    if metvar == '850T':
        c_start, c_stop, c_int = (1250,1635,25)
        contourm = map.contour(xi,yi,merrareduced_contours,levels=np.arange(c_start,c_stop,c_int),colors=contour_c,linewidths=contour_w)
        plt.clabel(contourm,levels=np.arange(c_start,c_stop,c_int*2),fontsize=6,inline_spacing=1,colors='k',zorder=5,manual=False)
        plt.title(f'{plottitle[metvar]} and Z850',fontweight="bold",fontsize=12)
    else:
        contourm = map.contour(xi,yi,merrareduced,levels=np.arange(contourstart[metvar],highlims[metvar]+1,contourint[metvar],dtype=float),colors=contour_c,linewidths=contour_w)
        plt.clabel(contourm,levels=np.arange(contourstart[metvar],highlims[metvar]+1,contourint[metvar]*2,dtype=float),fontsize=6,inline_spacing=1,colors='k',zorder=5,manual=False)
        plt.title(f'{plottitle[metvar]}',fontweight="bold",fontsize=12)
        
    #add yuba shape
    #map.readshapefile(os.path.join(ws_directory,f'{watershed}'), watershed,linewidth=0.8,color='r')
    plt.scatter(-120.9,39.5,color='tomato',edgecolors='r',marker='*',linewidths=0.8,zorder=6)
    
    #cbar_ax = fig.add_axes([0.15,0.04,0.7,0.0216]) #bottom colorbar
    if metvar == '850T':
        cbar = map.colorbar(colorm, location='right', pad="6%",ticks=np.arange(cbarstart[metvar],highlims[metvar]+1,cbarint[metvar]),extend='min')
        #cbar = map.colorbar(colorm, location='bottom', pad="10%",ticks=np.arange(cbarstart[metvar],highlims[metvar]+1,cbarint[metvar]),extend='min')
    else:
        cbar = map.colorbar(colorm, location='right', pad="6%",ticks=np.arange(cbarstart[metvar],highlims[metvar]+1,cbarint[metvar]))
        #cbar = map.colorbar(colorm, location='bottom', pad="10%",ticks=np.arange(cbarstart[metvar],highlims[metvar]+1,cbarint[metvar]))
    cbar.ax.tick_params(labelsize=9)
    cbar.set_label(cbarlabs[metvar],fontsize=9,fontweight='bold')
    
    #SHOW MAP
    save_dir='I:\\Emma\\FIROWatersheds\\Figures\\Composites'
    os.chdir(save_dir)
    plt.savefig(f'{metvar}_{percentile}_Composite.png',bbox_inches='tight',pad_inches=0,dpi=300)
    plt.show()