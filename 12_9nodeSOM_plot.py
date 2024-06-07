# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 13:49:49 2023

@author: emrus2

Plots SOM results in a 12-node SOM

UPDATED 6/12/2023
"""
#%% IMPORT MODULES
#import netCDF4 as nc #if this doesn't work, try to reinstall in anaconda prompt using
    #conda install netcdf4
import numpy as np
import os
import matplotlib.pyplot as plt
# from matplotlib.colors import ListedColormap
# import matplotlib.cm as cm
import matplotlib.transforms as mtransforms
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
from mpl_toolkits.basemap import Basemap #installed using 
    #conda install -c anaconda basemap
import scipy.io

#define watershed and directory for plotting
#watershed = 'UpperYuba'
#ws_directory = f'I:\\Emma\\FIROWatersheds\\Data\\WatershedShapefiles\\California\\{watershed}\\'

#%% IMPORT SOM DATA

# define metvar
metvar = 'IVT'
numpatterns = 12
percentile = 90

# change directory and import SOM data from .mat file
mat_dir='I:\\Emma\\FIROWatersheds\\Data\\SOMs\\SomOutput'
os.chdir(mat_dir)
soms = scipy.io.loadmat(f'{metvar}_{percentile}_{numpatterns}_sompatterns.mat')
patterns = soms['pats']
pat_prop = np.squeeze(soms['pat_prop'])
pat_freq = np.squeeze(soms['pat_freq'])
gridlat = np.squeeze(soms['lat'])
gridlon = np.squeeze(soms['lon'])

# define lat, lon region of data for plotting
latmin, latmax = (15.5,65.5)
lonmin, lonmax = (-170.25,-105.75)


#%% IMPORT AVERAGE PRECIP DATA
precip = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_{numpatterns}NodeAssignAvergagePrecip.npy')
precip_rounded = np.round(a=precip,decimals=1)

#%% DETERMINE MAX AND MIN VALIUES
zmax = 0
zmin = 1E8

for i, arr in enumerate(patterns):
    merrareduced = arr.reshape(len(gridlat),len(gridlon))
    
    #determine zmax and zmin for all days
    highlim = np.amax(arr)
    lowlim = np.amin(arr)
    #print(lowlim,highlim)
    if highlim > zmax:
        zmax = highlim
    if lowlim < zmin:
        zmin = lowlim

print(f'Lowest Value:{zmin} \nHighest Value:{zmax}')

#%%
lowlims, highlims = (0,701)
contourstart, contourint = (0,75)
cbarstart, cbarint = (0,100)
colormap = 'gnuplot2_r'
cbarlabs = 'kg $\mathregular{m^{-1}}$ $\mathregular{s^{-1}}$'
plottitle = metvar

#%% PLOT NODES from MATLAB
if numpatterns == 9:
    #create subplot for mapping multiple timesteps
    fig = plt.figure(figsize=(7.5,5.5))
    # fig.suptitle(f'{plottitle[metvar]} SOMs',fontsize=13,fontweight="bold",y=0.9875)
    
    for i, arr in enumerate(patterns):
        # reshape merra data for plotting
        merrareduced = arr.reshape(len(gridlat),len(gridlon))
        # define percentage of node assignment
        perc_assigned = round(pat_prop[i]*100,1)
        # define average precip
        precipavg = precip_rounded[i]
        
        #MAP DESIRED VARIABLE
        #convert lat and lon into a 2D array
        lon, lat = np.meshgrid(gridlon,gridlat) 
        #define area threshold for basemap
        area_thresh = 1E4
        #create equidistant cylindrical projection basemap
        map = Basemap(projection='cyl',llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,\
                  urcrnrlon=lonmax,resolution='l',area_thresh=area_thresh)
        xi, yi = map(lon,lat)
        ax = fig.add_subplot(3,3,i+1)
        sublabel_loc = mtransforms.ScaledTranslation(4/72, -4/72, fig.dpi_scale_trans)
        ax.text(x=0.0, y=1.0, s=i+1, transform=ax.transAxes + sublabel_loc,
            fontsize=9, fontweight='bold', verticalalignment='top', 
            bbox=dict(facecolor='1', edgecolor='none', pad=1.5),zorder=3)
        ax.text(x=0.75, y=1.0, s=f'{perc_assigned}%', transform=ax.transAxes + sublabel_loc,
            fontsize=9, fontweight='bold', verticalalignment='top', color = 'red',
            bbox=dict(facecolor='1', edgecolor='none', pad=1.5),zorder=3)
        ax.text(x=0.4, y=1.0, s=precipavg, transform=ax.transAxes + sublabel_loc,
            fontsize=9, fontweight='bold', verticalalignment='top', color = 'blue',
            bbox=dict(facecolor='1', edgecolor='none', pad=1.5),zorder=3)
    
    
        #create colormap of MERRA2 data
        colorm = map.pcolor(xi,yi,merrareduced,shading='auto',cmap=colormap[metvar],vmin=lowlims[metvar],vmax=highlims[metvar],zorder=1)
        
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
        contourm = map.contour(xi,yi,merrareduced,colors=contour_c,linewidths=contour_w,levels=np.arange(contourstart[metvar],highlims[metvar]+1,contourint[metvar]),zorder=2)
        plt.clabel(contourm,levels=np.arange(contourstart[metvar],highlims[metvar]+1,contourint[metvar]*2),fontsize=6,inline_spacing=1,colors='k',zorder=2,manual=False)
            
        #add yuba shape
        #map.readshapefile(os.path.join(ws_directory,f'{watershed}'), watershed,linewidth=0.8,color='r')
        plt.scatter(-120.9,39.5,color='tomato',edgecolors='r',marker='*',linewidths=0.8,zorder=4)
        
    #CUSTOMIZE SUBPLOT SPACING
    fig.subplots_adjust(left=0.05,right=0.89,bottom=0.021, top=0.955,hspace=0.05, wspace=0.05) #bottom colorbar
    #fig.add_axis([left,bottom, width,height])
    cbar_ax = fig.add_axes([0.904,0.05,0.025,0.88]) #bottom colorbar
    cbar = fig.colorbar(colorm, cax=cbar_ax,ticks=np.arange(cbarstart[metvar],highlims[metvar]+1,cbarint[metvar]),orientation='vertical')
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label(cbarlabs[metvar],fontsize=8.5,labelpad=0.5,fontweight='bold')
    
        
    #SHOW MAP
    save_dir='I:\\{your}'
    os.chdir(save_dir)
    # plt.savefig(f'{metvar}_{percentile}_{numpatterns}SOM.png',dpi=300)
    plt.show()
    
#%%
if numpatterns == 12:
    #create subplot for mapping multiple timesteps
    fig = plt.figure(figsize=(7.2,6.7))
    #fig.suptitle(f'{plottitle} SOMs',fontsize=13,fontweight="bold",y=0.9875)

    for i, arr in enumerate(patterns):
        # reshape merra data for plotting
        merrareduced = arr.reshape(len(gridlat),len(gridlon))
        # define percentage of node assignment
        perc_assigned = str(round(pat_prop[i]*100,1))
        # define average precip
        precipavg = precip_rounded[i]
        # define precip colors
        precipavgnorm = 0.8-(((precipavg-min(precip_rounded))/(max(precip_rounded)+7-min(precip_rounded))))
        precip_col = (precipavgnorm,1.0,precipavgnorm)
        # define frequency colors
        freqnorm = 0.8-(((pat_freq[i]-min(pat_freq))/(max(pat_freq)+7-min(pat_freq))))
        patfreq_col = (1,freqnorm,freqnorm)
        #MAP DESIRED VARIABLE
        #convert lat and lon into a 2D array
        lon, lat = np.meshgrid(gridlon,gridlat) 
        #define area threshold for basemap
        area_thresh = 1E4
        #create equidistant cylindrical projection basemap
        map = Basemap(projection='cyl',llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,\
                  urcrnrlon=lonmax,resolution='l',area_thresh=area_thresh)
        xi, yi = map(lon,lat)
        ax = fig.add_subplot(4,3,i+1)
        sublabel_loc = mtransforms.ScaledTranslation(4/72, -4/72, fig.dpi_scale_trans)
        ax.text(x=0.0, y=1.0, s=i+1, transform=ax.transAxes + sublabel_loc,
            fontsize=10, fontweight='bold', verticalalignment='top', 
            bbox=dict(facecolor='1', edgecolor='none', pad=1.5),zorder=3)
        if len(perc_assigned) == 4:
            xloc = 0.73
        else:
            xloc = 0.77
        ax.text(x=0.37, y=1.0, s=precipavg, transform=ax.transAxes + sublabel_loc,
            fontsize=9, fontweight='bold', verticalalignment='top', color = 'k',
            bbox=dict(facecolor=precip_col, edgecolor='none', pad=1.5),zorder=3)
        ax.text(x=xloc, y=1.0, s=f'{perc_assigned}%', transform=ax.transAxes + sublabel_loc,
            fontsize=9, fontweight='bold',verticalalignment='top', color = 'k',
            bbox=dict(facecolor=patfreq_col, edgecolor='none', pad=1.5),zorder=3)

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
        contourm = map.contour(xi,yi,merrareduced,colors=contour_c,linewidths=contour_w,levels=np.arange(contourstart,highlims+1,contourint),zorder=2)
        plt.clabel(contourm,levels=np.arange(contourstart,highlims+1,contourint*2),fontsize=6,inline_spacing=1,colors='k',zorder=2,manual=False)
            
        #add yuba shape
        plt.scatter(-120.9,39.5,color='w',marker='*',linewidths=0.7,zorder=4)

        
    #CUSTOMIZE SUBPLOT SPACING
    fig.subplots_adjust(left=0.05,right=0.9,bottom=0.026, top=0.985,hspace=0.05, wspace=0.05) #bottom colorbar
    #fig.add_axis([left,bottom, width,height])
    cbar_ax = fig.add_axes([0.91,0.05,0.025,0.9]) #bottom colorbar
    cbar = fig.colorbar(colorm, cax=cbar_ax,ticks=np.arange(cbarstart,highlims+1,cbarint),orientation='vertical')
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label(cbarlabs,fontsize=8.5,labelpad=0.5,fontweight='bold')

        
    #SHOW MAP
    save_dir='I:\\Emma\\FIROWatersheds\\Figures\\SOMs'
    os.chdir(save_dir)
    plt.savefig(f'{metvar}_{percentile}_{numpatterns}SOM.png',dpi=300)
    plt.show()