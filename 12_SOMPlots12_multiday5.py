-*- coding: utf-8 -*-
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
import matplotlib.transforms as mtransforms
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
from mpl_toolkits.basemap import Basemap #installed using 
    #conda install -c anaconda basemap
import scipy.io

#define watershed and directory for plotting
watershed = 'UpperYuba'
ws_directory = f'I:\\Emma\\FIROWatersheds\\Data\\WatershedShapefiles\\California\\{watershed}\\'
#%% IMPORT SOM DATA

# define metvar
metvar = 'IVT'
numpatterns = 12
percentile = 90
daysprior = 4
clusters = daysprior + 1

# change directory and import SOM data from .mat file
mat_dir='I:\\Emma\\FIROWatersheds\\Data\\SOMs\\SomOutput'
os.chdir(mat_dir)
soms = scipy.io.loadmat(f'IVT_{percentile}_{numpatterns}sompatterns_{clusters}d.mat')
patterns = soms['pats']
pat_prop = np.squeeze(soms['pat_prop'])
pat_freq = np.squeeze(soms['pat_freq'])
gridlat = np.squeeze(soms['lat'])
gridlon = np.squeeze(soms['lon'])

# define lat, lon region of data for plotting
latmin, latmax = (15.5,65.5)
lonmin, lonmax = (-170.25,-105.75)

#%% IMPORT AVERAGE PRECIP DATA
precip = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_{numpatterns}NodeAssignAvergagePrecip_{clusters}d.npy')
precip_rounded = np.round(a=precip,decimals=1)

#%% DETERMINE MAX AND MIN VALIUES
zmax = 0
zmin = 1E8

for i, arr in enumerate(patterns):
    merrareduced = arr.reshape(clusters,len(gridlat),len(gridlon))
    
    #determine zmax and zmin for all days
    highlim = np.amax(arr)
    lowlim = np.amin(arr)
    #print(lowlim,highlim)
    if highlim > zmax:
        zmax = highlim
    if lowlim < zmin:
        zmin = lowlim

print(f'Lowest Value:{zmin} \nHighest Value:{zmax}')
#%% DEFINE PLOTTING VARIABLES

lowlims, highlims = (0,763)
contourstart, contourint = (0,75)
cbarstart, cbarint = (0,100)
colormap = 'gnuplot2_r'
cbarlabs = 'kg $\mathregular{m^{-1}}$ $\mathregular{s^{-1}}$'
plottitle = metvar

#%% PLOT NODES from MATLAB

#create subplot for mapping multiple timesteps
fig = plt.figure(figsize=(10,19))
#fig.suptitle(f'{plottitle} SOMs',fontsize=13,fontweight="bold",y=0.9875)
for i, arr in enumerate(patterns):
    place = 1
    # reshape merra data for plotting
    merrareduced = arr.reshape(clusters,len(gridlat),len(gridlon))
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
    for j in range(len(merrareduced)):
        merraplot = merrareduced[j]
        #convert lat and lon into a 2D array
        lon, lat = np.meshgrid(gridlon,gridlat) 
        #define area threshold for basemap
        area_thresh = 1E4
        #create equidistant cylindrical projection basemap
        map = Basemap(projection='cyl',llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,\
                  urcrnrlon=lonmax,resolution='l',area_thresh=area_thresh)
        xi, yi = map(lon,lat)
        # plotloc = (i+1)+((place-1)*numpatterns)
        plotloc = 1 + (5*i) + j
        # create subplot of (rows,columns,location)
        ax = fig.add_subplot(numpatterns,clusters,plotloc)
        if i == 0:
            ax.set_title(f'Day {place-5}',fontsize=10,fontweight="bold",pad=1)  
            # ax.set_ylabel(f'Day {place-5}',fontsize=10,fontweight="bold",labelpad=0.3)
        sublabel_loc = mtransforms.ScaledTranslation(4/72, -4/72, fig.dpi_scale_trans)
        if place == 5:
            if len(perc_assigned) == 4:
                xloc = 0.695
            else:
                xloc = 0.745
            ax.text(x=0.005, y=1.0, s=precipavg, transform=ax.transAxes + sublabel_loc,
                fontsize=9, fontweight='bold', verticalalignment='top', color = 'k',
                bbox=dict(facecolor=precip_col, edgecolor='none', pad=1.5),zorder=3)
            ax.text(x=xloc, y=1.0, s=f'{perc_assigned}%', transform=ax.transAxes + sublabel_loc,
                fontsize=9, fontweight='bold',verticalalignment='top', color = 'k',
                bbox=dict(facecolor=patfreq_col, edgecolor='none', pad=1.5),zorder=3)
        if place == 1:
            # ax.set_title(f'{i+1}',fontsize=12,fontweight="bold",pad=2)        
            ax.set_ylabel(f'{i+1}',fontsize=12,fontweight="bold",labelpad=10,rotation=0)

        #create colormap of MERRA2 data
        colorm = map.pcolor(xi,yi,merraplot,shading='auto',cmap=colormap,vmin=lowlims,vmax=highlims,zorder=1)
        
        #define border color and thickness
        border_c = '0.4'
        border_w = 0.4
        #create map features
        map.drawcoastlines(color=border_c, linewidth=border_w)
        map.drawstates(color=border_c, linewidth=border_w)
        map.drawcountries(color=border_c, linewidth=border_w)
        gridlinefont = 9
        parallels = np.arange(20.,71.,20.)
        meridians = np.arange(-160.,-109.,20.)
        if place == clusters:
            map.drawparallels(parallels, labels=[0,1,0,0], fontsize=gridlinefont,color=border_c,linewidth=border_w)
            map.drawmeridians(meridians,color=border_c,linewidth=border_w)
            if i == numpatterns - 1:
                map.drawparallels(parallels, labels=[0,1,0,0], fontsize=gridlinefont,color=border_c,linewidth=border_w)
                map.drawmeridians(meridians, labels=[0,0,0,1], fontsize=gridlinefont,color=border_c,linewidth=border_w)
        elif i == numpatterns - 1:
            map.drawparallels(parallels, color=border_c,linewidth=border_w)
            map.drawmeridians(meridians, labels=[0,0,0,1], fontsize=gridlinefont,color=border_c,linewidth=border_w)
        else:
            map.drawparallels(parallels, color=border_c,linewidth=border_w)
            map.drawmeridians(meridians,color=border_c,linewidth=border_w)
        #define contour color and thickness
        contour_c = '0.1'
        contour_w = 0.7
        #create contour map
        contourm = map.contour(xi,yi,merraplot,colors=contour_c,linewidths=contour_w,levels=np.arange(contourstart,highlims+1,contourint),zorder=2)
        plt.clabel(contourm,levels=np.arange(contourstart,highlims+1,contourint*2),fontsize=6,inline_spacing=1,colors='k',zorder=2,manual=False)
            
        #add yuba shape
        plt.scatter(-120.9,39.5,color='w',marker='*',linewidths=0.7,zorder=4)
        place +=1
    
#CUSTOMIZE SUBPLOT SPACING
# fig.subplots_adjust(left=0.01,right=0.948,bottom=0.02, top=0.978,hspace=0.05, wspace=0.05) #bottom colorbar
# #fig.add_axis([left,bottom, width,height])
# cbar_ax = fig.add_axes([0.965,0.05,0.01,0.9]) #bottom colorbar
fig.subplots_adjust(left=0.035,right=0.96,bottom=0.046, top=0.991,hspace=0.05, wspace=0.05) #bottom colorbar
#fig.add_axis([left,bottom, width,height])
cbar_ax = fig.add_axes([0.05,0.023,0.9,0.014]) #bottom colorbar

cbar = fig.colorbar(colorm, cax=cbar_ax,ticks=np.arange(cbarstart,highlims+1,cbarint),orientation='horizontal')
cbar.ax.tick_params(labelsize=9)
cbar.set_label(cbarlabs,fontsize=9.5,labelpad=0.5,fontweight='bold')

    
#SHOW MAP
save_dir='I:\\Emma\\FIROWatersheds\\Figures\\SOMs'
os.chdir(save_dir)
plt.savefig(f'{metvar}_{percentile}_{numpatterns}SOM_{clusters}d_vert.png',dpi=300)
plt.show()