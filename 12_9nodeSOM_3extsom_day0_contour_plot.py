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
numpatterns = 9
percentile = 90
daysprior = 2
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
#latmin, latmax = (20.5,70.5)
#lonmin, lonmax = (-170.25,-105.75)
#latmin, latmax = (20.5,60.5)
latmin, latmax = (15.5,65.5)
lonmin, lonmax = (-170.25,-105.75)

#%% IMPORT AVERAGE PRECIP DATA
precip = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_{numpatterns}NodeAssignAvergagePrecip.npy')
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
lowlims = {'Z500':2850,'SLP':640,'IVT':0,'300W':0,'850T':246,'Z500Anom':0,'SLPAnom':0}
highlims = {'Z500':5700,'SLP':1000,'IVT':770,'300W':68,'850T':294,'Z500Anom':0,'SLPAnom':0}

contourstart = {'Z500':3000,'SLP':650,'IVT':0,'300W':7,'850T':250,'Z500Anom':-2.75,'SLPAnom':-2.75}
contourint = {'Z500':200,'SLP':50,'IVT':75,'300W':7,'850T':3,'Z500Anom':0.25,'SLPAnom':0.25}

cbarstart = {'Z500':3000,'SLP':650,'IVT':0,'300W':0,'850T':250,'Z500Anom':-2.5,'SLPAnom':-3}
cbarint = {'Z500':500,'SLP':50,'IVT':100,'300W':10,'850T':10,'Z500Anom':0.5,'SLPAnom':0.5}

colormap = {'Z500':'jet','SLP':'rainbow','IVT':'gnuplot2_r','300W':'hot_r','850T':'turbo','Z500Anom':'turbo','SLPAnom':'turbo'}
cbarlabs = {'Z500':'m','SLP':'hPa','IVT':'kg $\mathregular{m^{-1}}$ $\mathregular{s^{-1}}$','300W':'m/s','850T':'K','Z500Anom':r'$\mathbf{\sigma}$','SLPAnom':r'$\mathbf{\sigma}$'}
plottitle = {'Z500':'Z500','SLP':'SLP','IVT':'IVT','300W':'300 hPa Wind','850T':'850 hPa Temperature','Z500Anom':'Z500 Anomaly','SLPAnom':'SLP Anomaly'}


#%% PLOT NODES from MATLAB
contourcolor = {0:'#8F3A13',1:'#49DBCD',2:'0.15'}
contourline = {0:'dashdot',1:'dashed',2:'solid'}

#create subplot for mapping multiple timesteps
fig = plt.figure(figsize=(7.2,5))
#fig.suptitle(f'{plottitle[metvar]} SOMs',fontsize=13,fontweight="bold",y=0.9875)
for i, arr in enumerate(patterns):
    # reshape merra data for plotting    
    merrareduced = arr.reshape(clusters,len(gridlat),len(gridlon))
    # define percentage of node assignment
    perc_assigned = str(round(pat_prop[i]*100,1))
    # define average precip
    #precipavg = precip_rounded[i]
    #patfreq_col = str(1 / (pat_freq[i]/17))
    gbcolor = 1.1-(pat_freq[i]/44)
    patfreq_col = (1,gbcolor,gbcolor)
    # add subplot
    ax = fig.add_subplot(3,3,i+1)
    #convert lat and lon into a 2D array
    lon, lat = np.meshgrid(gridlon,gridlat) 
    #define area threshold for basemap
    area_thresh = 1E4
    #create equidistant cylindrical projection basemap
    map = Basemap(projection='cyl',llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,\
              urcrnrlon=lonmax,resolution='l',area_thresh=area_thresh)
    xi, yi = map(lon,lat)
    sublabel_loc = mtransforms.ScaledTranslation(4/72, -4/72, fig.dpi_scale_trans)
    ax.text(x=0.0, y=1.0, s=i+1, transform=ax.transAxes + sublabel_loc,
        fontsize=10, fontweight='bold', verticalalignment='top', 
        bbox=dict(facecolor='1', edgecolor='none', pad=1.5),zorder=3)
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
    #define contour color and thickness
    contour_c = '0.1'
    contour_w = 0.6

    #MAP DESIRED VARIABLE
    for j in range(len(merrareduced)):
        merraplot = merrareduced[j]
        if j == 2:
            colorm = map.pcolor(xi,yi,merraplot,shading='auto', cmap=colormap[metvar], \
                    vmin=lowlims[metvar],vmax=highlims[metvar],zorder=1)
        #create contour map
        else:
            contourm = map.contour(xi,yi,merraplot,colors=contourcolor[j], \
                    linewidths=contour_w,levels=np.arange(contourstart[metvar], \
                    highlims[metvar]+1,contourint[metvar]),zorder=j+2, linestyles=contourline[j])
        # plt.clabel(contourm,levels=np.arange(contourstart[metvar], \
        #         highlims[metvar]+1,contourint[metvar]*2),fontsize=6, \
        #         inline_spacing=1,colors=contourcolor[j],zorder=2,manual=False)
    if len(perc_assigned) == 4:
        ax.text(x=0.755, y=1.0, s=f'{perc_assigned}%', transform=ax.transAxes + sublabel_loc,
            fontsize=8, fontweight='bold',verticalalignment='top', color = 'k',
            bbox=dict(facecolor=patfreq_col, edgecolor='none', pad=1.5),zorder=3)
    else:
        ax.text(x=0.795, y=1.0, s=f'{perc_assigned}%', transform=ax.transAxes + sublabel_loc,
            fontsize=8, fontweight='bold',verticalalignment='top', color = 'k',
            bbox=dict(facecolor=patfreq_col, edgecolor='none', pad=1.5),zorder=3)
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
            
    #add yuba shape
    plt.scatter(-120.9,39.5,color='w',marker='*',linewidths=0.7,zorder=6)
    
#CUSTOMIZE SUBPLOT SPACING
fig.subplots_adjust(left=0.05,right=0.9,bottom=0.026, top=0.985,hspace=0.05, wspace=0.05) #bottom colorbar
#fig.add_axis([left,bottom, width,height])
cbar_ax = fig.add_axes([0.91,0.05,0.025,0.9]) #bottom colorbar
cbar = fig.colorbar(colorm, cax=cbar_ax,ticks=np.arange(cbarstart[metvar],highlims[metvar]+1,cbarint[metvar]),orientation='vertical')
cbar.ax.tick_params(labelsize=8)
cbar.set_label(cbarlabs[metvar],fontsize=8.5,labelpad=0.5,fontweight='bold')

    
#SHOW MAP
save_dir='I:\\Emma\\FIROWatersheds\\Figures\\SOMs'
os.chdir(save_dir)
# plt.savefig(f'{metvar}_{percentile}_{numpatterns}SOM_{clusters}d_sim.png',dpi=300)
plt.show()