# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 13:49:49 2023

@author: emrus2

Here we plot 3-day prior composites for each node

For a 9-node SOM

UPDATED 6/12/2023
"""
#%% IMPORT MODULES
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
#import matplotlib.transforms as mtransforms
from matplotlib.colors import ListedColormap
from mpl_toolkits.basemap import Basemap 
import os

#%% DEFINE VARIABLES
# define som pattern for plotting and met vars
for sompattern in range(2,10):
    #metvars = ['IVT','300W','Z500','SLP','SLPAnom]
    metvars = ['IVT','SLP','850TAnom','850W','Z500Anom','300W']
    #metvars = ['850W']
    merravar = {'Z500':'H','Z500Anom':'H','SLP':'SLP','SLPAnom':'SLP','850T':'T', \
                '850TAnom':'T','Z850':'H','IVT':'IVTmag','300W':'Windmag', \
                    '850W':'Windmag'}
    folderpath = (f'I:\\Emma\\FIROWatersheds\\Data\\SOMPreviousDaysComposites\\SOM{sompattern}')
    
    # define lon, lat bounds for plotting
    latmin, latmax = (15.5,65.5)
    lonmin, lonmax = (-170.25,-105.75)
    
    # define colorbar axes
    axwidth = 0.02
    axheight = 0.15
    ax_x = 0.887
    
    #%% LOOP THROUGH AND PLOT FIGURES
    fig = plt.figure(figsize=(7.4,7.5))
    i = 0
    for metvar in metvars:
        #print(metvar)
        for dayprev in range(3,-1,-1):
            if 'Anom' in metvar:
                red_metvar = metvar[:-4]
                filename = f'Node{sompattern}_{dayprev}dayprior_composite{red_metvar}.nc'
            else:
                filename = f'Node{sompattern}_{dayprev}dayprior_composite{metvar}.nc'
            filepath = os.path.join(folderpath,filename) 
            #COLLECT VARIABLE DATA FROM MERRA2 FILE
            #open the netcdf file in read mode
            gridfile = nc.Dataset(filepath,mode='r')
            #print(gridfile)
            gridlat = gridfile.variables['lat'][:]
            gridlon = gridfile.variables['lon'][:]
            merra = gridfile.variables[merravar[metvar]][:]
            gridfile.close()
            merra = np.squeeze(merra)
            
            # convert slp to hPa
            if metvar == 'SLP':
                merra = merra/100
                
            # calculate anomalies
            if 'Anom' in metvar:
                gridfile_mean = nc.Dataset(f'I:\\Emma\\FIROWatersheds\\Data\\WinterClimatologies\\{red_metvar}_Climatology_1980-2021_OCT-MAR.nc',mode='r')
                climmean = gridfile_mean.variables[merravar[metvar]][:]
                gridfile_mean.close()
                gridfile_std = nc.Dataset(f'I:\\Emma\\FIROWatersheds\\Data\\WinterClimatologies\\{red_metvar}_Climatology_Std_1980-2021_OCT-MAR.nc',mode='r')
                climstd = gridfile_mean.variables[merravar[metvar]][:]
                gridfile_std.close()
                merraanom = merra - climmean
                merra = merraanom / climstd
            #%% REDUCE LAT AND LON TO DESIRED AREA
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
            
            #print(np.amin(merrareduced),np.amax(merrareduced))
            
            #%% DETERMINE MAX AND MIN VALIUES
            # zmax = 0
            # zmin = 1E8
            
            # for i, arr in enumerate(som_composites):
                
            #     #determine zmax and zmin for all days
            #     highlim = np.nanmax(arr)
            #     lowlim = np.nanmin(arr)
            #     #print(lowlim,highlim)
            #     if highlim > zmax:
            #         zmax = highlim
            #     if lowlim < zmin:
            #         zmin = lowlim
            
            #print(f'Lowest Value:{zmin} \nHighest Value:{zmax}')
            
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
                lowanom, highanom = (-2.0, 1.15)
                newmap = center_colormap(lowanom, highanom, center=0)
            elif metvar == 'SLPAnom':
                lowanom, highanom = (-2.4, 0.9)
                newmap = center_colormap(lowanom, highanom, center=0)
            else:
                lowanom, highanom = (-1.3, 1.2)
                newmap = center_colormap(lowanom, highanom, center=0)
            lowlims = {'SLP':985,'IVT':0,'300W':0,'850T':252,'Z500Anom':lowanom, \
                       'Z850':1187,'SLPAnom':lowanom,'850TAnom':lowanom, \
                           '850W':0}
            highlims = {'SLP':1025,'IVT':575,'300W':54,'850T':293,'Z500Anom':highanom,\
                        'Z850':1548,'SLPAnom':highanom,'850TAnom':highanom, \
                            '850W':22}
            
            contourstart = {'SLP':990,'IVT':0,'300W':5,'850T':250,'Z500Anom':-1.75, \
                            'Z850':1190,'SLPAnom':-2.25,'850TAnom':-1.2,'850W':2}
            contourint = {'SLP':4,'IVT':75,'300W':5,'850T':2.5,'Z500Anom':0.25, \
                          'Z850':30,'SLPAnom':0.3,'850TAnom':0.15,'850W':3}
            
            cbarstart = {'SLP':990,'IVT':0,'300W':0,'850T':250,'Z500Anom':-2,\
                         'Z850':1200,'SLPAnom':-2.4,'850TAnom':-1,'850W':0}
            cbarint = {'SLP':10,'IVT':100,'300W':10,'850T':5,'Z500Anom':1,'Z850':50,\
                       'SLPAnom':0.4,'850TAnom':0.5,'850W':5}
            
            colormap = {'SLP':'rainbow','IVT':'gnuplot2_r','300W':'hot_r',\
                        '850T':'turbo','Z500Anom':newmap,'Z850':'turbo',\
                            'SLPAnom':newmap,'850TAnom':newmap,'850W':'hot_r'}
            cbarlabs = {'SLP':'SLP\n(hPa)','IVT':'IVT\n(kg $\mathregular{m^{-1}s^{-1}}$)',\
                        '300W':'V300 \n(m/s)','850T':'K','Z500Anom':'Z500 Anoms.\n'+r'($\mathbf{\sigma}$)','Z850':'m',\
                            'SLPAnom':r'$\mathbf{\sigma}$','850TAnom':'850T Anoms.\n'+r'($\mathbf{\sigma}$)',\
                                '850W':'V850 \n(m/s)'}
            # plottitle = {'SLP':'SLP','IVT':'IVT','300W':'300 hPa Wind',\
            #              '850T':'850 hPa Temperature','Z500Anom':'Z500 Anomaly',\
            #                  'Z850':'Z850','SLPAnom':'SLP Anomaly',\
            #                      '850TAnom':'850 hPa Temperature Anomaly', \
            #                          '850W':'850 hPa Wind'}
            cbarspacing = {'SLP':2.5,'IVT':5,'300W':12,'850T':0,'Z500Anom':10,\
                          'Z850':0,'SLPAnom':0,'850TAnom':3,'850W':12}
            #%% PLOT NODES from MATLAB
            
            #create subplot for mapping multiple timesteps
            #fig.suptitle(f'{plottitle[metvar]} Composites',fontsize=13,fontweight="bold",y=0.9875)
            
            #MAP DESIRED VARIABLE
            #convert lat and lon into a 2D array
            lon, lat = np.meshgrid(gridlonreduced,gridlatreduced) 
            #define area threshold for basemap
            area_thresh = 1E4
            #create equidistant cylindrical projection basemap
            map = Basemap(projection='cyl',llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,\
                      urcrnrlon=lonmax,resolution='l',area_thresh=area_thresh)
            xi, yi = map(lon,lat)
            ax = fig.add_subplot(len(metvars),4,i+1)
            if i in range(0,4):
                if i == 3:
                    ax.set_title(f'Day {dayprev}',fontsize=10,fontweight="bold",pad=2)                
                else:
                    ax.set_title(f'Day -{dayprev}',fontsize=10,fontweight="bold",pad=2)
    
            #sublabel_loc = mtransforms.ScaledTranslation(4/72, -4/72, fig.dpi_scale_trans)
            # ax.text(0.0, 1.0, i+1, transform=ax.transAxes + sublabel_loc,
            #     fontsize=10, fontweight='bold', verticalalignment='top', 
            #     bbox=dict(facecolor='1', edgecolor='none', pad=1.5),zorder=3)
            #create colormap of MERRA2 data
            colorm = map.pcolor(xi,yi,merrareduced,shading='auto',cmap=colormap[metvar],vmin=lowlims[metvar],vmax=highlims[metvar],zorder=1)
            #colorm = map.pcolor(xi,yi,merrareduced,shading='auto',cmap=colormap[metvar],zorder=1)
            
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
            if i in range(0,25,4):
                if i == 20:
                    map.drawparallels(parallels, labels=[1,0,0,0], fontsize=gridlinefont,color=border_c,linewidth=border_w)
                    map.drawmeridians(meridians, labels=[0,0,0,1], fontsize=gridlinefont,color=border_c,linewidth=border_w)
                else:
                    map.drawparallels(parallels, labels=[1,0,0,0], fontsize=gridlinefont,color=border_c,linewidth=border_w)
                    map.drawmeridians(meridians,color=border_c,linewidth=border_w)
            elif i in range(21,25):
                map.drawparallels(parallels, color=border_c,linewidth=border_w)
                map.drawmeridians(meridians, labels=[0,0,0,1], fontsize=gridlinefont,color=border_c,linewidth=border_w)
            else:
                map.drawparallels(parallels, color=border_c,linewidth=border_w)
                map.drawmeridians(meridians,color=border_c,linewidth=border_w)
            #define contour color and thickness
            contour_c = '0.1'
            contour_w = 0.7
            #create contour map
    
            contourm = map.contour(xi,yi,merrareduced,colors=contour_c,linewidths=contour_w,levels=np.arange(contourstart[metvar],highlims[metvar]+1,contourint[metvar]),zorder=2)
            #contourm = map.contour(xi,yi,merrareduced,colors=contour_c,linewidths=contour_w,zorder=2)                
            plt.clabel(contourm,levels=contourm.levels[::2],fontsize=6,inline_spacing=1,colors='k',zorder=2,manual=False)
                
            #add yuba shape
            #map.readshapefile(os.path.join(ws_directory,f'{watershed}'), watershed,linewidth=0.8,color='r')
            plt.scatter(-120.9,39.5,color='w',marker='*',linewidths=0.7,zorder=4)
                
            i += 1
            
    ##CUSTOMIZE SUBPLOT SPACING
    #fig.add_axis([left,bottom, width,height])
        ax_y = 0.831 - (0.1614*((i/4)-1))
        cbar_ax = fig.add_axes([ax_x,ax_y,axwidth,axheight]) #bottom colorbar
        cbar = fig.colorbar(colorm, cax=cbar_ax,ticks=np.arange(cbarstart[metvar],highlims[metvar]+1,cbarint[metvar]),orientation='vertical')
        #cbar = fig.colorbar(colorm, cax=cbar_ax,orientation='vertical')
        cbar.ax.tick_params(labelsize=8)
        #cbar_ax.text(1.3,0.5,cbarlabs[metvar],fontsize=8.5,fontweight='bold',rotation=90)
        cbar.set_label(cbarlabs[metvar],fontsize=8.5,labelpad=cbarspacing[metvar],fontweight='bold')
    
    fig.subplots_adjust(left=0.045,right=0.88,bottom=0.022, top=0.981,hspace=0.05, wspace=0.05) #bottom colorbar
    
        
    #SHOW MAP
    save_dir='I:\\Emma\\FIROWatersheds\\Figures\\SOMs\\Composites'
    os.chdir(save_dir)
    plt.savefig(f'Node{sompattern}_prevdays_composites.png',dpi=300)
    plt.show()