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
from datetime import datetime
import glob
#%% IMPORT EXTREME DAYS DATA
data_dir='I:\\Emma\\FIROWatersheds\\Data\\'
os.chdir(data_dir)
extremedays = np.load('ExtremePrecipitationDays_2022_23.npy',allow_pickle=True)


# define lat, lon region of data for plotting
latmin, latmax = (15.5,65.5)
lonmin, lonmax = (-170.25,-105.75)

# define percentile
percentile = 90

#%% CREATE ANOMALY MAP
def center_colormap(lowlim, highlim, center=0):
    dv = max(-lowlim, highlim) * 2
    N = int(256 * dv / (highlim-lowlim))
    bwr = cm.get_cmap('seismic', N)
    newcolors = bwr(np.linspace(0, 1, N))
    beg = int((dv / 2 + lowlim)*N / dv)
    end = N - int((dv / 2 - highlim)*N / dv)
    newmap = ListedColormap(newcolors[beg:end])
    return newmap

# create function to define anomaly max and mins
def anom_cm(metvar):
    if metvar == 'Z500Anom':
        lowanom, highanom = (-3.2, 2.3)
    elif metvar == 'SLPAnom':
        lowanom, highanom = (-3.55, 1.75)
    elif metvar == '850TAnom':
        lowanom, highanom = (-2.7, 1.85)
    elif metvar == '850QVECT':
        lowanom, highanom = (-600, 600)
    elif metvar == 'Z300Anom':
        lowanom, highanom = (-2.6,2.1)
    else:
        lowanom, highanom = (-0.4, 0.6)
    newmap = center_colormap(lowanom, highanom, center=0)
    return(lowanom,highanom,newmap)

#%% IMPORT MERRA2 DATA
# define metvar for plotting
metvar = '850T'
folder = 'I:\\MERRA2\\Daily_and_Subdaily\\'
metpath = {'IVT':'IVT_hourly_new\\','SLP':'Sea_level_pressure_3hourly\\', \
           'Z300':'300_hPa_Geopotential_Height_3hourly', \
           '300W':'East_and_North_wind_components_at_300_hPa\\', \
           '850W': 'East_and_North_wind_components_at_850_hPa\\', \
           '850T':'Temperature_at_850_hPa_3hourly'}
folderpath = os.path.join(folder,metpath[metvar])

#%% DEFINE MERRA2 VARIABLES NAMES
merravar = {'Z300':'H','SLP':'SLP','850T':'T','Z850':'H', \
            'Z300Anom':'H'}
    
# define plotting parameters for each metvar
lowlims = {'Z300':8200, \
           'SLP':964,'IVT':0,'300W':0,'850T':-32, \
           'Z500Anom':anom_cm('Z500Anom')[0], 'Z850':1114,'SLPAnom':anom_cm('SLPAnom')[0], \
           '850TAnom':anom_cm('850TAnom')[0],'850TADV':0,'850QVECT':0, \
           '850W':0,'Z300Anom':anom_cm('Z300Anom')[0]}
    
highlims = {'Z300':9750,'SLP':1045,'IVT':975,'300W':93,'850T':21, \
            'Z500Anom':anom_cm('Z500Anom')[1],'Z850':1595,'SLPAnom':anom_cm('SLPAnom')[1], \
            '850TAnom':anom_cm('850TAnom')[1],'850TADV':0,'850QVECT':0, \
            '850W':40,'Z300Anom':anom_cm('Z300Anom')[1]}

contourstart = {'Z300':8300,'SLP':966,'IVT':0,'300W':5,'850T':-30, \
                'Z500Anom':-2.8,'Z850':1120,'SLPAnom':-2.8, \
                '850TAnom':-2.4,'850TADV':-0.4,'850QVECT':-1000, \
                '850W':2,'Z300Anom':-2.4}
    
contourint = {'Z300':100,'SLP':4,'IVT':100,'300W':10,'850T':2.5, \
              'Z500Anom':0.4,'Z850':40,'SLPAnom':0.4, \
              '850TAnom':0.3,'850TADV':0.1,'850QVECT':100, \
              '850W':5,'Z300Anom':0.4}

cbarstart = {'Z300':8250,'SLP':970,'IVT':0,'300W':0,'850T':-30, \
             'Z500Anom':-3,'Z850':1150,'SLPAnom':-3, \
             '850TAnom':-2.5,'850TADV':-0.4,'850QVECT':-600, \
             '850W':0,'Z300Anom':-2.5}
    
cbarint = {'Z300':250,'SLP':15,'IVT':150,'300W':15,'850T':10, \
           'Z500Anom':1,'Z850':50,'SLPAnom':1, \
           '850TAnom':0.5,'850TADV':0.2,'850QVECT':200, \
           '850W':10,'Z300Anom':0.5}

colormap = {'Z300':'jet','SLP':'rainbow','IVT':'gnuplot2_r','300W':'hot_r', \
            '850T':'turbo','Z500Anom':anom_cm('Z500Anom')[2],'Z850':'turbo','SLPAnom':anom_cm('SLPAnom')[2], \
            '850TAnom':anom_cm('850TAnom')[2],'850TADV':'jet','850QVECT':'jet', \
                '850W':'hot_r','Z300Anom':anom_cm('Z300Anom')[2]}
    
cbarlabs = {'Z300':'m','SLP':'hPa','IVT':'kg/m/s', \
            '300W':'m/s','850T':u'\N{DEGREE SIGN}C','Z500Anom':r'$\mathbf{\sigma}$','Z850':'m', \
                'SLPAnom':r'$\mathbf{\sigma}$','850TAnom':r'$\mathbf{\sigma}$', \
                    '850TADV':u'\N{DEGREE SIGN}C/hr','850QVECT':'m/kgs', \
                    '850W':'m/s','Z300Anom':r'$\mathbf{\sigma}$'}
    
# plottitle = {'Z300':'Z500','SLP':'SLP','IVT':'IVT','300W':'300 hPa Wind', \
#              '850T':'850 hPa Temperature','Z500Anom':'Z500 Anomaly','Z850':'Z850', \
#                  'SLPAnom':'SLP Anomaly','850TAnom':'850 hPa Temperature Anomaly', \
#                      '850TADV':'Tadv','850W':'850 hPa Wind','Z300Anom':'Z300 Anomaly'}

#%%
fig = plt.figure(figsize=(8,4))

for i,arr in enumerate(extremedays):
    date = arr[0]
    dt = datetime.strptime(date,'%Y%m%d')
    precip = arr[1]
    filepath = glob.glob(os.path.join(folderpath,f'*.{date}.SUB.nc'))
    filename = filepath[0]
    
    #COLLECT VARIABLE DATA FROM MERRA2 FILE
    gridfile = nc.Dataset(filename,mode='r')
    # print(gridfile)
    gridlat = gridfile.variables['lat'][:]
    gridlon = gridfile.variables['lon'][:]
    if metvar == 'IVT':
        # collect U and V components
        Uvapor = gridfile.variables['UFLXQV'][:]
        Vvapor = gridfile.variables['VFLXQV'][:]
        # average vars. across timesteps for daily avg.
        Uvaporavg = Uvapor.mean(axis=0) 
        Vvaporavg = Vvapor.mean(axis=0)
        # convert U and V components to vector magnitude
        merra = np.sqrt(Uvaporavg**2 + Vvaporavg**2)
    elif 'W' in metvar or metvar == '850QVECT':
        # collect U and V components
        U = np.squeeze(gridfile.variables['U'][:])
        V = np.squeeze(gridfile.variables['V'][:])
        # convert U and V components to vector magnitude
        merra = np.sqrt(U**2 + V**2)
    # elif 'Anom' in metvar:
    #     merra = np.load(os.path.join(folderpath, \
    #             f'MERRA2_{metvar}_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST_{clusters}d.npy'))
    #     time = gridfile.variables['date'][:]
    else:
        merra = gridfile.variables[merravar[metvar]][:]
        merra = merra.mean(axis=0)
    gridfile.close()
    merra = np.squeeze(merra)
        
    if metvar == 'SLP':
        merra = merra/100 # convert to hPa
    if 'T' in metvar:
        merra = merra - 273.15
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
    precipround = round(float(precip),1)
    if len(str(precipround)) == 5:
        xloc = 0.7
    else:
        xloc = 0.75
    ax.text(xloc, 1.0, precipround, transform=ax.transAxes + sublabel_loc,
        fontsize=9, fontweight='bold', verticalalignment='top', 
        bbox=dict(facecolor=(0,1,0), edgecolor='none', pad=1.5),zorder=3)
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
    if 'T' in metvar:
        plt.rcParams['contour.negative_linestyle'] = 'solid'
    contourm = map.contour(xi,yi,merrareduced,colors=contour_c,linewidths=contour_w,levels=np.arange(contourstart[metvar],highlims[metvar]+1,contourint[metvar]),zorder=2)
    plt.clabel(contourm,levels=contourm.levels[::2],fontsize=6,inline_spacing=1,colors='k',zorder=2,manual=False)
        
    #add yuba shape
    #map.readshapefile(os.path.join(ws_directory,f'{watershed}'), watershed,linewidth=0.8,color='r')
    plt.scatter(-120.9,39.5,color='w',marker='*',linewidths=0.7,zorder=4)

fig.subplots_adjust(left=0.05,right=0.89,bottom=0.01, top=0.97,hspace=0.05, wspace=0.05) #bottom colorbar

#CUSTOMIZE SUBPLOT SPACING
#fig.add_axis([left,bottom, width,height])
cbar_ax = fig.add_axes([0.9,0.05,0.025,0.88]) #bottom colorbar
cbar = fig.colorbar(colorm, cax=cbar_ax,ticks=np.arange(cbarstart[metvar],highlims[metvar]+1,cbarint[metvar]),orientation='vertical')
cbar.ax.tick_params(labelsize=8)
cbar.set_label(cbarlabs[metvar],fontsize=8.5,labelpad=0.5,fontweight='bold')


#SHOW MAP
save_dir='I:\\Emma\\FIROWatersheds\\Figures\\Events'
os.chdir(save_dir)
plt.savefig(f'202223_ExtremeDays_{metvar}.png',dpi=300,bbox_inches='tight')
plt.show()