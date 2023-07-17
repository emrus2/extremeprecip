# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 13:49:49 2023

@author: emrus2

Here we will bin the days assigned to each IVT SOM together and then plot the
composite patterns of other variables, like SLP, Z500Anom, etc. 

For a 12-node SOM

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
import scipy.io
#%% IMPORT SOM DATA
# change directory and import SOM data from .mat file
mat_dir='I:\\Emma\\FIROWatersheds\\Data\\SOMs\\SomOutput'
numpatterns = 12
percentile = 90

os.chdir(mat_dir)
soms = scipy.io.loadmat(f'IVT_{percentile}_{numpatterns}_sompatterns.mat')
gridlat = np.squeeze(soms['lat'])
gridlon = np.squeeze(soms['lon'])
patterns = soms['pats']
pat_freq = np.squeeze(soms['pat_freq'])
pat_prop = np.squeeze(soms['pat_prop'])
asn_err = np.squeeze(soms['assignment'])
assignment = asn_err[:,0]

# define lat, lon region of data for plotting
latmin, latmax = (15.5,65.5)
lonmin, lonmax = (-170.25,-105.75)

#%% IMPORT MERRA2 DATA
# define metvar
metvars = ['SLP', '300W','Z500Anom','SLPAnom','Z850','850T','IVT']
metvars = ['IVT']
for metvar in metvars:
    #define composite location
    #metvar = input('Enter MERRA-2 Variable: Z500, SLP, 850T, 300W, or IVT \n')
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
    
    if metvar == 'IVT':
        Uvaporreduced = Uvapor[:,latind,:]
        Uvaporreduced = Uvaporreduced[:,:,lonind]
        Vvaporreduced = Vvapor[:,latind,:]
        Vvaporreduced = Vvaporreduced[:,:,lonind]
        
    #print(np.amin(merrareduced),np.amax(merrareduced))
    
    #%% CLUSTER ASSIGNED PATTERNS AND CALCULATE AVERAGE

    # create array to store 12 SOM composites
    som_composites = np.zeros((len(patterns),len(gridlatreduced),len(gridlonreduced)))
    U_composites = np.zeros((len(patterns),len(gridlatreduced),len(gridlonreduced)))
    V_composites = np.zeros((len(patterns),len(gridlatreduced),len(gridlonreduced)))
    # loop through all 12 som patterns
    for som in range(len(patterns[:,0])):
        # create array to store assigned days data
        som_merra = np.zeros((1,len(gridlatreduced),len(gridlonreduced)))
        U_merra = np.zeros((1,len(gridlatreduced),len(gridlonreduced)))
        V_merra = np.zeros((1,len(gridlatreduced),len(gridlonreduced)))
        # loop through all days
        for day,arr in enumerate(merrareduced):
            U_array = Uvaporreduced[day,:,:]
            V_array = Vvaporreduced[day,:,:]
            # add data to som_merra if day is assigned to node
            if assignment[day] == som + 1:
                som_merra = np.concatenate((som_merra,np.expand_dims(arr,axis=0)))
                U_merra = np.concatenate((U_merra,np.expand_dims(U_array,axis=0)))
                V_merra = np.concatenate((V_merra,np.expand_dims(V_array,axis=0)))
        # remove initial row of zeros
        som_merra = som_merra[1:,:,:]
        U_merra = U_merra[1:,:,:]
        V_merra = V_merra[1:,:,:]
        # confirm correct number of days assigned to node
        print(som+1,len(som_merra),pat_freq[som])
        # calculate the mean of assigned days
        som_mean = np.squeeze(np.mean(som_merra,axis=0))
        som_mean = np.squeeze(np.mean(som_merra,axis=0))
        som_mean = np.squeeze(np.mean(som_merra,axis=0))
        # append to array of composites
        som_composites[som] = som_mean
        som_mean = np.squeeze(np.mean(som_merra,axis=0))
        som_mean = np.squeeze(np.mean(som_merra,axis=0))
    
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
    
    if metvar == 'Z500Anom':
        lowanom, highanom = (-2.6, 1.4)
        newmap = center_colormap(lowanom, highanom, center=0)
    else:
        lowanom, highanom = (-2.8, 1.2)
        newmap = center_colormap(lowanom, highanom, center=0)
    #%% DEFINE PLOTTING VARIABLES
    if percentile == 95:
        if metvar == 'Z500Anom':
            lowanom, highanom = (-2.6, 1.4)
            newmap = center_colormap(lowanom, highanom, center=0)
        else:
            lowanom, highanom = (-2.8, 1.2)
            newmap = center_colormap(lowanom, highanom, center=0)
            
        lowlims = {'Z500':2850,'SLP':979,'IVT':0,'300W':0,'850T':248,'Z500Anom':lowanom,'Z850':1153,'SLPAnom':lowanom}
        highlims = {'Z500':5700,'SLP':1025,'IVT':800,'300W':65,'850T':293,'Z500Anom':highanom,'Z850':1555,'SLPAnom':highanom}
        
        contourstart = {'Z500':3000,'SLP':980,'IVT':0,'300W':5,'850T':250,'Z500Anom':-2.5,'Z850':1160,'SLPAnom':-2.75}
        contourint = {'Z500':200,'SLP':4,'IVT':75,'300W':5,'850T':2.5,'Z500Anom':0.25,'Z850':30,'SLPAnom':0.25}
        
        cbarstart = {'Z500':3000,'SLP':980,'IVT':0,'300W':0,'850T':250,'Z500Anom':-2.5,'Z850':1175,'SLPAnom':-2.5}
        cbarint = {'Z500':500,'SLP':5,'IVT':100,'300W':10,'850T':5,'Z500Anom':0.5,'Z850':50,'SLPAnom':0.5}

    elif percentile == 90:
        if metvar == 'Z500Anom':
            lowanom, highanom = (-1.9, 1.0)
            newmap = center_colormap(lowanom, highanom, center=0)
        else:
            lowanom, highanom = (-2.6, 1.0)
            newmap = center_colormap(lowanom, highanom, center=0)
            
        lowlims = {'Z500':2850,'SLP':984,'IVT':0,'300W':0,'850T':250,'Z500Anom':lowanom,'Z850':1175,'SLPAnom':lowanom}
        highlims = {'Z500':5700,'SLP':1022,'IVT':800,'300W':60,'850T':293,'Z500Anom':highanom,'Z850':1556,'SLPAnom':highanom}
        
        contourstart = {'Z500':3000,'SLP':985,'IVT':0,'300W':5,'850T':250,'Z500Anom':-1.75,'Z850':1180,'SLPAnom':-2.5}
        contourint = {'Z500':200,'SLP':4,'IVT':75,'300W':5,'850T':2.5,'Z500Anom':0.25,'Z850':30,'SLPAnom':0.25}
        
        cbarstart = {'Z500':3000,'SLP':985,'IVT':0,'300W':0,'850T':250,'Z500Anom':-1.5,'Z850':1175,'SLPAnom':-2.5}
        cbarint = {'Z500':500,'SLP':5,'IVT':100,'300W':10,'850T':5,'Z500Anom':0.5,'Z850':50,'SLPAnom':0.5}

    #create subplot for mapping multiple timesteps
    colormap = {'Z500':'jet','SLP':'rainbow','IVT':'gnuplot2_r','300W':'hot_r','850T':'turbo','Z500Anom':newmap,'Z850':'turbo','SLPAnom':newmap}
    cbarlabs = {'Z500':'m','SLP':'hPa','IVT':'kg $\mathregular{m^{-1}}$ $\mathregular{s^{-1}}$','300W':'m/s','850T':'K','Z500Anom':r'$\mathbf{\sigma}$','Z850':'m','SLPAnom':r'$\mathbf{\sigma}$'}
    plottitle = {'Z500':'Z500','SLP':'SLP','IVT':'IVT','300W':'300 hPa Wind','850T':'850 hPa Temperature','Z500Anom':'Z500 Anomaly','Z850':'Z850','SLPAnom':'SLP Anomaly'}
    
    #%% PLOT NODES from MATLAB
    
    #create subplot for mapping multiple timesteps
    fig = plt.figure(figsize=(7.5,7))
    fig.suptitle(f'{plottitle[metvar]} Composites',fontsize=13,fontweight="bold",y=0.9875)
    
    for i, arr in enumerate(som_composites):
        
        #MAP DESIRED VARIABLE
        #convert lat and lon into a 2D array
        lon, lat = np.meshgrid(gridlonreduced,gridlatreduced) 
        #define area threshold for basemap
        area_thresh = 1E4
        #create equidistant cylindrical projection basemap
        map = Basemap(projection='cyl',llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,\
                  urcrnrlon=lonmax,resolution='l',area_thresh=area_thresh)
        xi, yi = map(lon,lat)
        ax = fig.add_subplot(4,3,i+1)
        sublabel_loc = mtransforms.ScaledTranslation(4/72, -4/72, fig.dpi_scale_trans)
        ax.text(0.0, 1.0, i+1, transform=ax.transAxes + sublabel_loc,
            fontsize=9, fontweight='bold', verticalalignment='top', 
            bbox=dict(facecolor='1', edgecolor='none', pad=1.5),zorder=3)
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
        if i == 0 or i == 3 or i == 6:
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
        if metvar == 'IVT':
            interval = 7
            size = 600
            skip = (slice(None, None, interval), slice(None, None, interval))
            vectorm = map.quiver(xi[skip],yi[skip],ewindnew[skip],nwindnew[skip],pivot='mid',scale=size, scale_units='inches',headlength=3.4,headwidth=2,color='b',width=0.005)
        else:
            contourm = map.contour(xi,yi,arr,colors=contour_c,linewidths=contour_w,levels=np.arange(contourstart[metvar],highlims[metvar]+1,contourint[metvar]),zorder=2)
            plt.clabel(contourm,levels=np.arange(contourstart[metvar],highlims[metvar]+1,contourint[metvar]*2),fontsize=6,inline_spacing=1,colors='k',zorder=2,manual=False)
            
        #add yuba shape
        #map.readshapefile(os.path.join(ws_directory,f'{watershed}'), watershed,linewidth=0.8,color='r')
        plt.scatter(-120.9,39.5,color='tomato',edgecolors='r',marker='*',linewidths=0.8,zorder=4)
        
    #CUSTOMIZE SUBPLOT SPACING
    fig.subplots_adjust(left=0.05,right=0.89,bottom=0.02, top=0.96,hspace=0.05, wspace=0.05) #bottom colorbar
    #fig.add_axis([left,bottom, width,height])
    cbar_ax = fig.add_axes([0.904,0.05,0.025,0.88]) #bottom colorbar
    cbar = fig.colorbar(colorm, cax=cbar_ax,ticks=np.arange(cbarstart[metvar],highlims[metvar]+1,cbarint[metvar]),orientation='vertical')
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label(cbarlabs[metvar],fontsize=8.5,labelpad=0.5,fontweight='bold')
    
        
    #SHOW MAP
    save_dir='I:\\Emma\\FIROWatersheds\\Figures\\SOMs\\Composites'
    os.chdir(save_dir)
    plt.savefig(f'{metvar}_{percentile}_{numpatterns}_SOM_composite.png',dpi=300)
    plt.show()