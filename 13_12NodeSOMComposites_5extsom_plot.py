# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 13:49:49 2023

@author: emrus2

Here we will bin the days assigned to each IVT SOM together and then plot the
composite patterns of other variables, like SLP, Z500Anom, etc. 

For a 9-node SOM

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
from datetime import datetime, timedelta

#%% IMPORT SOM DATA
# change directory and import SOM data from .mat file
mat_dir='I:\\Emma\\FIROWatersheds\\Data'
os.chdir(mat_dir)
numpatterns = 12
percentile = 90
daysprior = 4
clusters = daysprior + 1

os.chdir(mat_dir)
soms = scipy.io.loadmat(f'SOMs\\SomOutput\\IVT_{percentile}_{numpatterns}sompatterns_{clusters}d.mat')
gridlat = np.squeeze(soms['lat'])
gridlon = np.squeeze(soms['lon'])
patterns = soms['pats']
pat_freq = np.squeeze(soms['pat_freq'])
pat_prop = np.squeeze(soms['pat_prop'])
asn_err = np.squeeze(soms['assignment'])
assignment = asn_err[:,0]

extremeevents = list(np.load(f'I:\\Emma\\FIROWatersheds\\Data\\ExtremeEvents_{clusters}d.npy', \
                             allow_pickle=True))
extremeassign = list(zip(extremeevents,assignment))

# define lat, lon region of data for plotting
latmin, latmax = (15.5,65.5)
lonmin, lonmax = (-170.25,-105.75)

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
# define metvar
# metvars = ['IVT','300W','Z500Anom','SLP','SLPAnom','Z850','850W','850TAnom']
metvars=['Z300Anom']
for metvar in metvars:
    os.chdir('I:\\Emma\\FIROWatersheds\\Data\\DailyMERRA2')
    #metvar = '300W'
    
    # #define composite location
    # #metvar = input('Enter MERRA-2 Variable: Z500, SLP, 850T, 300W, or IVT \n')
    if 'Anom' in metvar:
        folderpath = metvar[:-4]
        filename = f'MERRA2_{folderpath}_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST_{clusters}d.nc'
    else:
        folderpath = metvar
        filename = f'MERRA2_{metvar}_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST_{clusters}d.nc'
    filepath = os.path.join(folderpath,filename)
    
    
    #COLLECT VARIABLE DATA FROM MERRA2 FILE
    merravar = {'Z500':'H','SLP':'SLP','850T':'T','Z850':'H', \
                '850TADV':'__xarray_dataarray_variable__' ,\
                'Z300Anom':'H'}
    #open the netcdf file in read mode
    gridfile = nc.Dataset(filepath,mode='r')
    print(gridfile)
    gridlat = gridfile.variables['lat'][:]
    gridlon = gridfile.variables['lon'][:]
    if metvar == 'IVT':
        U = gridfile.variables['UFLXQV'][:]
        V = gridfile.variables['VFLXQV'][:]
        merra = np.sqrt(U**2 + V**2)
    elif 'W' in metvar or metvar == '850QVECT':
        U = np.squeeze(gridfile.variables['U'][:])
        V = np.squeeze(gridfile.variables['V'][:])
        merra = np.sqrt(U**2 + V**2)
        time = gridfile.variables['date'][:]
    elif 'Anom' in metvar:
        merra = np.load(os.path.join(folderpath, \
                f'MERRA2_{metvar}_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST_{clusters}d.npy'))
        time = gridfile.variables['date'][:]
    else:
        merra = gridfile.variables[merravar[metvar]][:]
    time = gridfile.variables['date'][:]
    time = [datetime.strptime('19800105','%Y%m%d') + timedelta(days=int(j)) for j in time]
    timestr = [datetime.strftime(day,'%Y%m%d') for day in time]
    
    
    #date = [datetime.strptime('198001090030','%Y%m%d%H%M') + timedelta(minutes = t) for t in time]
    merra = np.squeeze(merra)
    gridfile.close()
        
    if metvar == 'SLP':
        merra = merra/100
        
    # calculate Q vector convergence
    if metvar == '850QVECT':
        Ugrad = np.ufunc.reduce(np.add,np.gradient(U))
        Vgrad = np.ufunc.reduce(np.add,np.gradient(V))
        merra = Ugrad + Vgrad
        # merra = np.ufunc.reduce(np.add,np.gradient(U)) + np.ufunc.reduce(np.add,np.gradient(V))
        merra *= 10**15
        # redefine as merra
        
    
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
    
    #%% INCLUDE CALCULATION OF IVT VECTORS
    
    allsomcomposites = []
    
    if metvar == 'IVT' or 'W' in metvar or metvar == '850QVECT':
        # reduce u and v components
        Ureduced = U[:,latind,:]
        Ureduced = Ureduced[:,:,lonind]
        Vreduced = V[:,latind,:]
        Vreduced = Vreduced[:,:,lonind]
        
        # identify days assigned to som
        for som in range(len(patterns)):
            # arrays to store composites for each SOM 
            U_composites = np.zeros((clusters,len(gridlatreduced),len(gridlonreduced)))
            V_composites = np.zeros((clusters,len(gridlatreduced),len(gridlonreduced)))
        
            somdates = np.array([arr[0] for i,arr in enumerate(extremeassign) \
                                 if assignment[i] == int(som + 1)])
            
            # loop through each day of the 5 day som and calculate average
            for dates in range(clusters):
                # reduce to date of 5-day som
                somday = somdates[:,dates]
                # convert to strings
                somdaystr = [datetime.strftime(day,'%Y%m%d') for day in somday]
                # find index in time variable for som days
                somindx = [i for i,day in enumerate(timestr) if day in somdaystr]
                # find merra data for those indices
                UNode = [arr for i,arr in enumerate(Ureduced) if i in somindx]
                # calculate average of merra values
                UMean = np.squeeze(np.mean(UNode,axis=0))
                U_composites[dates]= UMean
                    
                VNode = [arr for i,arr in enumerate(Vreduced) if i in somindx]
                VMean = np.squeeze(np.mean(VNode,axis=0))
                V_composites[dates]= VMean
                
            # if metvar != '850QVECT':
            som_composites = np.sqrt(U_composites**2 + V_composites**2)
            allsomcomposites.append(som_composites)
    
    else:
        for som in range(len(patterns)):
            # arrays to store composites for each SOM 
            som_composites = np.zeros((clusters,len(gridlatreduced),len(gridlonreduced)))
            somdates = np.array([arr[0] for i,arr in enumerate(extremeassign) \
                                 if assignment[i] == int(som + 1)])
            
            # loop through each day of the 5 day som and calculate average
            for dates in range(clusters):
                # reduce to date of 5-day som
                somday = somdates[:,dates]
                # convert to strings
                somdaystr = [datetime.strftime(day,'%Y%m%d') for day in somday]
                # find index in time variable for som days
                somindx = [i for i,day in enumerate(timestr) if day in somdaystr]
                # find merra data for those indices
                MerraNode = [arr for i,arr in enumerate(merrareduced) if i in somindx]
                # calculate average of merra values
                MerraMean = np.squeeze(np.mean(MerraNode,axis=0))
                som_composites[dates]= MerraMean
            allsomcomposites.append(som_composites)
            
    #%% DETERMINE MAX AND MIN VALIUES
    zmax = 0
    zmin = 1E8
    
    for i, arr in enumerate(allsomcomposites):
        #determine zmax and zmin for all days
        highlim = np.nanmax(arr)
        lowlim = np.nanmin(arr)
        #print(lowlim,highlim)
        if highlim > zmax:
            zmax = highlim
        if lowlim < zmin:
            zmin = lowlim
    
    print(f'Lowest Value:{zmin} \nHighest Value:{zmax}')
        
    #%% DEFINE PLOTTING VARIABLES
        
    lowlims = {'Z500':2850, \
               'SLP':975,'IVT':0,'300W':0,'850T':244, \
               'Z500Anom':anom_cm('Z500Anom')[0], 'Z850':1114,'SLPAnom':anom_cm('SLPAnom')[0], \
               '850TAnom':anom_cm('850TAnom')[0],'850TADV':0,'850QVECT':0, \
               '850W':0,'Z300':anom_cm('Z300Anom')[0]}
        
    highlims = {'Z500':5700,'SLP':1034,'IVT':763,'300W':75,'850T':293, \
                'Z500Anom':anom_cm('Z500Anom')[1],'Z850':1595,'SLPAnom':anom_cm('SLPAnom')[1], \
                '850TAnom':anom_cm('850TAnom')[1],'850TADV':0,'850QVECT':0, \
                '850W':25,'Z300Anom':anom_cm('Z300Anom')[1]}
    
    contourstart = {'Z500':3000,'SLP':978,'IVT':0,'300W':4,'850T':250, \
                    'Z500Anom':-2.8,'Z850':1120,'SLPAnom':-2.8, \
                    '850TAnom':-2.4,'850TADV':-0.4,'850QVECT':-1000, \
                    '850W':2,'Z300Anom':-2.4}
        
    contourint = {'Z500':200,'SLP':4,'IVT':75,'300W':8,'850T':2.5, \
                  'Z500Anom':0.4,'Z850':40,'SLPAnom':0.4, \
                  '850TAnom':0.3,'850TADV':0.1,'850QVECT':100, \
                  '850W':3,'Z300Anom':0.4}
    
    cbarstart = {'Z500':3000,'SLP':975,'IVT':0,'300W':0,'850T':250, \
                 'Z500Anom':-3,'Z850':1150,'SLPAnom':-3, \
                 '850TAnom':-2.5,'850TADV':-0.4,'850QVECT':-600, \
                 '850W':0,'Z300Anom':-2.5}
        
    cbarint = {'Z500':500,'SLP':10,'IVT':100,'300W':10,'850T':5, \
               'Z500Anom':1,'Z850':50,'SLPAnom':1, \
               '850TAnom':0.5,'850TADV':0.2,'850QVECT':200, \
               '850W':4,'Z300Anom':0.5}
    
    colormap = {'Z500':'jet','SLP':'rainbow','IVT':'gnuplot2_r','300W':'hot_r', \
                '850T':'turbo','Z500Anom':anom_cm('Z500Anom')[2],'Z850':'turbo','SLPAnom':anom_cm('SLPAnom')[2], \
                '850TAnom':anom_cm('850TAnom')[2],'850TADV':'jet','850QVECT':'jet', \
                    '850W':'hot_r','Z300Anom':anom_cm('Z300Anom')[2]}
        
    cbarlabs = {'Z500':'m','SLP':'hPa','IVT':'kg $\mathregular{m^{-1}}$ $\mathregular{s^{-1}}$', \
                '300W':'m/s','850T':'K','Z500Anom':r'$\mathbf{\sigma}$','Z850':'m', \
                    'SLPAnom':r'$\mathbf{\sigma}$','850TAnom':r'$\mathbf{\sigma}$', \
                        '850TADV':u'\N{DEGREE SIGN}C/hr','850QVECT':'m/kgs', \
                        '850W':'m/s','Z300Anom':r'$\mathbf{\sigma}$'}
        
    plottitle = {'Z500':'Z500','SLP':'SLP','IVT':'IVT','300W':'300 hPa Wind', \
                 '850T':'850 hPa Temperature','Z500Anom':'Z500 Anomaly','Z850':'Z850', \
                     'SLPAnom':'SLP Anomaly','850TAnom':'850 hPa Temperature Anomaly', \
                         '850TADV':'Tadv','850W':'850 hPa Wind','Z300Anom':'Z300 Anomaly'}
    #%% PLOT NODES from MATLAB
    
    #create subplot for mapping multiple timesteps
    fig = plt.figure(figsize=(10,19))
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
    # loop through som composites
    
    for i, sommaps in enumerate(allsomcomposites):
        place = 1
    
        # loop through each SOM day
        for j,arr in enumerate(sommaps):
        # add subplot
            plotloc = 1 + (5*i) + j
            ax = fig.add_subplot(numpatterns,clusters,plotloc)
            if i == 0:
                ax.set_title(f'Day {place-5}',fontsize=10,fontweight="bold",pad=1)  
            sublabel_loc = mtransforms.ScaledTranslation(4/72, -4/72, fig.dpi_scale_trans)
            # add som number titles
            if place == 1:
                ax.set_ylabel(f'{i+1}',fontsize=12,fontweight="bold",labelpad=10,rotation=0)
            
            #create colormap of MERRA2 data
            colorm = map.pcolor(xi,yi,arr,shading='auto',cmap=colormap[metvar], \
                                vmin=lowlims[metvar],vmax=highlims[metvar],zorder=1)
        
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
                map.drawparallels(parallels, labels=[0,1,0,0], fontsize=gridlinefont, \
                                  color=border_c,linewidth=border_w)
                map.drawmeridians(meridians,color=border_c,linewidth=border_w)
                if i == numpatterns - 1:
                    map.drawparallels(parallels, labels=[0,1,0,0], \
                        fontsize=gridlinefont,color=border_c,linewidth=border_w)
                    map.drawmeridians(meridians, labels=[0,0,0,1],  \
                        fontsize=gridlinefont,color=border_c,linewidth=border_w)
            elif i == numpatterns - 1:
                map.drawparallels(parallels, color=border_c,linewidth=border_w)
                map.drawmeridians(meridians, labels=[0,0,0,1], \
                   fontsize=gridlinefont,color=border_c,linewidth=border_w)
            else:
                map.drawparallels(parallels, color=border_c,linewidth=border_w)
                map.drawmeridians(meridians,color=border_c,linewidth=border_w)
    
            #define contour color and thickness
            contour_c = '0.1'
            contour_w = 0.7
            #create contour map
            if metvar == '850QVECT':
                U_arrs = U_composites[i,:,:]
                V_arrs = V_composites[i,:,:]
                interval = 2
                size = 1000
                skip = (slice(None, None, interval), slice(None, None, interval))
                vectorm = map.quiver(xi[skip],yi[skip],U_arrs[skip],V_arrs[skip],color='darkgreen')
                vectorm = map.quiver(xi[skip],yi[skip],U_arrs[skip],V_arrs[skip], \
                         pivot='mid',scale=size, scale_units='inches',headlength=3.4, \
                             headwidth=2,color='g',width=0.005)
            elif metvar == '850TADV':
                contourm = map.contour(xi,yi,arr,colors=contour_c,linewidths=contour_w, \
                        levels=np.arange(-0.3,0.4,contourint[metvar]),zorder=2)
                plt.clabel(contourm,levels=contourm.levels[::2],fontsize=6, \
                           inline_spacing=1,colors='k',zorder=2,manual=False)
            else:
                contourm = map.contour(xi,yi,arr,colors=contour_c,linewidths=contour_w, \
                        levels=np.arange(contourstart[metvar],highlims[metvar]+1, \
                                         contourint[metvar]),zorder=2)
                plt.clabel(contourm,levels=contourm.levels[::2],fontsize=6, \
                           inline_spacing=1,colors='k',zorder=2,manual=False)
                
            #add yuba shape
            #map.readshapefile(os.path.join(ws_directory,f'{watershed}'), watershed,linewidth=0.8,color='r')
            plt.scatter(-120.9,39.5,color='w',marker='*',linewidths=0.7,zorder=4)
            place += 1
            
    #CUSTOMIZE SUBPLOT SPACING
    fig.subplots_adjust(left=0.035,right=0.96,bottom=0.046, top=0.991,hspace=0.05, wspace=0.05) #bottom colorbar
    #fig.add_axis([left,bottom, width,height])
    cbar_ax = fig.add_axes([0.05,0.023,0.9,0.014]) #bottom colorbar
    cbar = fig.colorbar(colorm, cax=cbar_ax,ticks=np.arange(cbarstart[metvar],highlims[metvar]+1, \
                                        cbarint[metvar]),orientation='horizontal')
    cbar.ax.tick_params(labelsize=9)
    cbar.set_label(cbarlabs[metvar],fontsize=9.5,labelpad=0.5,fontweight='bold')
    
        
    #SHOW MAP
    save_dir='I:\\Emma\\FIROWatersheds\\Figures\\SOMs\\Composites'
    os.chdir(save_dir)
    plt.savefig(f'{numpatterns}node_{metvar}_{clusters}dSOM_composite_vert.png',dpi=300)
    plt.show()