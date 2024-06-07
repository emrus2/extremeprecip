# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 13:49:49 2023

@author: emrus2

Here we will bin the days assigned to each IVT SOM together and then plot the
composite patterns of other variables, like SLP, Z500Anom, etc. 

For a 9-node SOM

UPDATED 6/12/2023

need to add 850 winds and temperature (and change 500 to 300 hPa anomalies)
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
import matplotlib.patches as patches
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
from mpl_toolkits.basemap import Basemap #installed using 
    #conda install -c anaconda basemap
import scipy.io
from datetime import datetime, timedelta

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
    if metvar == 'Z300Anom':
        lowanom, highanom = (-3.1,2.3)
    elif metvar == 'SLPAnom':
        lowanom, highanom = (-3.7, 1.9)
    elif metvar == '850TAnom':
        lowanom, highanom = (-3.1, 2.1)
    newmap = center_colormap(lowanom, highanom, center=0)
    return(lowanom,highanom,newmap)

#%% IMPORT SOM DATA
# change directory and import SOM data from .mat file
mat_dir='I:\\Emma\\FIROWatersheds\\Data\\SOMs\\SomOutput\\'
os.chdir(mat_dir)
numpatterns = 12
percentile = 90
clusters = 3

# import SOM data
os.chdir(mat_dir)
soms = scipy.io.loadmat(f'IVT_{percentile}_{numpatterns}sompatterns_{clusters}d.mat')
gridlat = np.squeeze(soms['lat'])
gridlon = np.squeeze(soms['lon'])
patterns = soms['pats']
pat_freq = np.squeeze(soms['pat_freq'])
pat_prop = np.squeeze(soms['pat_prop'])
asn_err = np.squeeze(soms['assignment'])
assignment = asn_err[:,0]

# import list of arrays with extreme event dates
extremeevents = list(np.load(f'I:\\Emma\\FIROWatersheds\\Data\\ExtremeEvents_{clusters}d.npy', \
                             allow_pickle=True))
# add 2 previous days for composites (becomes 1x5 arrays in list)
extremeevents_exp = [np.concatenate(([a[0]-timedelta(days=2)],[a[0]-timedelta(days=1)], a)) \
                     for a in extremeevents]
# create array with extreme event arrays and node assignment
extremeassign = list(zip(extremeevents_exp,assignment))

# define lat, lon region of data for plotting
latmin, latmax = (15.5,65.5)
lonmin, lonmax = (-170.25,-105.75)

#%% LOAD IN DENORMALIZED DATA
denordata = np.load(os.path.join(mat_dir,f'IVT_{percentile}_{numpatterns}sompatterns_{clusters}d_denormalized.npy'))

#%%
# change directory
#%% REDUCE NODE PATTERN
for node in range(3,numpatterns+1):
    os.chdir('I:\\Emma\\FIROWatersheds\\Data\\DailyMERRA2') 
# for node in range(1,2):
    # node = 1
    som = node - 1 # reduce by 1 for indexing purposes
    # reduce to som pattern for node
    nodepat = denordata[som,:]
    
    #%% OPEN MERRA DATA FOR EACH METVAR
    # define metvars of interest
    metvars =['IVT','300W','Z300Anom','SLP','SLPAnom','850TAnom','850W']
    
    #create empty list to store composite arrays for each var
    allsomcomposites = []
    
    #loop through each metvar for each node
    for metvar in metvars:
        # #define composite location
        if 'Anom' in metvar:
            # folder path is name of var (minus Anom)
            folderpath = metvar[:-4]
        else:
            # folder path is name of metvar
            folderpath = metvar
        filename = f'MERRA2_{folderpath}_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST_5d.nc'
        filepath = os.path.join(folderpath,filename)
        
        #COLLECT VARIABLE DATA FROM MERRA2 FILE
        merravar = {'Z500':'H','SLP':'SLP','850T':'T','Z850':'H'}
        #open the netcdf file in read mode
        gridfile = nc.Dataset(filepath,mode='r')
        # print(gridfile.variables) #date starts at 1980-01-05
        gridlat = gridfile.variables['lat'][:]
        gridlon = gridfile.variables['lon'][:]
        if metvar == 'IVT':
            U = gridfile.variables['UFLXQV'][:]
            V = gridfile.variables['VFLXQV'][:]
            merra = np.sqrt(U**2 + V**2)
        elif 'W' in metvar:
            U = np.squeeze(gridfile.variables['U'][:])
            V = np.squeeze(gridfile.variables['V'][:])
            merra = np.sqrt(U**2 + V**2)
        elif 'Anom' in metvar:
            merra = np.load(os.path.join(folderpath, \
                    f'MERRA2_{metvar}_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST_5d.npy'))
        else:
            merra = gridfile.variables[merravar[metvar]][:]
        time = gridfile.variables['date'][:] # in days since 1980-01-05
        # convert time to datetimes
        time = [datetime.strptime('19800105','%Y%m%d') + timedelta(days=int(j)) for j in time]
        # convert datetimes to strings
        timestr = [datetime.strftime(day,'%Y%m%d') for day in time]
        # remove any size 1 dimension if present
        merra = np.squeeze(merra)
        gridfile.close() #close nc file
            
        # convert from Pa to hPa
        if metvar == 'SLP':
            merra = merra/100 
        
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
        merrareduced = merra[:,latind,:]
        merrareduced = merrareduced[:,:,lonind]
            
        #%% CALCULATE COMPOSITES FOR METVARS
        if metvar == 'IVT' or 'W' in metvar:
            if metvar =='IVT':
                # only need days -4 and -3 composites (will use SOM output for other days)
                endind = clusters-1
            else:
                # other vars need all 5 days of composites
                endind = clusters+2
            # reduce u and v components
            Ureduced = U[:,latind,:]
            Ureduced = Ureduced[:,:,lonind]
            Vreduced = V[:,latind,:]
            Vreduced = Vreduced[:,:,lonind]
            
            # identify days assigned to som
            # arrays to store composites for each SOM 
            U_composites = np.zeros((endind,len(gridlatreduced),len(gridlonreduced)))
            V_composites = np.zeros((endind,len(gridlatreduced),len(gridlonreduced)))
        
            # array of dates for all events for the node assignment
            somdates = np.array([arr[0] for i,arr in enumerate(extremeassign) \
                                 if assignment[i] == int(som + 1)])
            
            # loop through each day of the 5 day som and calculate average
            for dates in range(endind):
                # reduce to date of 5-day som
                somday = somdates[:,dates]
                # convert to strings
                somdaystr = [datetime.strftime(day,'%Y%m%d') for day in somday]
                # find index in time variable for som days
                somindx = [i for i,day in enumerate(timestr) if day in somdaystr]
                # find merra data for those indices
                UNode = [arr for i,arr in enumerate(Ureduced) if i in somindx]
                VNode = [arr for i,arr in enumerate(Vreduced) if i in somindx]
                # calculate average of merra values
                UMean = np.squeeze(np.mean(UNode,axis=0))
                VMean = np.squeeze(np.mean(VNode,axis=0))
                # add to list of node values
                U_composites[dates]= UMean
                V_composites[dates]= VMean
            # calculate the average vector intensity from avg U and V component values
            som_composites = np.sqrt(U_composites**2 + V_composites**2)
            # if IVT, add node pattern as the data for the last 3 days
            if metvar == 'IVT':
                # reshape nodepattern into correct array shape
                nodepat = nodepat.reshape(clusters,len(gridlatreduced),len(gridlonreduced))
                # add nodepat to the end of som_composites array (should have len of 5)
                som_composites = np.concatenate((som_composites,nodepat),axis=0)
            # append array of average variable 5-day data to empty list
            allsomcomposites.append(som_composites)
        
        else:
            # arrays to store composites for each SOM 
            som_composites = np.zeros((clusters+2,len(gridlatreduced),len(gridlonreduced)))
            somdates = np.array([arr[0] for i,arr in enumerate(extremeassign) \
                                 if assignment[i] == int(som + 1)])            
            
            # loop through each day of the 5 day som and calculate average
            for dates in range(clusters+2):
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
    # for i, arr in enumerate(allsomcomposites):
    #     zmax = 0
    #     zmin = 1E8
    #     #determine zmax and zmin for all days
    #     highlim = np.nanmax(arr)
    #     lowlim = np.nanmin(arr)
    #     #print(lowlim,highlim)
    #     if highlim > zmax:
    #         zmax = highlim
    #     if lowlim < zmin:
    #         zmin = lowlim
    
    #     print(f'Variable:{metvars[i]} \nLowest Value:{zmin} \nHighest Value:{zmax}')
    
    #%% DEFINE PLOTTING VARIABLES
      
    lowlims = {'IVT': 0, \
               '300W': 0, \
               'Z300Anom': -3.1, \
               'SLP': 974, \
               'SLPAnom':anom_cm('SLPAnom')[0], \
               '850TAnom': anom_cm('850TAnom')[0], \
               '850W': 0 }
        
    highlims = {'IVT': 875, \
               '300W': 76, \
               'Z300Anom': 2.3, \
               'SLP': 1035, \
               'SLPAnom':anom_cm('SLPAnom')[1], \
               '850TAnom': anom_cm('850TAnom')[1], \
               '850W': 24 }
    
    contourstart = {'IVT': 0, \
               '300W': False, \
               'Z300Anom': -2.7, \
               'SLP': False, \
               'SLPAnom':-3.6, \
               '850TAnom': False, \
               '850W': False }
        
    contourint = {'IVT': 100, \
               '300W': False, \
               'Z300Anom': 0.3, \
               'SLP': False, \
               'SLPAnom':0.3, \
               '850TAnom': False, \
               '850W': False }
    
    cbarstart = {'IVT': 0, \
               '300W': 0, \
               'Z300Anom': False, \
               'SLP': 975, \
               'SLPAnom': False, \
               '850TAnom': -3, \
               '850W': False }
        
    cbarint = {'IVT': 250, \
               '300W': 25, \
               'Z300Anom': False, \
               'SLP': 20, \
               'SLPAnom': False, \
               '850TAnom': 1.5, \
               '850W': False }
    
    colormap = {'IVT': 'gnuplot2_r', \
               '300W': 'hot_r', \
               'Z300Anom': False , \
               'SLP': 'rainbow', \
               'SLPAnom': False, \
               '850TAnom': anom_cm('850TAnom')[2], \
               '850W': False }
    
    cbarlabs = {'IVT': 'kg/m/s', \
               '300W': 'm/s', \
               'Z300Anom': False , \
               'SLP': 'hPa', \
               'SLPAnom': False, \
               '850TAnom': r' ${\sigma}$', \
               '850W': False }
  
    plottitle = {'IVT':'a) IVT', \
                 '300W':'b) 300 hPa Wind Speed and Z300 Anomalies', \
                 'Z300Anom':'b) 300 hPa Wind Speed and Z300 Anomalies', \
                 'SLP':'c) Sea Level Pressure and Anomalies', \
                 '850TAnom':'d) 850 hPa Temperature Anomalies and 850 hPa Wind'}
        
    cbar_x, cbar_len, cbar_wid = (0.928,0.21,0.014)
    
    cbar_y = {'IVT': 0.757, \
               '300W': 0.511, \
               'Z300Anom': False , \
               'SLP': 0.266, \
               'SLPAnom': False, \
               '850TAnom': 0.022, \
               '850W': False }
        
    cbarpad = {'IVT': 6, \
               '300W': 12, \
               'Z300Anom': False , \
               'SLP': 1, \
               'SLPAnom': False, \
               '850TAnom': 2.75, \
               '850W': False }
    
    #%% CREATE FIGURE
    #create subplot for mapping multiple timesteps
    fig = plt.figure(figsize=(10,6.2)) #width, height

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
    ivals = [0,1,3,5]
    placecounter = 0
    for i in ivals:
        sommaps = allsomcomposites[i]
        metvar = metvars[i]
        place = 1
        # grab next som for plotting if multiple
        # variables plotted on same map
        if i in range(1,6,2):
            sommaps2 = allsomcomposites[i+1]
            metvar2 = metvars[i+1]
        # loop through each SOM day
        for j,arr in enumerate(sommaps):
            # grab next som patterns for plotting
            if i in range(1,6,2):
                arr2 = sommaps2[j]
            # add subplot
            plotloc = 1 + (5*placecounter) + j
            print(plotloc)
            ax = fig.add_subplot(4,clusters+2,plotloc)
            # add variable title
            if place == 1:
                    ax.set_title(plottitle[metvar],fontsize=10,pad=3,loc='left',y=1.0)
            sublabel_loc = mtransforms.ScaledTranslation(4/72, -4/72, fig.dpi_scale_trans)
            # add day labels
            if i == 0:
                ax.text(x=0.468, y=1.188, s=f'Day {place-5}', transform=ax.transAxes + sublabel_loc,
                        fontsize=10, fontweight='bold', color = 'k',
                        verticalalignment='top', horizontalalignment='center',
                        bbox=dict(facecolor='none', edgecolor='none', pad=1.5),zorder=3)


    
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
            if place == 1:
                map.drawparallels(parallels, labels=[1,0,0,0], fontsize=gridlinefont, \
                                  color=border_c,linewidth=border_w)
                map.drawmeridians(meridians,color=border_c,linewidth=border_w)
                if i == len(allsomcomposites) - 2:
                    map.drawparallels(parallels, labels=[1,0,0,0], \
                        fontsize=gridlinefont,color=border_c,linewidth=border_w)
                    map.drawmeridians(meridians, labels=[0,0,0,1],  \
                        fontsize=gridlinefont,color=border_c,linewidth=border_w)
            elif i == len(allsomcomposites) - 2:
                map.drawparallels(parallels, color=border_c,linewidth=border_w)
                map.drawmeridians(meridians, labels=[0,0,0,1], \
                   fontsize=gridlinefont,color=border_c,linewidth=border_w)
            else:
                map.drawparallels(parallels, color=border_c,linewidth=border_w)
                map.drawmeridians(meridians,color=border_c,linewidth=border_w)
    
            #define contour color and thickness
            contour_c = '0.1'
            contour_w = 0.7
            #create contour map for specific numbers 
            # (IVT regular)
            if i == 0:
                contourm = map.contour(xi,yi,arr,colors=contour_c,linewidths=contour_w, \
                        levels=np.arange(contourstart[metvar],highlims[metvar]+1, \
                                         contourint[metvar]),zorder=2)
                plt.clabel(contourm,levels=contourm.levels[::2],fontsize=6, \
                           inline_spacing=1,colors='k',zorder=2,manual=False)
            # (Z300 Anoms and SLP Anoms)
            elif i in range(1,4,2):
                contourm = map.contour(xi,yi,arr2,colors=contour_c,linewidths=contour_w, \
                        levels=np.arange(contourstart[metvar2],highlims[metvar2]+1, \
                                         contourint[metvar2]),zorder=2)
                plt.clabel(contourm,levels=contourm.levels[::2],fontsize=6, \
                           inline_spacing=1,colors='k',zorder=2,manual=False)
            # (850 Winds - Vector Map)
            elif i == 5:
                U_arrs = U_composites[j,:,:]
                V_arrs = V_composites[j,:,:]
                interval = 6
                size = 140
                skip = (slice(None, None, interval), slice(None, None, interval))
                vectorm = map.quiver(xi[skip],yi[skip],U_arrs[skip],V_arrs[skip], \
                          pivot='mid',scale=size, scale_units='inches', \
                             color='0.1',width=0.005,alpha=0.8,zorder=5)
            
            #add yuba shape
            plt.scatter(-120.9,39.5,color='w',marker='*',linewidths=0.7,zorder=4)
            place += 1
        #fig.add_axis([left,bottom, width,height])    
        cbar_ax = fig.add_axes([cbar_x,cbar_y[metvar],cbar_wid,cbar_len]) #bottom colorbar
        cbar = fig.colorbar(colorm, cax=cbar_ax,ticks=np.arange(cbarstart[metvar],highlims[metvar]+1, \
                                        cbarint[metvar]),orientation='vertical')
        cbar.ax.tick_params(labelsize=9)
        cbar.set_label(cbarlabs[metvar],fontsize=9.5,labelpad=cbarpad[metvar])
           
        placecounter += 1
        
    # add rectangle around IVT SOMs
    xrect, yrect = (0.389,0.7462) # bottom left corner
    widthrect, heightrect = (0.5345,0.225)
    fig.add_artist(plt.Rectangle((xrect, yrect), widthrect, heightrect, edgecolor='red', \
                    linewidth=2, fill=False))
    
    #CUSTOMIZE SUBPLOT SPACING
    fig.subplots_adjust(left=0.035,right=0.92,bottom=0.02, top=0.965,hspace=0.15, wspace=0.05) #bottom colorbar
    
    #SHOW MAP
    save_dir='I:\\Emma\\FIROWatersheds\\Figures\\SOMs\\Composites'
    os.chdir(save_dir)
    plt.savefig(f'{numpatterns}node_{node}Node_{clusters}d_composites.png',dpi=300)
    plt.show()
