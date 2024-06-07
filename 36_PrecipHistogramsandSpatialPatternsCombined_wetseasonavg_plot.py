# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 13:49:49 2023

@author: emrus2

Precipitation histograms for each node

UPDATED 7/11/2023
"""
#%% IMPORT MODULES
import netCDF4 as nc #if this doesn't work, try to reinstall in anaconda prompt using
    #conda install netcdf4
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import matplotlib.transforms as mtransforms
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
from mpl_toolkits.basemap import Basemap #installed using 
    #conda install -c anaconda basemap
import scipy.io

#%% IMPORT SOM DATA
# define metvar
metvar = 'IVT'
numpatterns = 12
percentile = 90
clusters = 5

# change directory and import SOM data from .mat file
mat_dir='I:\\Emma\\FIROWatersheds\\Data\\SOMs\\SomOutput'
os.chdir(mat_dir)
soms = scipy.io.loadmat(f'{metvar}_{percentile}_{numpatterns}sompatterns_{clusters}d.mat')
pat_prop = np.squeeze(soms['pat_prop'])
pat_freq = np.squeeze(soms['pat_freq'])
asn_err = np.squeeze(soms['assignment'])
assignment = asn_err[:,0]

#%% IMPORT WATERSHED SHAPEFILE
watershed = 'UpperYuba'
ws_directory = f'I:\\Emma\\FIROWatersheds\\Data\\WatershedShapefiles\\California\\{watershed}\\'

#%% IMPORT AVERAGE PRECIP DATA
avgprecip = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_{numpatterns}NodeAssignAvergagePrecip_{clusters}d.npy')
avgprecip_rounded = np.round(a=avgprecip,decimals=1)

medprecip = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_{numpatterns}NodeAssignMedianPrecip_{clusters}d.npy')
medprecip_rounded = np.round(a=medprecip,decimals=1)

allprecip = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_{numpatterns}NodeAssignAllPrecip_{clusters}d.npy',allow_pickle=True)

#%% DEFINE SPATIAL MAPPING PARAMETERS
IVT = True # plot IVT vectors
# define lat, lon region of data for plotting
latmin, latmax = (37.5,41.5)
lonmin, lonmax = (-124.5,-119.5)
#define area threshold for basemap
area_thresh = 1E4

#%% IMPORT NETCDF DATA
#IMPORT NETCDF DATA
#define NC location
filepath = ('I:\\Emma\\FIROWatersheds\\Data\\Gridmet\\GRIDMET_pr_Yuba_Extremes90_DailyAnomProp_1980-2021_WINTERDIST_wetseason.nc')
#open the netcdf file in read mode
gridfile = nc.Dataset(filepath,mode='r')
print(gridfile)
gridlat = gridfile.variables['lat'][:]
gridlon = gridfile.variables['lon'][:]
precip = np.squeeze(gridfile.variables['precipitation_amount'][:])
days = gridfile.variables['day'][:]
gridfile.close()
#dates = [datetime.strftime(datetime(1900,1,1)+timedelta(days=days[n]),"%Y%m%d") for n in range(days.shape[0])]

#REDUCE VARIABLES TO DESIRED AREA
#reduce lat
latlims = np.logical_and(gridlat > latmin, gridlat < latmax)
latind = np.where(latlims)[0]
latreduced = gridlat[latind]
#reduce lon
lonlims = np.logical_and(gridlon > lonmin, gridlon < lonmax)
lonind = np.where(lonlims)[0]
lonreduced = gridlon[lonind]
#reduce pressure
merrareduced = precip[:,latind,:]
merrareduced = merrareduced[:,:,lonind]

#convert lat and lon into a 2D array
lon, lat = np.meshgrid(lonreduced,latreduced)

#%% CLUSTER ASSIGNED PATTERNS AND CALCULATE AVERAGE

# create array to store 12 SOM composites
som_composites = np.zeros((numpatterns,len(latreduced),len(lonreduced)))

# loop through all 12 som patterns
for som in range(numpatterns):
    # create array to store assigned days data
    som_merra = np.zeros((1,len(latreduced),len(lonreduced)))
    # loop through all days
    for day,arr in enumerate(merrareduced):
        print(np.amax(arr))
        # add data to som_merra if day is assigned to node
        if assignment[day] == float(som + 1):
            som_merra = np.ma.concatenate((som_merra,np.expand_dims(arr,axis=0)))
            print(np.amax(som_merra))
    # remove initial row of zeros
    som_merra = som_merra[1:,:,:]
    # calculate the mean of assigned days
    som_mean = np.squeeze(np.mean(som_merra,axis=0))
    # append to array of composites
    som_composites[som] = som_mean

#%% DETERMINE MAX AND MIN VALIUES
zmax = 0
zmin = 1E8

som_composites[som_composites==0.] = np.nan

for i, arr in enumerate(som_composites):
    
    #determine zmax and zmin for all days
    highlim = np.nanmax(arr)
    lowlim = np.nanmin(arr)
    if highlim > zmax:
        zmax = highlim
    if lowlim < zmin:
        zmin = lowlim

print(f'Lowest Value:{zmin} \nHighest Value:{zmax}')

#%% IMPORT IVT DATA
if IVT == True:
    # IMPORT MERRA2 DATA
    metvar = 'IVT'
    merrapath = f'I:\\Emma\\FIROWatersheds\\Data\\DailyMERRA2\\{metvar}'
    merraname = f'MERRA2_{metvar}_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST.nc'
    merrafile = os.path.join(merrapath,merraname)

    gridmerra = nc.Dataset(merrafile,mode='r')
    print(gridmerra)
    gridlatm = gridmerra.variables['lat'][:]
    gridlonm = gridmerra.variables['lon'][:]
    Uvapor = gridmerra.variables['UFLXQV'][:]
    Vvapor = gridmerra.variables['VFLXQV'][:]
    merra = np.sqrt(Uvapor**2 + Vvapor**2)
    gridmerra.close()

    #INCLUDE CALCULATION OF IVT VECTORS
    latlimsm = np.logical_and(gridlatm > latmin, gridlatm < latmax)
    latindm = np.where(latlimsm)[0]
    latreducedm = gridlatm[latindm]
    #reduce lon
    lonlimsm = np.logical_and(gridlonm > lonmin, gridlonm < lonmax)
    lonindm = np.where(lonlimsm)[0]
    lonreducedm = gridlonm[lonindm]
    #reduce pressure
    Uvaporreduced = Uvapor[:,latindm,:]
    Uvaporreduced = Uvaporreduced[:,:,lonindm]
    Vvaporreduced = Vvapor[:,latindm,:]
    Vvaporreduced = Vvaporreduced[:,:,lonindm]

    lonm, latm = np.meshgrid(lonreducedm,latreducedm)

    # create array to store 12 SOM composites
    U_composites = np.zeros((numpatterns,len(latreducedm),len(lonreducedm)))
    V_composites = np.zeros((numpatterns,len(latreducedm),len(lonreducedm)))
    # loop through all 12 som patterns
    for som in range(numpatterns):
        # create array to store assigned days data
        U_merra = np.zeros((1,len(latreducedm),len(lonreducedm)))
        V_merra = np.zeros((1,len(latreducedm),len(lonreducedm)))
        # loop through all days
        for day,arr in enumerate(merrareduced):
            U_array = Uvaporreduced[day,:,:]
            V_array = Vvaporreduced[day,:,:]
            # add data to som_merra if day is assigned to node
            if assignment[day] == float(som + 1):
                U_merra = np.concatenate((U_merra,np.expand_dims(U_array,axis=0)))
                V_merra = np.concatenate((V_merra,np.expand_dims(V_array,axis=0)))
        # remove initial row of zeros
        U_merra = U_merra[1:,:,:]
        V_merra = V_merra[1:,:,:]
        # confirm correct number of days assigned to node
        #print(som+1,len(U_merra),pat_freq[som])
        # calculate the mean of assigned days
        U_mean = np.squeeze(np.mean(U_merra,axis=0))
        V_mean = np.squeeze(np.mean(V_merra,axis=0))
        # append to array of composites
        U_composites[som] = U_mean
        V_composites[som] = V_mean

#%% PLOTTING PRECIPITATION HISTOGRAMS

## SET PLOT A) PARAMETERS
colors = ('tomato','cornflowerblue','lightgreen','darkorchid','gold','lightblue','plum','mediumseagreen','indianred','royalblue','grey',(.9,0,.9))
binlist = [np.arange(50,240,10)] # define bins
width = 0.8  # the width of the bars

# define number of columns and rows for subplot
numcols = 33
numrows = 9
# set heights for each row
arows = 1
seprows = 0.25
brows = 1
rowheights = [arows,arows,arows,arows,seprows,brows,brows,brows,brows]
# create figure
fig, axs = plt.subplots(numrows, numcols, figsize=(7, 13), height_ratios=rowheights)
# remove all axes so they aren't seen under plots
for i in range(numrows):
    for j in range(numcols):
        axs[i,j].axis('off')

# set plotting details  
ymin = 0     
ymax = 20
yint = 4 
ylabels = np.arange(ymin,ymax-1,yint)

xmin = 45 
xmax = 215
xint = 25  
xlabels = np.arange(xmin+5,xmax,xint)      
        
## PLOT A) DATA
fig.suptitle('a)',fontsize=12,fontweight="bold",y=0.995,x=0.03) # add title
for i,node in enumerate(allprecip):
    locations = [1,2,3,5,6,7,9,10,11,13,14,15]
    #add subplot
    startloc = 11*i + 1
    endloc = 11*(i+1)
    ax = fig.add_subplot(numrows,numcols,(startloc,endloc))
    precipavg = avgprecip_rounded[i]
    precipmed = medprecip_rounded[i]
    #create colors
    precipmednorm = 0.9-(((precipmed-min(medprecip_rounded))/(max(medprecip_rounded)+7-min(medprecip_rounded))))
    precipmed_col = (1,precipmednorm,1)
    precipavgnorm = 0.82-(((precipavg-min(avgprecip_rounded))/(max(avgprecip_rounded)+7-min(avgprecip_rounded))))
    precipavg_col = (precipavgnorm,1,precipavgnorm)
    # define label location and text
    sublabel_loc = mtransforms.ScaledTranslation(4/72, -4/72, fig.dpi_scale_trans)
    ax.text(x=0.0, y=1.0, s=i+1, transform=ax.transAxes + sublabel_loc,
        fontsize=11, fontweight='bold', verticalalignment='top',
        bbox=dict(facecolor='1', alpha = 0.8, edgecolor='none', pad=1.5),zorder=3)
    ax.text(x=0.8, y=1.0, s=precipavg, transform=ax.transAxes + sublabel_loc,
        fontsize=9, fontweight='bold', verticalalignment='top', color = 'k',
        bbox=dict(facecolor=precipavg_col, edgecolor='none', pad=1.5),zorder=3)
    ax.text(x=0.8, y=0.87, s=precipmed, transform=ax.transAxes + sublabel_loc,
        fontsize=9, fontweight='bold', verticalalignment='top', color = 'k',
        bbox=dict(facecolor=precipmed_col, edgecolor='none', pad=1.5),zorder=3)
    # plot data histogram
    freqs = plt.hist(node,color=colors[i],label=i+1,align='mid',bins=binlist[0])
    # plot mean and median line markers
    plt.axvline(x=precipmed,color=(1,0,1))
    plt.axvline(x=precipavg,color=(0,1,0))    
    # set plot mins and maxs
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(xmin,xmax)
    # add x and y labels
    if i in range(0,7,3):
        ax.set_yticks(ylabels)
        ax.set_xticks(xlabels, [])
        if i == 3:
            ax.set_ylabel('Days',fontweight='bold',fontsize=10)
            ax.yaxis.set_label_coords(-0.12,-0.1)
    elif i == numpatterns-2 or i == numpatterns-1:
        ax.set_yticks(ylabels,[],fontsize=9)
        ax.set_xticks(xlabels,fontsize=9)
        if i == numpatterns-2:
            ax.set_xlabel('Precipitation (mm)',fontweight='bold',fontsize=10)
    elif i == numpatterns-3:
        ax.set_yticks(ylabels,fontsize=9)
        ax.set_xticks(xlabels,fontsize=9)
    else:
        ax.set_yticks(ylabels,[])
        ax.set_xticks(xlabels, [])
    ax.tick_params(direction='in',which='both',axis='y')

## PLOT B) DATA
lowlim = 2
highlim = 24

colors = ['lightyellow',"yellow",'greenyellow',"limegreen","lightseagreen",'royalblue','mediumblue','#7400E0','#B800E0','#E0ADB1','lavenderblush'] #,'mediumorchid','#A600E0','pink']
colormap = ListedColormap(colors)

for i, arr in enumerate(som_composites):
    print(np.nanmax(arr))
    #create equidistant cylindrical projection basemap
    map = Basemap(projection='cyl',llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,\
              urcrnrlon=lonmax,resolution='l',area_thresh=area_thresh)
    xi, yi = map(lon,lat)
    locations = [166,176,186,199,209,219,232,242,252,265,275,285]
    startloc2 = locations[i]
    endloc2 = locations[i] + 9
    ax = fig.add_subplot(numrows,numcols,(startloc2,endloc2))
    
    if i == 0:
        ax.set_title('b)',fontsize=12,fontweight="bold",y=0.9,x=-0.16)
    sublabel_loc = mtransforms.ScaledTranslation(4/72, -4/72, fig.dpi_scale_trans)
    ax.text(0.0, 1.0, i+1, transform=ax.transAxes + sublabel_loc,
        fontsize=11, fontweight='bold', verticalalignment='top', 
        bbox=dict(facecolor='1', edgecolor='none', alpha=0.8, pad=1.5),zorder=11)
    #create colormap of MERRA2 data
    colorm = map.pcolor(xi,yi,arr,shading='auto',cmap=colormap,vmin=lowlim,vmax=highlim,zorder=1)
    
    #define border color and thickness
    border_c = '0.4'
    border_w = 0.4
    #create map features
    map.drawcoastlines(color=border_c, linewidth=border_w)
    map.drawstates(color=border_c, linewidth=border_w)
    map.drawcountries(color=border_c, linewidth=border_w)
    gridlinefont = 9
    parallels = np.arange(38.,42.,1.)
    meridians = np.arange(-124.,-119.,2.)
    if i in range(0,7,3):
        map.drawparallels(parallels, labels=[1,0,0,0], fontsize=gridlinefont,color=border_c,linewidth=border_w)
        map.drawmeridians(meridians,color=border_c,linewidth=border_w)
    elif i == numpatterns-2 or i == numpatterns-1:
        map.drawparallels(parallels, color=border_c,linewidth=border_w)
        map.drawmeridians(meridians, labels=[0,0,0,1], fontsize=gridlinefont,color=border_c,linewidth=border_w)
    elif i == numpatterns-3:
        map.drawparallels(parallels, labels=[1,0,0,0], fontsize=gridlinefont,color=border_c,linewidth=border_w)
        map.drawmeridians(meridians, labels=[0,0,0,1], fontsize=gridlinefont,color=border_c,linewidth=border_w)
    else:
        map.drawparallels(parallels, color=border_c,linewidth=border_w)
        map.drawmeridians(meridians,color=border_c,linewidth=border_w)
    #define contour color and thickness
    contour_c = '0.1'
    contour_w = 0.7
    
    if IVT == True:
        #plot IVT vectors
        U_arrs = U_composites[i,:,:]
        V_arrs = V_composites[i,:,:]
        interval = 1
        size = 1500
        skip = (slice(None, None, interval), slice(None, None, interval))
        xim, yim = map(lonm,latm)
        #vectorm = map.quiver(xi2[skip],yi2[skip],U_arrs[skip],V_arrs[skip],color='darkgreen')
        vectorm = map.quiver(xim[skip],yim[skip],U_arrs[skip],V_arrs[skip],pivot='mid', \
                             scale=size, scale_units='inches',headlength=5,headwidth=3, \
                                 color='k',width=0.007,alpha=0.7,zorder=10)

    #add yuba shape
    map.readshapefile(os.path.join(ws_directory,f'{watershed}'), watershed,color='k',linewidth=0.6)
    
    # fig.subplots_adjust(left=0.05,right=0.908,bottom=0.02, top=0.975,hspace=0.05, wspace=0.05) #bottom colorbar
cbar_ax = fig.add_axes([0.911,0.02,0.025,0.46]) #[xloc, yloc, xwidth, yheight]
cbar = fig.colorbar(colorm, cax=cbar_ax,ticks=np.arange(4,highlim+1,2),orientation='vertical')
cbar.ax.tick_params(labelsize=9)
cbar.set_label('Proportion of Wet Season Average',fontsize=10,labelpad=3,fontweight='bold')

    
#CUSTOMIZE SUBPLOT SPACING
fig.subplots_adjust(left=0.07,right=0.985,bottom=0.014, top=0.995,hspace=0.05, wspace=0.5) #bottom colorbar
# fig.subplots_adjust(left=0.065,right=0.985,bottom=0.082, top=0.968,hspace=0.05, wspace=0.05) #bottom colorbar


save_dir='I:\\Emma\\FIROWatersheds\\Figures\\NodeHistograms'
os.chdir(save_dir)
plt.savefig(f'PrecipHistograms_and_SpatialPrecipPatterns_wetseason_{clusters}d_{numpatterns}node',dpi=300)
plt.show()
