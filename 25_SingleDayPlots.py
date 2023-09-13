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
#%%
# define lat, lon region of data for plotting
latmin, latmax = (15.5,65.5)
lonmin, lonmax = (-170.25,-105.75)
percentile = '90'

#%% IMPORT EXTREME PRECIPITATION DATE DATA
os.chdir('I:\\Emma\\FIROWatersheds\\Data')
dates = np.load('90Percentile_ExtremeDays.npy')

#%% IMPORT MERRA2 DATA
# define metvar
#metvars = ['IVT','300W','Z500Anom','SLP','SLPAnom','Z850','850T','850TAnom']
metvar = '850QVECT'
#define composite location
#metvar = input('Enter MERRA-2 Variable: Z500, SLP, 850T, 300W, or IVT \n')
if metvar == 'Z500Anom':
    folderpath = 'DailyMERRA2\\Z500'
    filename = f'MERRA2_Z500_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST.nc'
elif metvar == 'SLPAnom':
    folderpath = 'DailyMERRA2\\SLP'
    filename = f'MERRA2_SLP_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST.nc'
elif metvar == '850TAnom':
    folderpath = 'DailyMERRA2\\850T'
    filename = f'MERRA2_850T_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST.nc'
else:
    folderpath = f'DailyMERRA2\\{metvar}'
    filename = f'MERRA2_{metvar}_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST_updated.nc'
filepath = os.path.join(folderpath,filename)


#COLLECT VARIABLE DATA FROM MERRA2 FILE
merravar = {'Z500':'H','SLP':'SLP','850T':'T','Z850':'H','850TADV':'__xarray_dataarray_variable__'}
#open the netcdf file in read mode
gridfile = nc.Dataset(filepath,mode='r')
print(gridfile)
gridlat = gridfile.variables['lat'][:]
gridlon = gridfile.variables['lon'][:]
if metvar == 'IVT':
    Uvapor = gridfile.variables['UFLXQV'][:]
    Vvapor = gridfile.variables['VFLXQV'][:]
    merra = np.sqrt(Uvapor**2 + Vvapor**2)
elif metvar == '300W' or metvar == '850QVECT':
    U = gridfile.variables['U'][:]
    V = gridfile.variables['V'][:]
    merra = np.sqrt(U**2 + V**2)
elif 'Anom' in metvar:
    merra = np.load(os.path.join(folderpath,f'MERRA2_{metvar}_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST.npy'))
else:
    merra = gridfile.variables[merravar[metvar]][:]
gridfile.close()

merra = np.squeeze(merra)

if metvar == 'SLP':
    merra = merra/100
    
if metvar == '850TADV':
    filename2 = f'DailyMERRA2\\850QVECT\\MERRA2_850UQVECT_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST.nc' 
    gridfile2 = nc.Dataset(filename2,mode='r')
    Uvect = gridfile.variables[merravar[metvar]][:]
    gridfile2.close()
    filename3 = f'DailyMERRA2\\850QVECT\\MERRA2_850VQVECT_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST.nc' 
    gridfile3 = nc.Dataset(filename3,mode='r')
    Vvect = gridfile.variables[merravar[metvar]][:]
    gridfile3.close()
    merra = merra * 3600

# if metvar == '850QVECT':
#     Ugrad = np.ufunc.reduce(np.add,np.gradient(U),axis=2)
#     Vgrad = np.ufunc.reduce(np.add,np.gradient(V),axis=1)
#     merra = Ugrad + Vgrad
#     # merra = np.ufunc.reduce(np.add,np.gradient(U)) + np.ufunc.reduce(np.add,np.gradient(V))
#     merra *= 10**15

def divergence(F):
    """ Computes divergence of vector field 
    f: array -> vector field components [Fx,Fy,Fz,...]
    sp: array -> spacing between points in respecitve directions [spx, spy,spz,...]
    """
    num_dims = len(F)
    return np.ufunc.reduce(np.add, [np.gradient(F[i], axis=i) for i in range(num_dims)])
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
#reduce winds
Ureduced = U[:,latind,:]
Ureduced = Ureduced[:,:,lonind]
Vreduced = V[:,latind,:]
Vreduced = Vreduced[:,:,lonind]
#print(np.amin(merrareduced),np.amax(merrareduced))

#%%
# define day of interest
for day in np.arange(260):
    merraday = merrareduced[day,:,:] # K/s
    Uday = Ureduced[day,:,:]
    Vday = Vreduced[day,:,:]
    F = np.array([Uday,Vday])
    merraday = divergence(F)

    #print(np.nanmin(merraday),np.nanmax(merraday))

    # DETERMINE MAX AND MIN VALIUES
    zmax = 0
    zmin = 1E8
    
    #determine zmax and zmin for all days
    highlim = np.nanmax(merraday)
    lowlim = np.nanmin(merraday)
    #print(lowlim,highlim)
    if highlim > zmax:
        zmax = highlim
    if lowlim < zmin:
        zmin = lowlim
    
    #print(f'Lowest Value:{zmin} \nHighest Value:{zmax}')
    zmin, zmax = (-2,2)
    newmap = center_colormap(zmin, zmax, center=0)
    # PLOT
    
    #create subplot for mapping multiple timesteps
    fig = plt.figure(figsize=(7.2,5))
    plt.title(dates[day])
    #MAP DESIRED VARIABLE
    #convert lat and lon into a 2D array
    lon, lat = np.meshgrid(gridlonreduced,gridlatreduced) 
    #define area threshold for basemap
    area_thresh = 1E4
    #create equidistant cylindrical projection basemap
    map = Basemap(projection='cyl',llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,\
              urcrnrlon=lonmax,resolution='l',area_thresh=area_thresh)
    xi, yi = map(lon,lat)
    #create colormap of MERRA2 data
    #colorm = map.pcolor(xi,yi,merraday,cmap=newmap,vmin=zmin,vmax=zmax)
    colorm = map.pcolor(xi,yi,merraday,cmap=newmap)    
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
    contour_w = 0.7
    #create contour map
    #contourm = map.contour(xi,yi,merraday,colors=contour_c,linewidths=contour_w)
    #plt.clabel(contourm,levels=contourm.levels[::2],fontsize=6,inline_spacing=1,\
    #               colors='k',zorder=2,manual=False)
    plt.scatter(-120.9,39.5,color='w',marker='*',linewidths=0.7,zorder=4)
        
    interval = 5
    size = 2
    skip = (slice(None, None, interval), slice(None, None, interval))
    vectorm = map.quiver(xi[skip],yi[skip],Uday[skip],Vday[skip],color='darkgreen')
    #vectorm = map.quiver(xi[skip],yi[skip],U_arrs[skip],V_arrs[skip],pivot='mid',scale=size, scale_units='inches',headlength=3.4,headwidth=2,color='g',width=0.005)

    #CUSTOMIZE SUBPLOT SPACING
    cbar_ax = fig.add_axes([0.905,0.05,0.025,0.9]) #bottom colorbar
    cbar = fig.colorbar(colorm, cax=cbar_ax,orientation='vertical',extend='both')
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label(label='Temperature Advection (K/s)',fontsize=8.5,labelpad=0.5,fontweight='bold')
    
        
    #SHOW MAP
    save_dir='I:\\Emma\\FIROWatersheds\\Figures\\'
    os.chdir(save_dir)
    #plt.savefig(f'{metvar}_{day}',dpi=300)
    plt.show()