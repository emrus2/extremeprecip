# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 15:35:52 2024

@author: emrus2
"""
#%% IMPORT PACKAGES
import os
import netCDF4 as nc
import numpy as np
from matplotlib.colors import ListedColormap
import matplotlib.cm as cm
import matplotlib.transforms as mtransforms
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
from mpl_toolkits.basemap import Basemap #installed using 
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import shapefile 

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


#%% 
# import extreme dates of upper yuba
extremedates = np.load('I:/Emma/FIROWatersheds/Data/90Percentile_ExtremeDays.npy',allow_pickle=True)

#%% OPEN 850TAnom Data
percentile = 90
folderpath = '850T'
filename = f'MERRA2_850T_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST_updated.nc'
filepath = os.path.join(folderpath,filename)
merravar = 'T'
    
os.chdir('I:\\Emma\\FIROWatersheds\\Data\\DailyMERRA2')
gridfile = nc.Dataset(filepath,mode='r')
print(gridfile.variables)
gridlat = gridfile.variables['lat'][:]
gridlon = gridfile.variables['lon'][:]
time = gridfile.variables['date'][:]
gridfile.close()

# import anomaly data
merra = np.load(os.path.join(folderpath,f'MERRA2_850TAnom_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST_updated.npy'))

# convert time from hours to datetimes
# time = [datetime.strptime('19800105','%Y%m%d') + timedelta(days=int(j)) for j in time]
# timestr = [datetime.strftime(day,'%Y%m%d') for day in time]
#%% REDUCE CATALOG TO EXTREME DATES
# find extreme dates in AR catalog dataset
#reduce dates
# datelims = np.isin(timestr,extremedates) # find dates that are present in extreme dates list
# dateind = np.where(datelims)[0] # find those indices
# datereduced = timestr[dateind] # reduce dates to indices of those presnt in extreme dates list
# dates_red = sorted(list(set(datereduced))) # sort dates back to chronological order
    
#%% REDUCE AR CATALOG DATA TO WATERSHED REGION
# import watershed data
watershed_name= 'UpperYuba'
ws_directory = f'I:\\Emma\\FIROWatersheds\\Data\\WatershedShapefiles\\California\\{watershed_name}\\'
ws_filename = (f'{watershed_name}.shp')
ws_path = os.path.join(ws_directory,ws_filename)
shp = shapefile.Reader(ws_path, crs="epsg:4326")
feature = shp.shapes()

# Extract array of lat/lon coordinates
coords = np.squeeze(np.array([s.points for s in shp.shapes()]))
shp_lon = np.array(coords)[:,0]
shp_lat = np.array(coords)[:,1]

# set gridcells of watershed region
latminws, latmaxws = (38.5,40.5)
lonminws, lonmaxws = (-122.5,-119.9)

# reduce lat
latlimsws = np.logical_and(gridlat > latminws, gridlat < latmaxws)
latindws = np.where(latlimsws)[0]
gridlatws = gridlat[latindws]
# reduce lon
lonlimsws = np.logical_and(gridlon > lonminws, gridlon < lonmaxws)
lonindws = np.where(lonlimsws)[0]
gridlonws = gridlon[lonindws]
#%%
# reduce AR data
tempwatershed = merra[:,latindws,:]
tempwatershed = tempwatershed[:,:,lonindws]

#%% IDENTIFY AR PRESENCE IN WATERSHED REGION
# add up AR frequency across watershed gridcells
watershedmean = np.nanmean(tempwatershed,axis=(1,2))
watershedstr = [str(val) for val in watershedmean]
extremedateslist = [str(i) for i in extremedates]
# combine ar presence data with dates data
dates_temp = np.array(list(zip(extremedateslist,watershedstr)))

# save as numpy array
savefolder = 'I:/Emma/Classes/SnowHydrology/'
np.savetxt(os.path.join(savefolder,'850TempAvgAnomalies_90_Percentile_Days.csv'), \
           dates_temp,delimiter=",",fmt='%s')

#%% PLOT AR FREQUENCY TO CONFIRM ALGORITHM
for i in range(len(tempwatershed)):
# for i in dated:
    latmin, latmax = (25,55)
    lonmin, lonmax = (-150,-110)
    
    latmin, latmax = (38.5,40.5)
    lonmin, lonmax = (-122.5,-119.9)

    
    #reduce lat
    latlims = np.logical_and(gridlat > latmin, gridlat < latmax)
    latind = np.where(latlims)[0]
    gridlatreduced = gridlat[latind]
    #reduce lon
    lonlims = np.logical_and(gridlon > lonmin, gridlon < lonmax)
    lonind = np.where(lonlims)[0]
    gridlonreduced = gridlon[lonind]
    
    # datareduced = merra[i,:,:]
    # datareduced = datareduced[latind,:]
    # datareduced = datareduced[:,lonind]
    
    datareduced = tempwatershed[i,:,:]

    #plot
    fig = plt.figure(figsize=(7.2,5.1))
    plt.title(f'{watershedsum[i]}')
    #MAP DESIRED VARIABLE
    #convert lat and lon into a 2D array
    lon, lat = np.meshgrid(gridlonreduced,gridlatreduced)  
    #define area threshold for basemap
    area_thresh = 1E4
    #create equidistant cylindrical projection basemap
    map = Basemap(projection='cyl',llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,\
              urcrnrlon=lonmax,resolution='c',area_thresh=area_thresh)
    xi, yi = map(lon,lat)
    border_c = '0.4'
    border_w = 0.4
    #create map features
    try:
        map.drawcoastlines(color=border_c, linewidth=border_w,zorder=2)
    except:
        pass
    map.drawstates(color=border_c, linewidth=border_w,zorder=2)
    map.drawcountries(color=border_c, linewidth=border_w,zorder=2)
    
    parallels = np.arange(20.,71.,2.)
    meridians = np.arange(-160.,-109.,2.)
    map.drawparallels(parallels, labels=[1,0,0,0])
    map.drawmeridians(meridians, labels=[0,0,0,1])
    
    map.readshapefile(os.path.join(ws_directory,f'{watershed}'), watershed,zorder=6,color='k')
    plt.scatter(360-120.9,39.5,color='b',marker='*',linewidths=0.7,zorder=4)
    
    
    #create colormap of MERRA2 data
    colorm = map.pcolor(xi,yi,datareduced,shading='auto',cmap=anom_cm('850TAnom')[2], \
             vmin=anom_cm('850TAnom')[0],vmax=anom_cm('850TAnom')[1])           
    #define border color and thickness
    plt.colorbar()
    plt.show()