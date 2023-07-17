# -*- coding: utf-8 -*-
"""
Created on Sun Jan 15 13:33:53 2023

Saves a netcdf file of gridmet data reduced to region around dataset
Can be used for several different meteorological variables of gridMET

@author: emrus2
"""
#%%
#REDUCE NCDATA TO WATERSHED REGION
""" 
Creates a smaller netcdf file of the data around the wateshed
    need to adjust the lat and lon limits
"""
    
#IMPORT MODULES
import os
import netCDF4 as nc
import numpy as np

#DEFINE WATERSHED NAME
watershed = 'UpperYuba'

#IMPORT NETCDF DATA
year = '1997'
#metvars = ['pr','rmax','rmin','vpd','tmmn','tmmx','vs']
metvars = ['vs']
varnames = ['precipitation_amount','relative_humidity','relative_humidity', \
            'air_temperature','air_temperature','wind_speed']
unitnames = ['mm','%','%','K','K','m/s']

i = 5
for metvar in metvars:
    nc_filename = (f'{metvar}_{year}.nc')
    if metvar == 'tmmn' or metvar == 'tmmx':
        nc_directory = 'I:\\GRIDMET\\temp\\'
    elif metvar == 'vs':
        nc_directory = 'I:\\GRIDMET\\vs_(10m_wind_speed)\\'
    else:
        nc_directory = f'I:\\GRIDMET\\{metvar}\\'
    filepath = os.path.join(nc_directory,nc_filename)
    gridfile = nc.Dataset(filepath,mode='r')
    print(gridfile)
    
    gridlat = gridfile.variables['lat'][:]
    gridlon = gridfile.variables['lon'][:]
    gridday = gridfile.variables['day'][:]
    gridcrs = gridfile.variables['crs'][:]
    precip = gridfile.variables[varnames[i]][:,:,:] #time x height x lat x lon
    gridfile.close()
    
    #REDUCE VARIABLES TO DESIRED AREA
    #convert height to a 3D array
    precipsm = np.squeeze(precip)
    #define lat lon restrictions
    latmin = 38 #latmin = 29.75
    latmax = 42 #latmax = 41.75
    lonmin = -124 #lonmin = -130.75
    lonmax = -119 #lonmax = -110.25
    
    #reduce lat
    latlims = np.logical_and(gridlat > latmin, gridlat < latmax)
    latind = np.where(latlims)[0]
    latreduced = gridlat[latind]
    #reduce lon
    lonlims = np.logical_and(gridlon > lonmin, gridlon < lonmax)
    lonind = np.where(lonlims)[0]
    lonreduced = gridlon[lonind]
    #reduce pressure
    precipreduced = precipsm[:,latind,:]
    precipreduced = precipreduced[:,:,lonind]
    
    
    #SAVE NEW DATA TO NETCDF
    #define new file name and directory
    newnc_file = f'GRIDMET_{metvar}_{year}.nc'
    newnc_dir = 'I:\\Emma\\FIROWatersheds\\Data\\GridmetWatershed\\'
    #create new nc file with designated filename and location
    newnc = nc.Dataset(os.path.join(newnc_dir,newnc_file),'w',format='NETCDF4')
    
    #create dimensions of nc file
    day = newnc.createDimension('day', None)
    lat = newnc.createDimension('lat', len(latreduced))
    lon = newnc.createDimension('lon',len(lonreduced))
    crs = newnc.createDimension('crs', len(gridcrs))
    
    #create attributes of nc file
    newnc.title='Sea Level Pressure over Upper Yuba Watershed'
    newnc.subtitle='Restricted to region over Upper Yuba'
    
    #create variables of nc file
    latitude = newnc.createVariable('lat', 'f8', ('lat',))  
    latitude.units = 'degrees_north'
    longitude = newnc.createVariable('lon', 'f8', ('lon',))
    latitude.units = 'degrees_east'
    day = newnc.createVariable('day', 'f8', ('day',))
    day.units = 'days since 1900-01-01 00:00:00'
    crs = newnc.createVariable('crs', 'i2', ('crs',))
    precipitation_amount = newnc.createVariable(varnames[i], 'i2', ('day','lat', 'lon'))
    precipitation_amount.units = unitnames[i]
    
    #input data to nc file
    latitude[:] = latreduced
    longitude[:] = lonreduced
    day[:] = gridday
    crs[:] = gridcrs
    precipitation_amount[:,:,:] = precipreduced
    
    print(newnc)
    newnc.close()
    
    i += 1

#%%
#REDUCE NCDATA TO WATERSHED REGION
""" 
Creates a smaller netcdf file of the data around the wateshed
    need to adjust the lat and lon limits
"""
    
#IMPORT MODULES
import os
import netCDF4 as nc
import numpy as np

#DEFINE WATERSHED NAME
watershed = 'UpperYuba'

#IMPORT NETCDF DATA
year = '1997'
metvar = 'temp'

nc_filename = (f'{metvar}_{year}.nc')
nc_directory = f'I:\\GRIDMET\\{metvar}\\'
filepath = os.path.join(nc_directory,nc_filename)
gridfile = nc.Dataset(filepath,mode='r')
print(gridfile)
gridlat = gridfile.variables['lat'][:]
gridlon = gridfile.variables['lon'][:]
gridday = gridfile.variables['day'][:]
gridcrs = gridfile.variables['crs'][:]
precip = gridfile.variables['metvar'][:,:,:] #time x height x lat x lon
gridfile.close()



#REDUCE VARIABLES TO DESIRED AREA
#convert height to a 3D array
precipsm = np.squeeze(precip)
#define lat lon restrictions
latmin = 38 #latmin = 29.75
latmax = 42 #latmax = 41.75
lonmin = -124 #lonmin = -130.75
lonmax = -119 #lonmax = -110.25

#reduce lat
latlims = np.logical_and(gridlat > latmin, gridlat < latmax)
latind = np.where(latlims)[0]
latreduced = gridlat[latind]
#reduce lon
lonlims = np.logical_and(gridlon > lonmin, gridlon < lonmax)
lonind = np.where(lonlims)[0]
lonreduced = gridlon[lonind]
#reduce pressure
precipreduced = precipsm[:,latind,:]
precipreduced = precipreduced[:,:,lonind]


#SAVE NEW DATA TO NETCDF
#define new file name and directory
newnc_file = f'GRIDMET_{metvar}_{year}.nc'
newnc_dir = 'I:\\Emma\\FIROWatersheds\\Data\\GridmetWatershed\\'
#create new nc file with designated filename and location
newnc = nc.Dataset(os.path.join(newnc_dir,newnc_file),'w',format='NETCDF4')

#create dimensions of nc file
day = newnc.createDimension('day', None)
lat = newnc.createDimension('lat', len(latreduced))
lon = newnc.createDimension('lon',len(lonreduced))
crs = newnc.createDimension('crs', len(gridcrs))

#create attributes of nc file
newnc.title='Sea Level Pressure over Upper Yuba Watershed'
newnc.subtitle='Restricted to region over Upper Yuba'

#create variables of nc file
latitude = newnc.createVariable('lat', 'f8', ('lat',))  
latitude.units = 'degrees_north'
longitude = newnc.createVariable('lon', 'f8', ('lon',))
latitude.units = 'degrees_east'
day = newnc.createVariable('day', 'f8', ('day',))
day.units = 'days since 1900-01-01 00:00:00'
crs = newnc.createVariable('crs', 'i2', ('crs',))
precipitation_amount = newnc.createVariable('precipitation_amount', 'i2', ('day','lat', 'lon'))
precipitation_amount.units = 'mm'

#input data to nc file
latitude[:] = latreduced
longitude[:] = lonreduced
day[:] = gridday
crs[:] = gridcrs
precipitation_amount[:,:,:] = precipreduced

print(newnc)
newnc.close()

#%%
#RESTRICT DATA TO ONLY WATERSHED

"""
Masks all other nc data other than what is included in the shapefile

This for some reason only works on work laptop, issue importing rioxarray otherwise
"""

#IMPORT MODULES
import os
import geopandas as gpd
import xarray
from shapely.geometry import mapping
import rioxarray

#DEFINE WATERSHED
watershed= 'UpperYuba'

year = '1997'
#metvars = ['pr','rmax','rmin','vpd']
metvars = ['vs']

for metvar in metvars:
    #IMPORT UPPER YUBA SHAPEFILE
    #define CA Watersheds shapefile location
    ws_filename = (f'{watershed}.shp')
    ws_directory = f'I:\\Emma\\FIROWatersheds\\Data\\WatershedShapefiles\\California\\{watershed}\\'
    ws_path = os.path.join(ws_directory,ws_filename)
    upperyuba = gpd.read_file(ws_path, crs="epsg:4326")
    
    
    #IMPORT NETCDF DATA
    newnc_file = f'GRIDMET_{metvar}_{year}.nc'
    newnc_dir = 'I:\\Emma\\FIROWatersheds\\Data\\GridmetWatershed\\'
    nc_path = os.path.join(newnc_dir,newnc_file)
    NCdata = xarray.open_dataarray(nc_path)
    NCdata.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
    NCdata.rio.write_crs("epsg:4326", inplace=True)
    
    #REDUCE NC DATA TO SHAPEFILE OF UPPERYUBA and SAVE AS NEW NC FILE
    #NCclipped = NCdata.rio.clip(upperyuba.geometry.apply(mapping), upperyuba.crs, drop=False)
    NCclipped = NCdata.rio.clip(upperyuba.geometry.apply(mapping), upperyuba.crs, drop=True)
    NCclipped.to_netcdf(os.path.join(newnc_dir,f'GRIDMET_{metvar}_{watershed}_{year}.nc'))
