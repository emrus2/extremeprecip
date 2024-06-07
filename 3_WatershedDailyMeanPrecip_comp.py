# -*- coding: utf-8 -*-
"""
Created on Sun Jan 15 13:33:53 2023

@author: emrus2
"""
"""
Masks data outside of watershed and calculates spatial average each day. 
Saves to dataset, but needs to be saved multiple times,
 otherwise the data is too big.

Saved Files:
    'I:\\Emma\\FIROWatersheds\\Data\\Gridmet\\UpperYuba\\DailyMeans'
        'UpperYubadailymeans_1979_1996.nc'
        'UpperYubadailymeans_1996_1997.nc'
        'UpperYubadailymeans_1998_2016.nc'
        'UpperYubadailymeans_2017_2022.nc'
"""

#RESTRICT DATA TO ONLY WATERSHED USING MASKING
"""
works on single files
"""

#IMPORT MODULES
import os
import geopandas as gpd
import xarray as xr
from shapely.geometry import mapping
import matplotlib.pyplot as plt
import glob
# import rioxarray as rio

#DEFINE WATERSHED
watershed_name= 'UpperYuba'
#DEFINE YEAR OF PRECIP DATA
year = "1979"


#IMPORT WATERSHED SHAPEFILE
#define CA Watersheds shapefile location
ws_filename = (f'{watershed_name}.shp')
ws_directory = f'I:\\Emma\\FIROWatersheds\\Data\\WatershedShapefiles\\California\\{watershed_name}\\'
ws_path = os.path.join(ws_directory,ws_filename)
upperyuba = gpd.read_file(ws_path, crs="epsg:4326")


#IMPORT NETCDF DATA
nc_filename = (f'pr_{year}.nc')
nc_directory = 'I:\\GRIDMET\\pr\\'
nc_path = os.path.join(nc_directory,nc_filename)
NCdata = xr.open_dataarray(nc_path)
NCdata.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
NCdata.rio.write_crs("epsg:4326", inplace=True)


#REDUCE NC DATA TO SHAPEFILE OF UPPERYUBA and SAVE AS NEW NC FILE
watershed = NCdata.rio.clip(upperyuba.geometry.apply(mapping), upperyuba.crs, drop=False)
#newnc_dir = 'I:\\Emma\\FIROWatersheds\\Data\\GridmetWatershed\\'
#watershed.to_netcdf(os.path.join(newnc_dir,f'pr_{year}_{watershed_name}'))

#REDUCE DATA TO DAILY MEANS >2mm
daily_means = watershed.mean(dim=['lat','lon'],skipna=True)
daily_means_thresh = daily_means.where(daily_means > 2)
# daily_means_thresh = np.where(daily_means, daily_means>2)

#%%

""" works to combine multiple files limited to watershed """
#IMPORT MODULES
import os
import geopandas as gpd
import xarray as xr
from shapely.geometry import mapping

#DEFINE WATERSHED
watershed_name= 'UpperYuba'
#DEFINE YEAR OF PRECIP DATA
#IMPORT WATERSHED SHAPEFILE
#define CA Watersheds shapefile location
ws_filename = (f'{watershed_name}.shp')
ws_directory = f'I:\\Emma\\FIROWatersheds\\Data\\WatershedShapefiles\\California\\{watershed_name}\\'
ws_path = os.path.join(ws_directory,ws_filename)
upperyuba = gpd.read_file(ws_path, crs="epsg:4326")

#IMPORT NETCDF DATA
year = "1979"
nc_filename = (f'pr_{year}.nc')
nc_directory = 'I:\\GRIDMET\\pr\\'
nc_path = os.path.join(nc_directory,nc_filename)
NCdata = xr.open_dataarray(nc_path)
NCdata.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
NCdata.rio.write_crs("epsg:4326", inplace=True)
watershed = NCdata.rio.clip(upperyuba.geometry.apply(mapping), upperyuba.crs, drop=False)
daily_means = watershed.mean(dim=['lat','lon'],skipna=True)
daily_means_thresh1 = daily_means.where(daily_means > 2)

newnc_dir = f'I:\\Emma\\FIROWatersheds\\Data\\Gridmet\\{watershed_name}\\DailyMeans'

for year in range(1980,1998):
    print(daily_means_thresh1)
    year = str(year)
    print(year)
    try:
        nc_filename = (f'pr_{year}.nc')
        nc_directory = 'I:\\GRIDMET\\pr\\'
        nc_path = os.path.join(nc_directory,nc_filename)
        NCdata = xr.open_dataarray(nc_path)
        NCdata.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
        NCdata.rio.write_crs("epsg:4326", inplace=True)
        watershed = NCdata.rio.clip(upperyuba.geometry.apply(mapping), upperyuba.crs, drop=False)
        daily_means = watershed.mean(dim=['lat','lon'],skipna=True)
        daily_means_thresh = daily_means.where(daily_means > 2)
        daily_means_thresh1 = xr.concat([daily_means_thresh1,daily_means_thresh],dim='day')
    except:
        daily_means_thresh1.to_netcdf(os.path.join(newnc_dir,f'{watershed_name}dailymeans_1979_{year-1}.nc'))    


daily_means_thresh1.to_netcdf(os.path.join(newnc_dir,f'{watershed_name}dailymeans_1979_{year}.nc'))

#%%
""" works to combine multiple files limited to watershed """
#IMPORT MODULES
import os
import geopandas as gpd
import xarray as xr
from shapely.geometry import mapping

#DEFINE WATERSHED
watershed_name= 'UpperYuba'
#DEFINE YEAR OF PRECIP DATA
#IMPORT WATERSHED SHAPEFILE
#define CA Watersheds shapefile location
ws_filename = (f'{watershed_name}.shp')
ws_directory = f'I:\\Emma\\FIROWatersheds\\Data\\WatershedShapefiles\\California\\{watershed_name}\\'
ws_path = os.path.join(ws_directory,ws_filename)
upperyuba = gpd.read_file(ws_path, crs="epsg:4326")

#IMPORT NETCDF DATA
year = "1996"
nc_filename = (f'pr_{year}.nc')
nc_directory = 'I:\\GRIDMET\\pr\\'
nc_path = os.path.join(nc_directory,nc_filename)
NCdata = xr.open_dataarray(nc_path)
NCdata.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
NCdata.rio.write_crs("epsg:4326", inplace=True)
watershed = NCdata.rio.clip(upperyuba.geometry.apply(mapping), upperyuba.crs, drop=False)
daily_means = watershed.mean(dim=['lat','lon'],skipna=True)
daily_means_thresh1 = daily_means.where(daily_means > 2)

newnc_dir = f'I:\\Emma\\FIROWatersheds\\Data\\Gridmet\\{watershed_name}\\DailyMeans'

for year in range(1997,1998):
    print(daily_means_thresh1)
    year = str(year)
    print(year)
    try:
        nc_filename = (f'pr_{year}.nc')
        nc_directory = 'I:\\GRIDMET\\pr\\'
        nc_path = os.path.join(nc_directory,nc_filename)
        NCdata = xr.open_dataarray(nc_path)
        NCdata.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
        NCdata.rio.write_crs("epsg:4326", inplace=True)
        watershed = NCdata.rio.clip(upperyuba.geometry.apply(mapping), upperyuba.crs, drop=False)
        daily_means = watershed.mean(dim=['lat','lon'],skipna=True)
        daily_means_thresh = daily_means.where(daily_means > 2)
        daily_means_thresh1 = xr.concat([daily_means_thresh1,daily_means_thresh],dim='day')
    except:
        daily_means_thresh1.to_netcdf(os.path.join(newnc_dir,f'{watershed_name}dailymeans_1996_{year}.nc'))    


daily_means_thresh1.to_netcdf(os.path.join(newnc_dir,f'{watershed_name}dailymeans_1996_{year}.nc'))

#%%
""" works to combine multiple files limited to watershed """
#IMPORT MODULES
import os
import geopandas as gpd
import xarray as xr
from shapely.geometry import mapping

#DEFINE WATERSHED
watershed_name= 'UpperYuba'
#DEFINE YEAR OF PRECIP DATA
#IMPORT WATERSHED SHAPEFILE
#define CA Watersheds shapefile location
ws_filename = (f'{watershed_name}.shp')
ws_directory = f'I:\\Emma\\FIROWatersheds\\Data\\WatershedShapefiles\\California\\{watershed_name}\\'
ws_path = os.path.join(ws_directory,ws_filename)
upperyuba = gpd.read_file(ws_path, crs="epsg:4326")

#IMPORT NETCDF DATA
year = "1998"
nc_filename = (f'pr_{year}.nc')
nc_directory = 'I:\\GRIDMET\\pr\\'
nc_path = os.path.join(nc_directory,nc_filename)
NCdata = xr.open_dataarray(nc_path)
NCdata.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
NCdata.rio.write_crs("epsg:4326", inplace=True)
watershed = NCdata.rio.clip(upperyuba.geometry.apply(mapping), upperyuba.crs, drop=False)
daily_means = watershed.mean(dim=['lat','lon'],skipna=True)
daily_means_thresh1 = daily_means.where(daily_means > 2)

newnc_dir = f'I:\\Emma\\FIROWatersheds\\Data\\Gridmet\\{watershed_name}\\DailyMeans'

for year in range(1999,2017):
    print(daily_means_thresh1)
    year = str(year)
    print(year)
    try:
        nc_filename = (f'pr_{year}.nc')
        nc_directory = 'I:\\GRIDMET\\pr\\'
        nc_path = os.path.join(nc_directory,nc_filename)
        NCdata = xr.open_dataarray(nc_path)
        NCdata.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
        NCdata.rio.write_crs("epsg:4326", inplace=True)
        watershed = NCdata.rio.clip(upperyuba.geometry.apply(mapping), upperyuba.crs, drop=False)
        daily_means = watershed.mean(dim=['lat','lon'],skipna=True)
        daily_means_thresh = daily_means.where(daily_means > 2)
        daily_means_thresh1 = xr.concat([daily_means_thresh1,daily_means_thresh],dim='day')
    except:
        daily_means_thresh1.to_netcdf(os.path.join(newnc_dir,f'{watershed_name}dailymeans_1998_{year}.nc'))    


daily_means_thresh1.to_netcdf(os.path.join(newnc_dir,f'{watershed_name}dailymeans_1998_{year}.nc'))

#%%
""" works to combine multiple files limited to watershed """
#IMPORT MODULES
import os
import geopandas as gpd
import xarray as xr
from shapely.geometry import mapping

#DEFINE WATERSHED
watershed_name= 'UpperYuba'
#DEFINE YEAR OF PRECIP DATA
#IMPORT WATERSHED SHAPEFILE
#define CA Watersheds shapefile location
ws_filename = (f'{watershed_name}.shp')
ws_directory = f'I:\\Emma\\FIROWatersheds\\Data\\WatershedShapefiles\\California\\{watershed_name}\\'
ws_path = os.path.join(ws_directory,ws_filename)
upperyuba = gpd.read_file(ws_path, crs="epsg:4326")

#IMPORT NETCDF DATA
year = "2017"
nc_filename = (f'pr_{year}.nc')
nc_directory = 'I:\\GRIDMET\\pr\\'
nc_path = os.path.join(nc_directory,nc_filename)
NCdata = xr.open_dataarray(nc_path)
NCdata.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
NCdata.rio.write_crs("epsg:4326", inplace=True)
watershed = NCdata.rio.clip(upperyuba.geometry.apply(mapping), upperyuba.crs, drop=False)
daily_means = watershed.mean(dim=['lat','lon'],skipna=True)
daily_means_thresh1 = daily_means.where(daily_means > 2)

newnc_dir = f'I:\\Emma\\FIROWatersheds\\Data\\Gridmet\\{watershed_name}\\DailyMeans'

for year in range(2018,2023):
    print(daily_means_thresh1)
    year = str(year)
    print(year)
    try:
        nc_filename = (f'pr_{year}.nc')
        nc_directory = 'I:\\GRIDMET\\pr\\'
        nc_path = os.path.join(nc_directory,nc_filename)
        NCdata = xr.open_dataarray(nc_path)
        NCdata.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
        NCdata.rio.write_crs("epsg:4326", inplace=True)
        watershed = NCdata.rio.clip(upperyuba.geometry.apply(mapping), upperyuba.crs, drop=False)
        daily_means = watershed.mean(dim=['lat','lon'],skipna=True)
        daily_means_thresh = daily_means.where(daily_means > 2)
        daily_means_thresh1 = xr.concat([daily_means_thresh1,daily_means_thresh],dim='day')
    except:
        daily_means_thresh1.to_netcdf(os.path.join(newnc_dir,f'{watershed_name}dailymeans_2017_{year}.nc'))    


daily_means_thresh1.to_netcdf(os.path.join(newnc_dir,f'{watershed_name}dailymeans_2017_{year}.nc'))
