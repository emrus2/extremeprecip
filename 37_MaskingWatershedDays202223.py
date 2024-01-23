# -*- coding: utf-8 -*-
"""
Created on Sun Jan 15 13:33:53 2023

@author: emrus2
"""
"""
Masks data outside of watershed for 2022-2023 files
The files saved:
    'UpperYubadailymeans_2022-2023.nc'
"""

#%%

""" works to combine multiple files limited to watershed """
#IMPORT MODULES
import os
import xarray as xr
from shapely.geometry import mapping
import shapefile

#DEFINE WATERSHED
watershed_name= 'UpperYuba'
#DEFINE YEAR OF PRECIP DATA
#IMPORT WATERSHED SHAPEFILE
#define CA Watersheds shapefile location
ws_filename = (f'{watershed_name}.shp')
ws_directory = f'I:\\Emma\\FIROWatersheds\\Data\\WatershedShapefiles\\California\\{watershed_name}\\'
ws_path = os.path.join(ws_directory,ws_filename)
upperyuba = shapefile.Reader(ws_path, crs="epsg:4326")

#IMPORT NETCDF DATA
year = "2022"
nc_filename = (f'pr_{year}.nc')
nc_directory = 'I:\\GRIDMET\\pr\\'
nc_path = os.path.join(nc_directory,nc_filename)
NCdata = xr.open_dataarray(nc_path)
NCdata.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
NCdata.rio.write_crs("epsg:4326", inplace=True)
watershed = NCdata.rio.clip(upperyuba.geometry.apply(mapping), upperyuba.crs, drop=False)
daily_means = watershed.mean(dim=['lat','lon'],skipna=True)
daily_means_thresh1 = daily_means.where(daily_means > 2)

newnc_dir = 'I:\\Emma\\FIROWatersheds\\Data\\GridmetWatershed\\'

for year in range(2023,2024):
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
        daily_means_thresh1.to_netcdf(os.path.join(newnc_dir,f'{watershed_name}dailymeans_2022_{year-1}.nc'))    


daily_means_thresh1.to_netcdf(os.path.join(newnc_dir,f'{watershed_name}dailymeans_2022_{year}.nc'))