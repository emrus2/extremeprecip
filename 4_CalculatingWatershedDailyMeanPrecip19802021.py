# -*- coding: utf-8 -*-
"""
Created on Sun Jan 15 13:33:53 2023

@author: emrus2
"""
"""
Now have a saved file of all the daily precipitation means.
Time: 1/1/1979 - 12/17/2022 (total of 16057 timesteps)
Precipitation values are the mean over the Upper Yuba watershed region.
Daily mean values less than 2 mm were removed from the dataset.
The 95th percentile of the daily mean distribution was calculated.
Extreme precipitation days greater than this 95th percentile were also saved

The file was then reduced to begin at 01-01-1980 
and end at 12-31-2021 to match MERRA2

The files saved:
    UpperYubadailymeanTOTAL.nc
    UpperYubadailymeanEXTREMES.nc
    UpperYubadailymeanEXTREMES_1980_2021.nc
"""

#%%
#IMPORT NETCDF DATA

import os
import xarray as xr

watershed = "UpperYuba"

newnc_dir = f'I:\\Emma\\FIROWatersheds\\Data\\Gridmet\\{watershed}\\DailyMeans'

file1 = f'{watershed}dailymeans_1979_1996.nc'
file15 = f'{watershed}dailymeans_1996_1997.nc'
file2 = f'{watershed}dailymeans_1998_2016.nc'
file3 = f'{watershed}dailymeans_2017_2022.nc'

nc_path1 = os.path.join(newnc_dir,file1)
nc_path15 = os.path.join(newnc_dir,file15)
nc_path2 = os.path.join(newnc_dir,file2)
nc_path3 = os.path.join(newnc_dir,file3)

NCdata1 = xr.open_dataarray(nc_path1)
NCdata15 = xr.open_dataarray(nc_path15)
NCdata2 = xr.open_dataarray(nc_path2)
NCdata3 = xr.open_dataarray(nc_path3)

#total netcdf file
NCdatatotal = xr.concat([NCdata1,NCdata15,NCdata2,NCdata3],dim='day')
#combined should be 16057 in length

NCdatatotal.to_netcdf(os.path.join(newnc_dir,f'{watershed}dailymeansTOTAL.nc'))

#percent_95 = NCdatatotal.quantile(0.95)
#print(percent_95)

#extremes = NCdatatotal.where(NCdatatotal>percent_95)
#NCdatatotal.plot.line(marker='o',lw=0,markersize=2)
#extremes.plot.line(marker='o',lw=0,markersize=2)
#xr.plot.hist(NCdatatotal)
#xr.plot.hist(extremes)

#extremes.to_netcdf(os.path.join(newnc_dir,f'{watershed}dailymeanEXTREMES.nc'))
#print(extremes)

#%%

""" Remove 1979 and 2021 to align with MERRA2 Data Availability """
import xarray as xr
import os

watershed_name = "UpperYuba"

filedir = 'I:\\Emma\\FIROWatersheds\\Data\\GridmetWatershed\\'
# file = f'{watershed_name}dailymeanEXTREMES.nc'
# filepath = os.path.join(filedir,file)
# NCdata = xr.open_dataarray(filepath)
# NCdata_reduced = NCdata[365:-351]
# print(NCdata)
# print(NCdata_reduced)
# NCdata_reduced.to_netcdf(os.path.join(filedir,f'{watershed_name}dailymeanEXTREMES_1980_2021.nc'))


file2 = f'{watershed_name}dailymeansTOTAL.nc'
filepath2 = os.path.join(filedir,file2)
NCdata = xr.open_dataarray(filepath2)
NCdata_reduced = NCdata[365:-351]
print(NCdata)
print(NCdata_reduced)
NCdata_reduced.to_netcdf(os.path.join(filedir,f'{watershed_name}dailymeansTOTAL_1980_2021.nc'))
