# -*- coding: utf-8 -*-
"""
Created on Wed May 10 15:48:01 2023

@author: emrus2

Creates a climatology of wet season daily mean precipitaiton in Upper Yuba Watershed
Also calculates the extreme precipitation days depending on desired percentile
Has been used for 95th and 90th percentile 

Saved files:
    'UpperYubaDailymeansWINTERTOTAL{percentile}_1980_2021'
    'UpperYubaDailymeanWINTEREXTREMES{percentile}_1980_2021'

"""
#import required packages
import os
import xarray as xr

#define filepath and file of daily mean total precipitation
watershed = "UpperYuba"
folderpath = f'I:\\Emma\\FIROWatersheds\\Data\\Gridmet\\{watershed}\\DailyMeans'
filename = f'{watershed}dailymeansTOTAL_1980_2021.nc'
filepath = os.path.join(folderpath, filename)

#open filepath with xarray
NCdatatotal = xr.open_dataarray(filepath)
print(NCdatatotal)

#reduce total data to wet season data (October-March)
winter = [1, 2, 3, 10, 11, 12]
NCdatawinter =  NCdatatotal.where(NCdatatotal['day.month'].isin(winter), drop=True)
print(NCdatawinter)

percentile = 90
#calculate 95th percentile threshold
percent_95 = NCdatawinter.quantile(percentile/100)
print(percent_95)

#Calculate extreme days as those exceeding 95th percentile
extremes_winter = NCdatawinter.where(NCdatawinter>percent_95)
print(extremes_winter)
NCdatawinter.plot.line(marker='o',lw=0,markersize=2)
extremes_winter.plot.line(marker='o',lw=0,markersize=2)

#save values to nc files
newfolder = f'I:\\Emma\\FIROWatersheds\\Data\\Gridmet\\{watershed}\\Extremes'
#NCdatawinter.to_netcdf(os.path.join(folderpath,f'{watershed}DailymeansWINTERTOTAL_1980_2021.nc'))
extremes_winter.to_netcdf(os.path.join(newfolder,f'{watershed}DailymeanWINTEREXTREMES{percentile}_1980_2021.nc'))

# #%%
# #COLLECT VARIABLE DATA FROM NC FILE
# #open the netcdf file in read mode
# gridfile = nc.Dataset(filepath,mode='r')
# print(gridfile)
# print(gridfile.variables['day'])
# day = gridfile.variables['day'][:]
# #days are days since 01/01/1979
# precip = gridfile.variables['precipitation_amount'][:]
# #define datetime for associated precipitations
# dates = np.array([datetime.strftime(datetime(1979,1,1)+timedelta(days=float(day[n])),"%Y%m%d") for n in range(precip.shape[0])])
# #length is 15341

# # find indices for wet season days
# winter_ind = []
# for i in range(len(dates)):
#     date = dates[i]
#     if int(date[5]) not in range (4,10):
#            winter_ind.append(i)
           
# #reduce precip and date data to winter indices
# precip_winter = precip[winter_ind]
# day_winter = day[winter_ind]
# dates_winter = dates[winter_ind]

# #define precipitation array as xarray
# precip_winterxr = xr.DataArray(precip_winter)     
      
# #calculate 95th percentile threshold
# percent_95 = precip_winterxr.quantile(0.95)
# print(percent_95)

# #Calculate extreme days as those exceeding 95th percentile
# extremes_winterxr = precip_winterxr.where(precip_winterxr>percent_95)
# print(extremes_winterxr)
# precip_winterxr.plot.line(marker='o',lw=0,markersize=2)
# extremes_winterxr.plot.line(marker='o',lw=0,markersize=2)

# #save values to nc files
# precip_winterxr.to_netcdf(os.path.join(folderpath,f'{watershed}dailymeansWINTERTOTAL_1980_2021.nc'))
# extremes_winterxr.to_netcdf(os.path.join(folderpath,f'{watershed}dailymeanWINTEREXTREMES_95_1980_2021.nc'))