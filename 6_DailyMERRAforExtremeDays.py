# -*- coding: utf-8 -*-
# file_name.py
"""
Created on Fri December 9 15:27:00 2022

@author: emrus2

Compares MERRA-2 file names to extreme precipitation dates
Creates file of MERRA-2 extreme precipitation days

Saved files:
    'MERRA2_Z500_Yuba_Extremes_Daily_1980-2021_WINTERDIST.nc'
    'MERRA2_SLP_Yuba_Extremes_Daily_1980-2021_WINTERDIST.nc'
    'MERRA2_300W_Yuba_Extremes_Daily_1980-2021_WINTERDIST.nc'
    'MERRA2_850T_Yuba_Extremes_Daily_1980-2021_WINTERDIST.nc'
    'MERRA2_IVT_Yuba_Extremes_Daily_1980-2021_WINTERDIST.nc'
    'MERRA2_Z850_Yuba_Extremes_Daily_1980-2021_WINTERDIST.nc'

UPDATED 4/30/2023

#still need to get a few files foor Z850
"""
#%% IMPORT EXTREME DATES
#import desired modules
import os
import numpy as np
import xarray as xr
import glob
import netCDF4 as nc
#import numpy as np
from datetime import datetime, timedelta

#import extreme values from GRIDMET precip data
watershed = "UpperYuba"
percentile = 90
gridmetpath = f'I:\\Emma\\FIROWatersheds\\Data\\Gridmet\\{watershed}\\Extremes'
gridmetfile = f"{watershed}DailymeanWINTEREXTREMES{percentile}_1980_2021.nc"
gridmet = os.path.join(gridmetpath,gridmetfile)
extremes = nc.Dataset(gridmet,mode='r')
print(extremes) #7655 in length
print(extremes.dimensions)
print(extremes.variables)
precipdays = extremes.variables['day'][:]
precip = extremes.variables['precipitation_amount'][:] #time x height x lat x lon
extremes.close()

#convert datetime to matching string of MERRA2 files
datestr = [datetime.strftime(datetime(1979,1,1)+timedelta(days=float(precipdays[n])),"%Y%m%d") for n in range(precip.shape[0])]

#create list of datetimes where precip is extreme
extremedates = []
precipvalues = []
#extremes are all where precip is not NaN in dataset
for n in range(precip.shape[0]):
    if precip[n]:
        extremedates.append(datestr[n])
        precipvalues.append(precip[n])
        
#np.save(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_ExtremeDays.npy',np.array(extremedates))
#np.save(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_ExtremeDaysPrecip.npy',np.array(precipvalues))

        #should be equal to 130 for 95th, 260 for 90th
#remove extreme days in summer (keep OCT-MAR)
extremestr = [date for date in extremedates if int(date[5]) not in range(4,10)]
    #should be equal to 260

#%% DEFINE EXTREME MERRA2 FILES
#define files of interest
#metvar = input('Enter MERRA-2 Variable: Z500, SLP, 300W, 850T, Z850 or IVT \n')
#metvars = ['Z500', 'SLP', '300W', '850T', 'Z850','IVT','850W']
metvars = ['850W']
for metvar in metvars:
    metpath = {'Z500':'500_hPa_Geopotential_Height_3hourly','SLP':'Sea_level_pressure_3hourly','IVT':'IVT_daily', \
               '300W':'East_and_North_wind_components_at_300_hPa','850T':'Temperature_at_850_hPa_3hourly', \
                   'Z850':'850_hPa_Geopotential_Height_3hourly','850W':'East_and_North_wind_components_at_850_hPa'}
    folderpath = f'I:\\MERRA2\\Daily_and_Subdaily\\{metpath[metvar]}'
    
    extremefiles = []
    extremefiles2 = []
    
    #define files of interest, in this case 1980-2021 years for ALL days
    if metvar == 'IVT':
        files = glob.glob(os.path.join(folderpath,"MERRA2.tavg1_2d_int_Nx.1*")) + glob.glob(os.path.join(folderpath,'MERRA2.tavg1_2d_int_Nx.200*')) + glob.glob(os.path.join(folderpath,'MERRA2.tavg1_2d_int_Nx.201[0-6]*'))
        files2 = glob.glob(os.path.join(folderpath,'MERRA2.tavg1_2d_int_Nx.201[7-9]*')) + glob.glob(os.path.join(folderpath,'MERRA2.tavg1_2d_int_Nx.202[0-1]*'))
    elif metvar == '850T' or metvar == '850W':
        files = glob.glob(os.path.join(folderpath,"MERRA2*.inst3_3d_asm_Np.1*.SUB.nc")) + glob.glob(os.path.join(folderpath,'MERRA2*.inst3_3d_asm_Np.20[0-1]*.SUB.nc')) + glob.glob(os.path.join(folderpath,'MERRA2*.inst3_3d_asm_Np.202[0-1]*.SUB.nc'))
    elif metvar == 'Z850':
        files = glob.glob(os.path.join(folderpath,"MERRA2.inst3_3d_asm_Np.1*.SUB.nc")) + glob.glob(os.path.join(folderpath,'MERRA2.inst3_3d_asm_Np.20[0-1]*.SUB.nc')) + glob.glob(os.path.join(folderpath,'MERRA2.inst3_3d_asm_Np.202[0-1]*.SUB.nc')) + glob.glob(os.path.join(folderpath,'MERRA2*.inst3_3d_asm_Np.202109*.SUB.nc'))
    else:
        files = glob.glob(os.path.join(folderpath,"MERRA2.inst3_3d_asm_Np.1*.SUB.nc")) + glob.glob(os.path.join(folderpath,'MERRA2.inst3_3d_asm_Np.20[0-1]*.SUB.nc')) + glob.glob(os.path.join(folderpath,'MERRA2.inst3_3d_asm_Np.202[0-1]*.SUB.nc'))
         #files length should be equal to 15341
    
    #compare list of strings to string of file
    for file in files:
        #add extreme files to list
        if (any(date in file for date in extremestr)):
            extremefiles.append(file)
            #extremefiles should be 130 in length (112 for IVT), 260 for 90th percentile
    if metvar == 'IVT':
        for file in files2:
            if (any(date in file for date in extremestr)):
                extremefiles2.append(file)
                #extremefiles2 should be 18 in length for IVT
            
    #%% CREATE EXTREMES DATASET
    #create dataset of extreme files
    if metvar == 'IVT':
        ds1 = xr.open_mfdataset(extremefiles, combine='nested', concat_dim='time')
        ds2 = xr.open_mfdataset(extremefiles2, combine='nested',concat_dim='time')
        ds1_red = ds1.drop_vars(['UFLXQL','VFLXQL'])
        ds = xr.concat([ds1_red,ds2],dim='time')
    else:
        ds = xr.open_mfdataset(extremefiles, combine='nested',concat_dim='time')
        #calculate daily average merra data
        ds = ds.groupby('time.date').mean() #this works great!!
        ds.squeeze() #should result in 162(260) dates
        ds['date'] = ds['date'].astype('datetime64')
        
    print(ds)
    
    #calculate composite mean for merra
    #total_means = np.squeeze(ds.mean(dim='time'))
    #print(total_means) #should be 1 time dimension
    
    #%% SAVE EXTREMES DATASET
    #save as file
    #define new files directory and names
    os.chdir(f'I:\\Emma\\FIROWatersheds\\Data\\DailyMERRA2\\{metvar}')
    savefile = f'MERRA2_{metvar}_Yuba_Extremes{percentile}_Daily_1980-2021_WINTERDIST.nc'
    #savefile2 = f'{metvar}_extremes_composite.nc'
    ds.to_netcdf(savefile)
    #total_means.to_netcdf(savefile2)

#%%
"""
For entire year climatology extremes
"""
# #import desired modules
# import os
# import xarray as xr
# import glob
# import netCDF4 as nc
# #import numpy as np
# from datetime import datetime, timedelta

# #import extreme values from GRIDMET precip data
# watershed = "UpperYuba"
# gridmetpath = f'I:\\Emma\\FIROWatersheds\\Data\\Gridmet\\{watershed}\\Extremes'
# gridmetfile = "DailymeanEXTREMES_1980_2021.nc"
# gridmet = os.path.join(gridmetpath,gridmetfile)
# extremes = nc.Dataset(gridmet,mode='r')
# print(extremes)
# print(extremes.dimensions)
# print(extremes.variables)
# precipdays = extremes.variables['day'][:]
# precip = extremes.variables['precipitation_amount'][:] #time x height x lat x lon
# extremes.close()

# #convert datetime to matching string of MERRA2 files
# datestr = [datetime.strftime(datetime(1979,1,1)+timedelta(days=float(precipdays[n])),"%Y%m%d") for n in range(precip.shape[0])]

# #create list of datetimes where precip is extreme
# extremestrall = []
# #extremes are all where precip is not NaN in dataset
# for n in range(precip.shape[0]):
#     if precip[n]:
#         extremestrall.append(datestr[n])
#         #should be equal to 170
# #remove extreme days in summer (keep OCT-MAR)
# extremestr = [date for date in extremestrall if int(date[5]) not in range(4,10)]
#     #should be equal to 162

# #%% DEFINE EXTREME MERRA2 FILES
# #define files of interest
# metvar = input('Enter MERRA-2 Variable: Z500, SLP, 300W, 850T, or IVT \n')
# metpath = {'Z500':'500_hPa_Geopotential_Height_3hourly','SLP':'Sea_level_pressure_3hourly','IVT':'IVT_daily','300W':'East_and_North_wind_components_at_300_hPa','850T':'Temperature_at_850_hPa_3hourly'}
# folderpath = f'I:\\MERRA2\\Daily_and_Subdaily\\{metpath[metvar]}'

# extremefiles = []
# extremefiles2 = []

# #define files of interest, in this case 1980-2021 years for ALL days
# if metvar == 'IVT':
#     files = glob.glob(os.path.join(folderpath,"MERRA2.tavg1_2d_int_Nx.1*")) + glob.glob(os.path.join(folderpath,'MERRA2.tavg1_2d_int_Nx.200*')) + glob.glob(os.path.join(folderpath,'MERRA2.tavg1_2d_int_Nx.201[0-6]*'))
#     files2 = glob.glob(os.path.join(folderpath,'MERRA2.tavg1_2d_int_Nx.201[7-9]*')) + glob.glob(os.path.join(folderpath,'MERRA2.tavg1_2d_int_Nx.202[0-1]*'))
# elif metvar == '850T':
#     files = glob.glob(os.path.join(folderpath,"MERRA2*.inst3_3d_asm_Np.1*.SUB.nc")) + glob.glob(os.path.join(folderpath,'MERRA2*.inst3_3d_asm_Np.20[0-1]*.SUB.nc')) + glob.glob(os.path.join(folderpath,'MERRA2*.inst3_3d_asm_Np.202[0-1]*.SUB.nc'))
# else:
#     files = glob.glob(os.path.join(folderpath,"MERRA2.inst3_3d_asm_Np.1*.SUB.nc")) + glob.glob(os.path.join(folderpath,'MERRA2.inst3_3d_asm_Np.20[0-1]*.SUB.nc')) + glob.glob(os.path.join(folderpath,'MERRA2.inst3_3d_asm_Np.202[0-1]*.SUB.nc'))
#      #files length should be equal to 15341

# #compare list of strings to string of file
# for file in files:
#     #add extreme files to list
#     if (any(date in file for date in extremestr)):
#         extremefiles.append(file)
#         #extremefiles should be 162 in length (149 for IVT)
# if metvar == 'IVT':
#     for file in files2:
#         if (any(date in file for date in extremestr)):
#             extremefiles2.append(file)
#             #extremefiles2 should be 21 in length for IVT
        
#  #%% CREATE EXTREMES DATASET
# #create dataset of extreme files
# if metvar == 'IVT':
#     ds1 = xr.open_mfdataset(extremefiles, combine='nested', concat_dim='time')
#     ds2 = xr.open_mfdataset(extremefiles2, combine='nested',concat_dim='time')
#     ds1_red = ds1.drop_vars(['UFLXQL','VFLXQL'])
#     ds = xr.concat([ds1_red,ds2],dim='time')
# else:
#     ds = xr.open_mfdataset(extremefiles, combine='nested',concat_dim='time')
#     #calculate daily average merra data
#     ds = ds.groupby('time.date').mean() #this works great!!
#     ds.squeeze() #should result in 162 dates
#     ds['date'] = ds['date'].astype('datetime64')
    
# print(ds)

# #calculate composite mean for merra
# #total_means = np.squeeze(ds.mean(dim='time'))
# #print(total_means) #should be 1 time dimension

# #%% SAVE EXTREMES DATASET
# #save as file
# #define new files directory and names
# os.chdir(f'I:\\Emma\\FIROWatersheds\\Data\\DailyMERRA2\\{metvar}')
# savefile = f'MERRA2_{metvar}_Yuba_Extremes_Daily_1980-2021_OCT-MAR.nc'
# #savefile2 = f'{metvar}_extremes_composite.nc'
# ds.to_netcdf(savefile)
# #total_means.to_netcdf(savefile2)