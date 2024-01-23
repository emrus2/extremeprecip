# -*- coding: utf-8 -*-
"""
Created on Wed May 10 15:48:01 2023

@author: emrus2

Computing days in the 2022-2023 winter that exceed the 1980-2021 90th
percentile threshold
"""
#%% IMPORT PACKAGES
import xarray as xr
import os
import numpy as np

#%% LOAD 90TH PERCENTILE THRESHOLD
savepath = 'I:\\Emma\\FIROWatersheds\\Data\\'
threshold = np.load(os.path.join(savepath,'90percentilethreshold.npy'))

#%% OPEN 2022-2023 DATA

# create dataset of 2022-2023 data
gridmetpath = 'I:\\GRIDMET\\pr\\'
# pr_22 = os.path.join(gridmetpath,'pr_2022.nc')
# pr_23 = os.path.join(gridmetpath,'pr_2023.nc')
# ds = xr.open_mfdataset((pr_22, pr_23), combine='nested',concat_dim='day') #365x2 = 730

pr_22 = xr.open_dataarray(os.path.join(gridmetpath,'pr_2022.nc'))
pr_23 = xr.open_dataarray(os.path.join(gridmetpath,'pr_2023.nc'))
ds = xr.concat([pr_22,pr_23],dim='day') #365x2 = 730
print(ds)

# datafolder = 'I:\\Emma\\FIROWatersheds\\Data\\Gridmet\\'
# ds.to_netcdf(os.path.join(datafolder,f'{watershed}pr_20222023.nc'))

#reduce total data to wet season data (October-March)
winter = [1, 2, 3, 10, 11, 12]
dswinter =  ds.where(ds['day.month'].isin(winter), drop=True) # 182x2 = 364
print(dswinter)

#Calculate extreme days as those exceeding 90th percentile
extremes_winter = dswinter.where(dswinter>threshold)
print(extremes_winter)
dswinter.plot.scatter()
extremes_winter.plot.scatter()

#save values to nc files
# newfolder = f'I:\\Emma\\FIROWatersheds\\Data\\Gridmet\\{watershed}\\Extremes'
# NCdatawinter.to_netcdf(os.path.join(folderpath,f'{watershed}DailymeansWINTERTOTAL_1980_2021.nc'))
# extremes_winter.to_netcdf(os.path.join(newfolder,f'{watershed}DailymeanWINTEREXTREMES{percentile}_1980_2021.nc'))