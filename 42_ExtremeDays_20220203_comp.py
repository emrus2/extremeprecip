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
from datetime import datetime,timedelta

#%% LOAD 90TH PERCENTILE THRESHOLD
savepath = 'I:\\Emma\\FIROWatersheds\\Data\\'
threshold = np.load(os.path.join(savepath,'90percentilethreshold.npy'))
threshold = float(threshold)

#%% OPEN 2022-2023 DATA
years = [2023,2024]
gridmetpath = 'I:\\Emma\\FIROWatersheds\\Data\\GRIDMET\\'
pr_22 = list(np.load(os.path.join(gridmetpath,f'UpperYubadailymeanprecip_{str(years[0])}.npy')))
pr_23 = list(np.load(os.path.join(gridmetpath,f'UpperYubadailymeanprecip_{str(years[1])}.npy')))
# combine pr data
pre_2223 = pr_22 + pr_23
# write dates data
dates = [datetime.strftime(datetime(years[0],1,1)+timedelta(days=n),"%Y%m%d") for n in range(len(pre_2223))]
# write months data
precipmonths = [item[4:6] for item in dates]

# combine dates,precip, and months data
precip = np.array(list(zip(dates,pre_2223,precipmonths)))
# reduce total data to wet season data (October-March)
winter = ['01', '02','03', '10', '11', '12']
precipwinter =  np.array([precip[i] for i in range(len(precip)) if precip[i][2] in winter])
# reduce data to those with nonzero values
precipreduced = np.array([precipwinter[i] for i in range(len(precipwinter)) if precipwinter[i][1]])
#Calculate extreme days as those exceeding 90th percentile
extremes_winter = np.array([precipreduced[i] for i in range(len(precipreduced)) if float(precipreduced[i][1]) > threshold])
# remove month data
extremes_winter = extremes_winter[:,0:2]

#save values to numpy file
np.save(os.path.join(savepath,'ExtremePrecipitationDays_2022_23.npy'),extremes_winter)
