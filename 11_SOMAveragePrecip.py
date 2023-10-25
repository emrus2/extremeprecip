# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 13:49:49 2023

@author: emrus2

Calculates the precipitation values for days assigned to each node
Calculates the average precipitation of days assigned to each node
Saves two .npy files, one with all precip and one with average precip

Saved files:
    'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_{numpatterns}NodeAssignAvergagePrecip.npy'
    'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_{numpatterns}NodeAssignAllPrecip.npy'

UPDATED 6/12/2023
"""
#%% IMPORT MODULES
import numpy as np
import os
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
import scipy.io

#%% IMPORT SOM DATA
numpatterns = 9
percentile = 90
clusters = 5

# change directory and import SOM data from .mat file
mat_dir='I:\\Emma\\FIROWatersheds\\Data\\SOMs\\SomOutput'
os.chdir(mat_dir)
soms = scipy.io.loadmat(f'IVT_{percentile}_{numpatterns}sompatterns_{clusters}d.mat')
asn_err = np.squeeze(soms['assignment'])
assignment = asn_err[:,0]

#%% IMPORT EXTREME PRECIPITATION DATA
precip = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_ExtremeDaysPrecip.npy')
    
#%% CLUSTER ASSIGNED PATTERNS AND CALCULATE AVERAGE PRECIP
# create array to store 12 SOM average precip
soms_averageprecip = []
soms_medianprecip = []
soms_allprecip = []
# loop through all 12 som patterns
for som in range(numpatterns):
    # list to store precip vals for som days
    somprecip = []
    for day,precipval in enumerate(precip):
        if assignment[day] == som + 1:
            somprecip.append(precipval)
    soms_allprecip.append(somprecip)
    som_average = np.mean(somprecip)
    som_median = np.median(somprecip)
    soms_averageprecip.append(som_average)
    soms_medianprecip.append(som_median)
    
# create numpy array for all precipitation values
arr_size = (max(len(lst) for lst in soms_allprecip))
allprecip_arr = np.empty((numpatterns,arr_size))
allprecip_arr[:] = np.nan
for i,j in enumerate(soms_allprecip):
    allprecip_arr[i][0:len(j)] = j
    
#%% SAVE DATA
os.chdir('I:\\Emma\\FIROWatersheds\\Data')
np.save(f'{percentile}Percentile_{numpatterns}NodeAssignAvergagePrecip_{clusters}d.npy',soms_averageprecip)
np.save(f'{percentile}Percentile_{numpatterns}NodeAssignMedianPrecip_{clusters}d.npy',soms_medianprecip)
np.save(f'{percentile}Percentile_{numpatterns}NodeAssignAllPrecip_{clusters}d.npy',allprecip_arr)
