# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 13:49:49 2023

@author: emrus2

Denormalizes the SOM output that was initially normalized in the data prep section

Saved Files:
   'I:\\Emma\\FIROWatersheds\\Data\\SOMs\\SomOutput' 
        f'IVT_90_12sompatterns_5d_denormalized'

UPDATED 6/12/2023
"""
#%% IMPORT MODULES
#import netCDF4 as nc #if this doesn't work, try to reinstall in anaconda prompt using
    #conda install netcdf4
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
from mpl_toolkits.basemap import Basemap #installed using 
    #conda install -c anaconda basemap
import scipy.io

#%%
#define watershed and directory for plotting
watershed = 'UpperYuba'
ws_directory = f'I:\\Emma\\FIROWatersheds\\Data\\WatershedShapefiles\\California\\{watershed}\\'
#%% IMPORT SOM DATA

# define metvar
metvar = 'IVT'
numpatterns = 12
percentile = 90
daysprior = 2
clusters = daysprior + 1

# change directory and import SOM data from .mat file
mat_dir='I:\\Emma\\FIROWatersheds\\Data\\SOMs\\SomOutput'
os.chdir(mat_dir)
soms = scipy.io.loadmat(f'IVT_{percentile}_{numpatterns}sompatterns_{clusters}d.mat')
patterns = soms['pats']
pat_prop = np.squeeze(soms['pat_prop'])
pat_freq = np.squeeze(soms['pat_freq'])
gridlat = np.squeeze(soms['lat'])
gridlon = np.squeeze(soms['lon'])

#%% LOAD DAILY STD DEVIATION DATA

#%% CONVERT DATA TO UNWEIGHTED VALUES BY LATITUDE
# define original weights computation
latweights = np.sqrt(np.cos((np.radians(gridlat))))
# identify final length of 1d
newpatterns = np.empty_like(patterns)
# loop through each som pattern (5day)
for i, arr in enumerate(patterns):
    print(i,arr)
    # reshape to 3d
    merrareduced = arr.reshape(clusters,len(gridlat),len(gridlon))
    # unweight by LATITUDE
    for n in range(merrareduced.shape[1]): #lats
        merrareduced[:,n,:] /= latweights[n]
    # reshape to 1d, clusters
    merraflat = merrareduced.reshape(-2) # remove 2 dimensions
    newpatterns[i] = merraflat

#%% SAVE NEW PATTERNS
np.save(os.path.join(mat_dir,f'IVT_{percentile}_{numpatterns}sompatterns_{clusters}d_denormalized'),newpatterns)