# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 13:49:49 2023

@author: emrus2

Calculates the annual composition of node patterns

UPDATED 6/12/2023
"""
#%% IMPORT MODULES
import numpy as np
import os
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
import scipy.io
#%% IMPORT SOM DATA
# change directory and import SOM data from .mat file
mat_dir='I:\\Emma\\FIROWatersheds\\Data\\SOMs\\SomOutput'
os.chdir(mat_dir)
numpatterns = 12
percentile = 90
soms = scipy.io.loadmat(f'IVT_{percentile}_{numpatterns}_sompatterns.mat')
pats = np.squeeze(soms['pats'])
asn_err = np.squeeze(soms['assignment'])
assignment = asn_err[:,0]

extremedays = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_ExtremeDays.npy')

assign_day = np.stack((assignment,extremedays),axis=1)
#np.save('I:\\Emma\\FIROWatersheds\\Data\\ExtremeDaysandNodeAssign.npy',assign_day)
    
#%% CLUSTER ASSIGNED PATTERNS AND CALCULATE AVERAGE
years = np.arange(1980,2022,1)

# create array to store 12 SOM composites
node_yearfreq = np.zeros((numpatterns,len(years)))
# loop through all 12 som patterns
for som in range(numpatterns):
    # list to store years where som is assigned
    yearlist = []
    for i,day in enumerate(assignment):
        # if day is assigned to som
        if assignment[i] == som+1:
            date = extremedays[i] 
            year = date[0:4]
            # add year to yearlist
            yearlist.append(year)
    yearlistint =  [int(y) for y in yearlist]
    print(yearlistint)
    # count the years in the yearlist
    for j,yearstr in enumerate(years):
        yearcount = yearlistint.count(yearstr)
        print(yearcount)
        node_yearfreq[som,j] = yearcount
        
# create total across years
totalnodesyear = np.sum(node_yearfreq,axis=0)

np.save(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}_{numpatterns}NodeAssignYearlyHistogram.npy',node_yearfreq)
np.save(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}_{numpatterns}NodeAssignTotalYearly.npy',totalnodesyear)
