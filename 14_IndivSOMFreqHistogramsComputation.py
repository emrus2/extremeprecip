# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 13:49:49 2023

@author: emrus2

Calculates the monthly frequency of each node, saves as:
    
    I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_{numpatterns}NodeAssignMonthlyHistogram.npy

UPDATED 6/12/2023
"""
#%% IMPORT MODULES
import numpy as np
import os
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
import scipy.io
from datetime import datetime
#%% IMPORT SOM DATA
# change directory and import SOM data from .mat file
mat_dir='I:\\Emma\\FIROWatersheds\\Data\\SOMs\\SomOutput'
numpatterns = 12
percentile = 90
clusters = 5

os.chdir(mat_dir)
soms = scipy.io.loadmat(f'IVT_{percentile}_{numpatterns}sompatterns_{clusters}d.mat')
pats = np.squeeze(soms['pats'])
asn_err = np.squeeze(soms['assignment'])
assignment = asn_err[:,0]

if clusters == 1:
    extremedays = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_ExtremeDays.npy')
    assign_day = np.stack((assignment,extremedays),axis=1)
    #np.save('I:\\Emma\\FIROWatersheds\\Data\\ExtremeDaysandNodeAssign.npy',assign_day)
    # np.save(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_ExtremeDaysand{numpatterns}NodeAssign.npy',assign_day)
else:
    # import extreme days list
    extremeevents = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\ExtremeEvents_{clusters}d.npy', \
                                 allow_pickle=True)
    # reduce to only last day (actual extreme day)
    extremedays = extremeevents[:,-1]
    # convert to list of strings
    extremedays = [datetime.strftime(day,'%Y%m%d') for day in extremedays]
    
#%% CLUSTER ASSIGNED PATTERNS AND CALCULATE AVERAGE
months = ['10','11','12','01','02','03']
numpatterns = len(pats[:,0])

# create array to store 12 SOM composites
node_monthfreq = np.zeros((numpatterns,6))
# loop through all 12 som patterns
for som in range(len(pats[:,0])):
    monthlist = []
    for i,day in enumerate(assignment):
        date = extremedays[i]        
        if assignment[i] == som+1:
            month = date[4:6]
            monthlist.append(month)
    print(monthlist)
    for count,monthstr in enumerate(months):
        monthcount = monthlist.count(monthstr)
        node_monthfreq[som,count] = monthcount

np.save(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_{numpatterns}NodeAssignMonthlyHistogram_{clusters}d.npy',node_monthfreq)