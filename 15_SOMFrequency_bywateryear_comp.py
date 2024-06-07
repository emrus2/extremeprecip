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
years = np.arange(1979,2022,1)

# create array to store 12 SOM composites
node_yearfreq = np.zeros((numpatterns,len(years)))
# loop through all 12 som patterns
for som in range(numpatterns):
    # list to store years where som is assigned
    datelist = []
    for i,day in enumerate(assignment):
        # if day is assigned to som
        if assignment[i] == som+1:
            date = extremedays[i] 
            yearmonth = date[0:6]
            #year = date[0:4]
            # add year to yearlist
            datelist.append(yearmonth)
            
    datelistint =  [int(y) for y in datelist]
    print(datelistint)
    # count the years in the yearlist
    yearcounts = []
    for count,year in enumerate(years):
        # define water year starts and ends
        yearstart = int(str(year) + '10')
        yearend = int(str(year+1) + '03')
        # count occurrence for water years (oct-mar)
        wycount = 0
        for day in datelistint:
            if day in range(yearstart,yearend+1):
                wycount +=1
            print(wycount)
            yearcounts.append(wycount)
        #yearcount = yearlistint.count(yearstr)
        #print(yearcount)
        node_yearfreq[som,count] = wycount
        
# create total across years
totalnodesyear = np.sum(node_yearfreq,axis=0)

np.save(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}_{numpatterns}NodeAssignYearlyHistogram_wateryear.npy',node_yearfreq)
np.save(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}_{numpatterns}NodeAssignTotalYearly_wateryear.npy',totalnodesyear)
