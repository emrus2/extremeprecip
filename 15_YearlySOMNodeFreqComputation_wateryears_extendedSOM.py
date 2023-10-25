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
from datetime import datetime
#%% IMPORT SOM DATA
# change directory and import SOM data from .mat file
mat_dir='I:\\Emma\\FIROWatersheds\\Data\\SOMs\\SomOutput'
os.chdir(mat_dir)

# define som type
numpatterns = 9
percentile = 90
clusters = 5

# import SOM data
soms = scipy.io.loadmat(f'IVT_{percentile}_{numpatterns}sompatterns_{clusters}d.mat')
pats = np.squeeze(soms['pats'])
asn_err = np.squeeze(soms['assignment'])
assignment = asn_err[:,0]

# import extreme days list
extremeevents = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\ExtremeEvents_{clusters}d.npy', \
                             allow_pickle=True)
# reduce to only last day (actual extreme day)
extremedays = extremeevents[:,-1]
# convert to list of strings
extremedays = [datetime.strftime(day,'%Y%m%d') for day in extremedays]
    
#%% CLUSTER ASSIGNED PATTERNS AND CALCULATE AVERAGE
years = np.arange(1979,2022,1)

# create array to store 9 SOM composites
node_yearfreq = np.zeros((numpatterns,len(years)))
# loop through all som patterns
for som in range(numpatterns):
    # list to store years where som is assigned
    datelist = []
    # loop through each day's assignment
    for i,day in enumerate(assignment):
        # if day is assigned to som
        if assignment[i] == som+1:
            # find date string
            date = extremedays[i] 
            # identify month and year of date
            yearmonth = date[0:6]
            # add year and month to yearlist
            datelist.append(yearmonth)
    
    # convert string list to integers
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
        print(yearcounts)
        node_yearfreq[som,count] = wycount
        
# create total across years
totalnodesyear = np.sum(node_yearfreq,axis=0)

np.save(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}_{numpatterns}NodeAssignYearlyHistogram_wateryear_{clusters}d.npy',node_yearfreq)
np.save(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}_{numpatterns}NodeAssignTotalYearly_wateryear_{clusters}d.npy',totalnodesyear)
