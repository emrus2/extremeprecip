# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 13:49:49 2023

@author: emrus2

Calculates the succession of patterns in each node, saves as .npy file
    
    I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_{numpatterns}NodeAssignMonthlyHistogram.npy

UPDATED 6/12/2023
"""
#%% IMPORT MODULES
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
import scipy.io
#%% IMPORT SOM DATA
metvar = 'IVT'
numpatterns = 9
percentile = 90

# change directory and import SOM data from .mat file
mat_dir='I:\\Emma\\FIROWatersheds\\Data\\SOMs\\SomOutput'
os.chdir(mat_dir)
soms = scipy.io.loadmat(f'{metvar}_{percentile}_{numpatterns}_sompatterns.mat')
pat_prop = np.squeeze(soms['pat_prop']) #number of patterns in each SOM
#%% CATEGORIZE EXTREME DAYS INTO EVENTS BASED ON SUCCESSION

# import extreme days and node assignment array
nodeassign = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_ExtremeDaysand{numpatterns}NodeAssign.npy')
# import extreme days and precipitation amount
extremeprecip = np.load(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_ExtremeDaysPrecip.npy')        
extremeprecip = extremeprecip.reshape(len(extremeprecip),1)

# combine two arrays
dateassign = np.concatenate((nodeassign,extremeprecip), axis=1)

n = 0 # event duration counter
# create empty lists to store all event information
alleventdates = []
alleventnodes = []
alleventprecip = []
# create empty lists to store event information
eventdates = []
eventnodes = []
eventprecip = []
# loop through each extreme day
for i, data in enumerate(dateassign):
    assignment = float(data[0])
    assignment_prev = float(dateassign[i-1,0])
    date = int(data[1])
    date_prev = int(dateassign[i-1,1])
    precip = float(data[2])
    precip_prev = float(dateassign[i-1,2])
    # if the event continues:
    if date == date_prev + 1:
        n += 1
        # add event information to event lists
        eventdates.append(date_prev)
        eventdates.append(date)
        eventnodes.append(assignment_prev)
        eventnodes.append(assignment)
        eventprecip.append(precip_prev)
        eventprecip.append(precip)
        # remove repeating dates
        if len(eventdates) >= 2:
            eventdates.pop()
            eventnodes.pop()
            eventprecip.pop()
    # if the event has ended (does not continue):
    else:
        # add event information to event lists
        eventdates.append(date_prev)
        eventnodes.append(assignment_prev)
        eventprecip.append(precip_prev)
        alleventdates.append(eventdates)
        alleventnodes.append(eventnodes)
        alleventprecip.append(eventprecip)
        # reset event duration counter and event lists
        n = 0
        eventdates = []
        eventnodes = []
        eventprecip = []

'''these have been saved'''
#np.save(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_{numpatterns}NodeClusteredEvents.npy',alleventdates)
#np.save(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_{numpatterns}NodeClusteredEventsNodes.npy',alleventnodes)
#np.save(f'I:\\Emma\\FIROWatersheds\\Data\\{percentile}Percentile_{numpatterns}NodeClusteredEventsPrecip.npy',alleventprecip)