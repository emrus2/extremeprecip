# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 10:31:15 2024

@author: emrus2

Calculates the AR frequency for each SOM node 
- AR catalog data only goes to 2014 - need more data for the result to be accurate
"""
#%% IMPORT PACKAGES
import numpy as np
import scipy
import os

#%% IMPORT DATA
folder = 'I:\\Emma\\FIROWatersheds\\Data\\'
numpatterns = 12
arpresence = np.load(os.path.join(folder,'ARPresence_90_Percentile_Days.npy')) # ar true/false data
soms = scipy.io.loadmat(os.path.join(folder,f'SOMs\\SomOutput\\IVT_90_{numpatterns}sompatterns_5d.mat')) 
assignment = soms['assignment'][:,0] # som assignment for each day
pat_freq = list(np.squeeze(soms['pat_freq'])) # pattern frequencies

# combine extreme days and assignments and save
# extremedays = np.load(os.path.join(folder,'90Percentile_ExtremeDays.npy'),allow_pickle=True) # extreme dates array
# extremedays_assign = np.array(list(zip(extremedays,assignment)))
# np.save(os.path.join(folder,'90Percentile_ExtremeDays_5d12NodeAssign.npy'),extremedays_assign)

# extreme days and som assignment
extremedays_assign = np.load(os.path.join(folder,'90Percentile_ExtremeDays_5d12NodeAssign.npy'),allow_pickle=True)


#%% CLUSTER ASSIGNED PATTERNS AND CALCULATE AR FREQUENCY
# create array to store 12 SOM ar presence true values
soms_arpresence = []

# loop through all 12 som patterns
for som in range(numpatterns):
    # list to store true ar vals for som days
    artrue = []
    for day,artag in enumerate(arpresence):
        # if assignment is som
        if assignment[day] == som + 1:
            # if ar is present
            if artag[1] == 'True':
                artrue.append(artag[1])
    soms_arpresence.append(artrue)
    
# calculate ar count for each som using length of ar presence lists
soms_arcount = [len(soms_arpresence[n]) for n, j in enumerate(soms_arpresence)]
# calculate ar frequency by dividing by total number of days assigned to som
soms_arfreq = [i / j * 100 for i, j in zip(soms_arcount, pat_freq)]

#%% SAVE FILES
np.save(os.path.join(folder,'ARFrequencies_{numpatterns}node_5d.npy'),soms_arfreq)