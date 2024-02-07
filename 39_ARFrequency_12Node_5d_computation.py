# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 10:31:15 2024

@author: emrus2

Calculates the AR frequency for each SOM node 
- AR catalog data only goes to 2020 - need more data for the result to be accurate
"""
#%% IMPORT PACKAGES
import numpy as np
import scipy
import os

#%% IMPORT DATA
# ar presence data
folder = 'I:\\Emma\\FIROWatersheds\\Data\\'
arpresence = np.load(os.path.join(folder,'ARPresence_90_Percentile_Days_kidmap.npy')) # ar true/false data
    # includes dates (col 1) and ar presence (col 2 - booleans)

# som assignment data
numpatterns = 12
soms = scipy.io.loadmat(os.path.join(folder,f'SOMs\\SomOutput\\IVT_90_{numpatterns}sompatterns_5d.mat')) 
assignment = soms['assignment'][:,0] # som assignment for each day
pat_freq = list(np.squeeze(soms['pat_freq'])) # pattern frequencies
# extreme days and som assignment
extremedays_assign = np.load(os.path.join(folder,'90Percentile_ExtremeDays_5d{numpatterns}NodeAssign.npy'),allow_pickle=True)

#%% CLUSTER ASSIGNED PATTERNS AND CALCULATE AR FREQUENCY
# create array to store 12 SOM ar presence true values
soms_arpresence = []

# loop through all 12 som patterns
for som in range(numpatterns):
    # list to store true ar vals for each som
    artrue = []
    for i,date_arbool in enumerate(arpresence):
        # if assignment is som
        if assignment[i] == som + 1:
            # if ar boolean is true (boolean in 2nd column [1])
            if date_arbool[1] == 'True':
                artrue.append(date_arbool[1])
    soms_arpresence.append(artrue)
    
# calculate ar count for each som using length of ar presence lists
soms_arcount = [len(soms_arpresence[n]) for n, j in enumerate(soms_arpresence)]
# calculate ar frequency by dividing by total number of days assigned to som
soms_arfreq = [i / j * 100 for i, j in zip(soms_arcount, pat_freq)]

#%% SAVE FILES
np.save(os.path.join(folder,f'ARFrequencies_{numpatterns}node_5d_to2021.npy'),soms_arfreq)
