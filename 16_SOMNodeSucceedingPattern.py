# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 13:49:49 2023

@author: emrus2

Here we will bin the days assigned to each IVT SOM together and then plot the
composite patterns of other variables, like SLP, Z500Anom, etc. 

UPDATED 6/12/2023
"""
#%% IMPORT MODULES
#import netCDF4 as nc #if this doesn't work, try to reinstall in anaconda prompt using
    #conda install netcdf4
import numpy as np
import os
#import matplotlib.pyplot as plt
#from matplotlib.colors import ListedColormap
#import matplotlib.cm as cm
#import matplotlib.transforms as mtransforms
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
#from mpl_toolkits.basemap import Basemap #installed using 
    #conda install -c anaconda basemap
#import scipy.io
#%% IMPORT SOM DATA
# change directory and import SOM data from .mat file
#mat_dir='I:\\Emma\\FIROWatersheds\\Data\\SOMs\\SomOutput'
#os.chdir(mat_dir)
#soms = scipy.io.loadmat('IVT_sompatterns_rd.mat')
#pats = np.squeeze(soms['pats'])
#asn_err = np.squeeze(soms['assignment'])
#assignment = asn_err[:,0]

assign_day = np.load('I:\\Emma\\FIROWatersheds\\Data\\ExtremeDaysandNodeAssign.npy')
assignment = assign_day[:,0]
#%% CLUSTER ASSIGNED PATTERNS AND CALCULATE AVERAGE
months = ['10','11','12','01','02','03']
numpatterns = 12

#NEED TO MAKE SURE IT IS THE SUCCEEDING DAY, NOT THE NEXT EXTREME DAY.....

# create array to store 12 SOM composites
node_monthfreq = np.zeros((numpatterns,6))
# loop through all 12 som patterns
for som in range(numpatterns):
    succeedlist = []
    for i,day in enumerate(assignment):
        if float(assignment[i]) == som+1:
            succeed = assignment[i+1]
            print(assignment[i+1],som+1)
            monthlist.append(succeed)
    print(monthlist)
    for count,monthstr in enumerate(months):
        monthcount = monthlist.count(monthstr)
        node_monthfreq[som,count] = monthcount

np.save('I:\\Emma\\FIROWatersheds\\Data\\NodeAssignMonthlyHistogram_rd.npy',node_monthfreq)

#%% PLOT HISTOGRAM DATA