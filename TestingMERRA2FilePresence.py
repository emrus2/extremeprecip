#%%
# -*- coding: utf-8 -*-
# file_name.py
"""
Created on Fri December 9 15:27:00 2022

@author: emrus2

Compares number of files in folder to Z500, to confirm all files are there
for desired meteorological variable

UPDATED 4/30/2023
"""

#%% DEFINE EXTREME MERRA2 FILES

import numpy as np
import glob
import os

#define files of interest
metvars = ['Z500', 'Z850']
for i in np.arange(202101,202113,1):
    metvar1 = 'Z500'
    metpath1 = '500_hPa_Geopotential_Height_3hourly'
    metvar2 = 'Z850'
    metpath2 = '850_hPa_Geopotential_Height_3hourly'
    folderpath1 = f'I:\\MERRA2\\Daily_and_Subdaily\\{metpath1}'
    folderpath2 = f'I:\\MERRA2\\Daily_and_Subdaily\\{metpath2}'
    files1 = glob.glob(os.path.join(folderpath1,f'MERRA2.inst3_3d_asm_Np.{i}*.SUB.nc'))
    files2 = glob.glob(os.path.join(folderpath2,f'MERRA2.inst3_3d_asm_Np.{i}*.SUB.nc'))
    difference = len(files1) - len(files2)
            #files = glob.glob(os.path.join(folderpath,"MERRA2.inst3_3d_asm_Np.1*.SUB.nc")) + glob.glob(os.path.join(folderpath,'MERRA2.inst3_3d_asm_Np.20[0-1]*.SUB.nc')) + glob.glob(os.path.join(folderpath,'MERRA2.inst3_3d_asm_Np.202[0-1]*.SUB.nc'))
            #files = glob.glob(os.path.join(folderpath,f'MERRA2.inst3_3d_asm_Np.{i}*.SUB.nc'))
            #files length should be equal to 15341
    print(i,difference)
