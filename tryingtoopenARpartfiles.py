# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 13:00:50 2024

@author: emrus2
"""

import numpy as np
import netCDF4 as nc

filepath = 'C:/Users/emrus2/Downloads/globalARcatalog_MERRA2_1980-2021_2_v3.0.nc'
gridfile = nc.Dataset(filepath,mode='r')
