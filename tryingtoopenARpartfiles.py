# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 13:00:50 2024

@author: emrus2
"""

import numpy as np
import netCDF4 as nc

filepath = 'I:/Emma/ARCatalogData/globalARcatalog_MERRA2_1980-2021_v3.0.nc'
gridfile = nc.Dataset(filepath,mode='r')
