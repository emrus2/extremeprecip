# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 15:56:16 2023

Calculate temperature advection and q-vectors for defined pressure level

# run script 6 

@author: emrus2
"""
#%% IMPORT MODULES
import os
#import numpy as np
import xarray as xr
#import glob
#import netCDF4 as nc
#import numpy as np
#from datetime import datetime, timedelta
# import metpy
import metpy.calc as mpcalc
from metpy.units import units
# import matplotlib.pyplot as plt

#%% IMPORT MERRA2 DATA

# define pressure level of interest
plev = '850'
plev_int = int(plev) # integer form for q-vector calculation

# hPa temperature
tempfolder =  f'I:\\Emma\\FIROWatersheds\\Data\\DailyMERRA2\\{plev}T'
tempfile = f'MERRA2_{plev}T_Yuba_Extremes90_Daily_1980-2021_WINTERDIST_updated.nc'
temppath = os.path.join(tempfolder,tempfile)
temps = xr.open_dataset(temppath).squeeze()
print(temps)

# hPa winds
windfolder = f'I:\\Emma\\FIROWatersheds\\Data\\DailyMERRA2\\{plev}W'
windfile = f'MERRA2_{plev}W_Yuba_Extremes90_Daily_1980-2021_WINTERDIST_updated.nc'
windpath = os.path.join(windfolder,windfile)
winds = xr.open_dataset(windpath)
winds = winds.squeeze()
print(winds)

#%% COMPUTE TEMPERATURE ADVECTION AND Q-VECTORS
# calculate advection
tempadvection = mpcalc.advection(scalar=temps.T,u=winds.U,v=winds.V,x_dim=-1,y_dim=-2,vertical_dim=None)
print(tempadvection.data.units)
tempadvection = tempadvection.metpy.convert_units('celsius/hour')

# calculate q-vector components
u_qvect, v_qvect = mpcalc.q_vector(u=winds.U,v=winds.V,temperature=temps.T,pressure=plev_int*units.hPa)
print(u_qvect.data.units == v_qvect.data.units)
#print(u_qvect.data.units)
u_qvect = u_qvect.rename('U')
v_qvect = v_qvect.rename('V')
qvects = xr.merge([u_qvect,v_qvect])
#%% SAVE ADVECTION AND Q-VECTOR DATASETS
# temperature advection
os.chdir(f'I:\\Emma\\FIROWatersheds\\Data\\DailyMERRA2\\{plev}TADV')
savefile = f'MERRA2_{plev}TADV_Yuba_Extremes90_Daily_1980-2021_WINTERDIST_updated.nc'
tempadvection.to_netcdf(savefile)

# U q-vector component
os.chdir(f'I:\\Emma\\FIROWatersheds\\Data\\DailyMERRA2\\{plev}QVECT')
savefileQ = f'MERRA2_{plev}QVECT_Yuba_Extremes90_Daily_1980-2021_WINTERDIST_updated.nc'
qvects.to_netcdf(savefileQ)


#%%

# #single variables
# tadv1 = tempadvection[0]
# uwind1, vwind1 =  Uwind[0], Vwind[0]
# u1,v1 = u_qvect[0], v_qvect[0]

# #start figure and set axis
# fig, ax = plt.subplots(figsize=(5, 5))

# # plot isotherms
# cs = ax.contour(temp.lon, temp.lat, temp.values[0], range(274, 300, 2), colors='tab:red',
#                 linestyles='dashed', linewidths=3)
# plt.clabel(cs, fmt='%d', fontsize=16)

# # plot temperature advection and convert units to Kelvin per 3 hours
# cf = ax.contourf(temp.lon, temp.lat, tempadvection[0].metpy.convert_units('kelvin/hour') * 3, range(-6, 7, 1),
#                   cmap=plt.cm.bwr, alpha=0.75)
# plt.colorbar(cf, pad=0, aspect=50)
# ax.barbs(temp.lon.values[::2], temp.lat.values[::2],
#           uwind1[::2, ::2], vwind1[::2, ::2],
#           color='black', length=6, alpha=0.5)
# ax.set(xlim=(260, 270), ylim=(30, 40))
# ax.set_title('Temperature Advection Calculation')

# # plot Q-vectors as arrows, every other arrow
# qvec = ax.quiver(temp.lon.values[::2], temp.lat.values[::2],
#                   u1[::2, ::2] * 1e13, v1[::2, ::2] * 1e13,
#                   color='black', scale=1000, alpha=0.5, width=0.01)

# plt.show()