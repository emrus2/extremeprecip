# -*- coding: utf-8 -*-
"""
Created on Sun Jan 15 13:33:53 2023

@author: emrus2
"""
"""
Masks data outside of watershed for 2022-2023 files
The files saved:
    'UpperYubadailymeans_2022-2023.nc'
"""

#%% IMPORT MODULES
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import os
import shapefile

#%% CREATE FUNCTIONS TO COMPARE SHAPEFILE AND NC FILE
def distance(x1, y1, x2, y2):
    """
    Calculate distance from (x1,y1) to (x2,y2)
    """
    return ((x1-x2)**2 + (y1-y2)**2)**0.5

def point_is_on_line(x, y, x1, y1, x2, y2):
    """
    Check whether point (x,y) is on line (x1,y1) to (x2,y2)
    """

    d1 = distance(x,  y,  x1, y1)
    d2 = distance(x,  y,  x2, y2)
    d3 = distance(x1, y1, x2, y2)

    eps = 1e-12
    return np.abs((d1+d2)-d3) < eps

def is_left(xp, yp, x0, y0, x1, y1):
    """
    Check whether point (xp,yp) is left of line segment ((x0,y0) to (x1,y1))
    returns:  >0 if left of line, 0 if on line, <0 if right of line
    """

    return (x1-x0) * (yp-y0) - (xp-x0) * (y1-y0)

def is_inside(xp, yp, x_set, y_set, size):
    """
    Given location (xp,yp) and set of line segments (x_set, y_set), determine
    whether (xp,yp) is inside polygon.
    """

    # First simple check on bounds
    if (xp < x_set.min() or xp > x_set.max() or yp < y_set.min() or yp > y_set.max()):
        return False

    wn = 0
    for i in range(size-1):

        # Second check: see if point exactly on line segment:
        if point_is_on_line(xp, yp, x_set[i], y_set[i], x_set[i+1], y_set[i+1]):
            return False

        # Calculate winding number
        if (y_set[i] <= yp):
            if (y_set[i+1] > yp):
                if (is_left(xp, yp, x_set[i], y_set[i], x_set[i+1], y_set[i+1]) > 0):
                    wn += 1
        else:
            if (y_set[i+1] <= yp):
                if (is_left(xp, yp, x_set[i], y_set[i], x_set[i+1], y_set[i+1]) < 0):
                    wn -= 1

    if wn == 0:
        return False
    else:
        return True

def calc_mask(mask, lon, lat, shp_lon, shp_lat):
    """
    Calculate mask of grid points which are inside `shp_lon, shp_lat`
    """
    for j in range(lat.size):    
        for i in range(lon.size):
            if is_inside(lon[i], lat[j], shp_lon, shp_lat, shp_lon.size):
                mask[:,j,i] = True
#%%
# IMPORT WATERSHED DATA
watershed_name= 'UpperYuba'
ws_directory = f'I:\\Emma\\FIROWatersheds\\Data\\WatershedShapefiles\\California\\{watershed_name}\\'
ws_filename = (f'{watershed_name}.shp')
ws_path = os.path.join(ws_directory,ws_filename)
shp = shapefile.Reader(ws_path, crs="epsg:4326")
feature = shp.shapes()


#IMPORT NETCDF DATA
year = "2024"
nc_filename = (f'pr_{year}.nc')
nc_directory = 'I:\\GRIDMET\\pr\\'
nc_path = os.path.join(nc_directory,nc_filename)
ncf = nc.Dataset(nc_path,mode='r')
nc_lon = ncf.variables['lon'][:]
nc_lat = ncf.variables['lat'][:]
nc_pr  = ncf.variables['precipitation_amount'][:]
ncf.close()

# Extract array of lat/lon coordinates:
coords = np.squeeze(np.array([s.points for s in shp.shapes()]))
shp_lon = np.array(coords)[:,0]
shp_lat = np.array(coords)[:,1]


# Calculate mask
mask = np.zeros_like(nc_pr, dtype=bool)
calc_mask(mask, nc_lon, nc_lat, shp_lon, shp_lat)

# Mask the data array
nc_pr_masked = np.ma.masked_where(~mask, nc_pr)
np.ma.asanyarray(nc_pr_masked)
# Calculate spatial average
daily_means = np.ma.masked_array.mean(nc_pr_masked,axis=(1,2)) #should be 365 (1d)
daily_means_thresh1 = np.ma.where(daily_means > 2,daily_means,'')
finalmeans = list(daily_means_thresh1)
#%%
# save data
savefolder = 'I:\\Emma\\FIROWatersheds\\Data\\GRIDMET'
np.save(os.path.join(savefolder,f'{watershed_name}dailymeanprecip_{year}'),finalmeans)

#%%
for i in range(365):
# Plot!
    plt.figure(figsize=(8,4))
    plt.title(i)
    plt.pcolormesh(nc_lon, nc_lat, nc_pr_masked[i,:,:])
    plt.xlim(-122, -120)
    plt.ylim(39, 40)
    plt.colorbar()
    plt.tight_layout()
    plt.show()