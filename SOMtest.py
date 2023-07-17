# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 13:49:49 2023

@author: emrus2

First draft of SOM analysis using MINISOM toolbox
Not used in analysis - used MATLAB SOM TOOLBOX instead

UPDATED 6/12/2023
"""
#%% IMPORT MODULES
import netCDF4 as nc #if this doesn't work, try to reinstall in anaconda prompt using
    #conda install netcdf4
import numpy as np
import os
import matplotlib.pyplot as plt
#from matplotlib.colors import ListedColormap
#import matplotlib.cm as cm
os.environ["PROJ_LIB"] = os.path.join(os.environ["CONDA_PREFIX"], "share", "proj")
from mpl_toolkits.basemap import Basemap #installed using 
    #conda install -c anaconda basemap
from minisom import MiniSom
import scipy.io

#define watershed and directory for plotting
watershed = 'UpperYuba'
ws_directory = f'I:\\Emma\\FIROWatersheds\\Data\\WatershedShapefiles\\California\\{watershed}\\'

#%% IMPORT MERRA2
#define MERRA2 data location
#metvars = ['Z500', 'SLP', '850T', '300W', 'IVT']
metvar = 'IVT'
folderpath = f'I:\\Emma\\FIROWatersheds\\Data\\DailyMERRA2\\{metvar}'
filename = f'MERRA2_{metvar}_Yuba_Extremes_Daily_1980-2021_WINTERDIST.nc'
filepath = os.path.join(folderpath,filename)

#COLLECT VARIABLE DATA FROM MERRA2 FILE
#merravar = {'Z500':'H','SLP':'SLP','850T':'T'}
#open the netcdf file in read mode
gridfile = nc.Dataset(filepath,mode='r')
print(gridfile)
gridlat = gridfile.variables['lat'][:]
gridlon = gridfile.variables['lon'][:]
Uvapor = gridfile.variables['UFLXQV'][:]
Vvapor = gridfile.variables['VFLXQV'][:]
merra = np.squeeze(np.sqrt(Uvapor**2 + Vvapor**2))
gridfile.close() 

#%% REDUCE TRAINING DATA TO REGION OF INTEREST
#reduce data to desired restrictions
#define lat lon restrictions
latmin, latmax = (20.5,70.5)
lonmin, lonmax = (-170.25,-105.75)
#reduce lat
latlims = np.logical_and(gridlat > latmin, gridlat < latmax)
latind = np.where(latlims)[0]
gridlatreduced = gridlat[latind]
#reduce lon
lonlims = np.logical_and(gridlon > lonmin, gridlon < lonmax)
lonind = np.where(lonlims)[0]
gridlonreduced = gridlon[lonind]
#reduce pressure
merrareduced = merra[:,latind,:]
merrareduced = merrareduced[:,:,lonind]

#%% PREPARE TRAINING DATA
#define data as numpy array (for minisom package)
data = np.asarray(merrareduced) #(130x99x103)
data_cp = np.copy(data)

#normalize by temporal standard deviation
#can wait to do this step, seems only important for multivariate SOMs
''' Loikith 2017 uses temporal std, while Aragon 2020
and Taylor 2023 use the spatial std to normalize '''
# for index, day in enumerate(data):
#     print(np.mean(day),np.std(day))
#     day -= np.mean(day)
#     day /= np.std(day)
        
#weight data based on area (square root of the cosine of latitude)
latweights = np.sqrt(np.cos((np.radians(gridlatreduced))))
for i in range(data.shape[1]): #lats
    data[:,i,:] *= latweights[i]

#confirm dataset has been correctly weighted
datatest = data_cp - data
print(np.all(datatest) == 0) #Should be FALSE
    
#flatten data from 3D to 2D (130x10197)
data = data.reshape(len(data),-1)

#%% SAVE DATA ARRAY TO OPEN IN MATLAB
mat_dir='I:\\Emma\\FIROWatersheds\\Scripts\\Matlab'
os.chdir(mat_dir)
datadict = {'lat': gridlatreduced, 'lon': gridlonreduced, 'IVT': data}
scipy.io.savemat('data.mat', mdict = datadict)
#%% TRAIN SOM
# som = MiniSom(x=4, y=3, input_len=data.shape[1], sigma=1.0, learning_rate=0.5, neighborhood_function='gaussian')
#     #sigma is the neighborhood radius
#     # initialization of 3x4 neuron SOM
#     # want to have a small learning rate
#     # want large number of initializations

# # radius sizes for initial and final training
# rad_init = 4
# rad_fin = 1
# # training iterations for initial and final training
# iter_init = 100
# iter_fin = 2000

# # train SOM with batch version
# som.train_batch(data=data, num_iteration=100, verbose=True)
# # print(som.winner(data[0]))

# # define quantization and topographic errors of data
# quanterr = som.quantization_error(data)
# topoerr = som.topographic_error(data)

 
# # find the best-matching unit (BMU) time series
# # in this module, this is the 'winning' SOM, where we can then find the weights
# #for item in enumerate(data):
# #    print(item) #item is tuple: with integer and array
#     #print(som.winner(data[item]))
# #som.winner(data[0])
#     #prints out the winning (i,j) node for the day of desired index
# mapped_data = som.win_map(data)
#     #seems to produce a dictionary of arrays and the associated winning SOM
# #weightings = som.get_weights()
#     #gets weights in a 4x3x10197 array
# #print(mapped_data)

# ''' mapped data seems to provide us with all the information needed, all the arrays
# that are assigned to that node are clustered into a dictionary.
# The days of these arrays are not present, but the number of days assigned
# to each node are available. However, how do we know what the node pattern
# looks like? '''
#%% PLOT NODES
# for i, key in enumerate(mapped_data):
#     assigned = mapped_data[key]
#     nodemean = np.mean(assigned, axis=0)
#     merrareduced = nodemean.reshape(len(gridlatreduced),len(gridlonreduced))

#     fig = plt.figure(figsize=(5.7,4))
#     colormap = 'gnuplot2_r'
#     cbarlabs = 'kg $\mathregular{m^{-1}}$ $\mathregular{s^{-1}}$'
#     plottitle = f'{key}'
    
#     #MAP DESIRED VARIABLE
#     #convert lat and lon into a 2D array
#     lon, lat = np.meshgrid(gridlonreduced,gridlatreduced) 
#     #define area threshold for basemap
#     area_thresh = 1E4
#     #create equidistant cylindrical projection basemap
#     map = Basemap(projection='cyl',llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,\
#               urcrnrlon=lonmax,resolution='l',area_thresh=area_thresh)
#     xi, yi = map(lon,lat)
#     #create colormap of MERRA2 data
#     colorm = map.pcolor(xi,yi,merrareduced,shading='auto',cmap=colormap, zorder=1)
    
#     #define border color and thickness
#     border_c = '0.4'
#     border_w = 0.4
#     #create map features
#     map.drawcoastlines(color=border_c, linewidth=border_w)
#     map.drawstates(color=border_c, linewidth=border_w)
#     map.drawcountries(color=border_c, linewidth=border_w)
#     gridlinefont = 9
#     map.drawparallels(np.arange(30.,71.,15.), labels=[1,0,0,0], fontsize=gridlinefont,color=border_c, linewidth=border_w)
#     map.drawmeridians(np.arange(-160.,-109.,20.), labels=[0,0,0,1], fontsize=gridlinefont,color=border_c, linewidth=border_w)
#     #define contour color and thickness
#     contour_c = '0.1'
#     contour_w = 0.7
#     #create contour map
#     contourm = map.contour(xi,yi,merrareduced,colors=contour_c,linewidths=contour_w)
#     plt.clabel(contourm,fontsize=6,inline_spacing=1,colors='k',zorder=5,manual=False)
#     plt.title(plottitle,fontweight="bold",fontsize=12)
        
#     #add yuba shape
#     #map.readshapefile(os.path.join(ws_directory,f'{watershed}'), watershed,linewidth=0.8,color='r')
#     plt.scatter(-120.9,39.5,color='tomato',edgecolors='r',marker='*',linewidths=0.8,zorder=6)
    
#     #cbar_ax = fig.add_axes([0.15,0.04,0.7,0.0216]) #bottom colorbar
#     cbar = map.colorbar(colorm, location='right', pad="6%")
#     #cbar = map.colorbar(colorm, location='bottom', pad="10%",ticks=np.arange(cbarstart[metvar],highlims[metvar]+1,cbarint[metvar]))
#     cbar.ax.tick_params(labelsize=9)
#     cbar.set_label(cbarlabs,fontsize=9)
    
#     #SHOW MAP
#     save_dir='I:\\Emma\\FIROWatersheds\\Figures\\SOMs'
#     os.chdir(save_dir)
#     plt.savefig(f'{metvar}_SOM_{key}.png',dpi=300)
#     plt.show()
    
#%% PLOT NODES from MATLAB
save_dir='I:\\Emma\\FIROWatersheds\\Scripts\\Matlab'
os.chdir(save_dir)
soms = scipy.io.loadmat('sompatterns.mat')
patterns = soms['pats']

zmax = 0
zmin = 1E8

for i, arr in enumerate(patterns):
    merrareduced = arr.reshape(len(gridlatreduced),len(gridlonreduced))
    
    #determine zmax and zmin for all days
    highlim = np.amax(merrareduced)
    lowlim = np.amin(merrareduced)
    print(lowlim,highlim)
    if highlim > zmax:
        zmax = highlim
    if lowlim < zmin:
        zmin = lowlim

print(f'Lowest Value:{zmin} \nHighest Value:{zmax}')

    fig = plt.figure(figsize=(5.7,4))
    colormap = 'gnuplot2_r'
    cbarlabs = 'kg $\mathregular{m^{-1}}$ $\mathregular{s^{-1}}$'
    plottitle = i + 1
    
    #MAP DESIRED VARIABLE
    #convert lat and lon into a 2D array
    lon, lat = np.meshgrid(gridlonreduced,gridlatreduced) 
    #define area threshold for basemap
    area_thresh = 1E4
    #create equidistant cylindrical projection basemap
    map = Basemap(projection='cyl',llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,\
              urcrnrlon=lonmax,resolution='l',area_thresh=area_thresh)
    xi, yi = map(lon,lat)
    #create colormap of MERRA2 data
    colorm = map.pcolor(xi,yi,merrareduced,shading='auto',cmap=colormap, zorder=1)
    
    #define border color and thickness
    border_c = '0.4'
    border_w = 0.4
    #create map features
    map.drawcoastlines(color=border_c, linewidth=border_w)
    map.drawstates(color=border_c, linewidth=border_w)
    map.drawcountries(color=border_c, linewidth=border_w)
    gridlinefont = 9
    map.drawparallels(np.arange(30.,71.,15.), labels=[1,0,0,0], fontsize=gridlinefont,color=border_c, linewidth=border_w)
    map.drawmeridians(np.arange(-160.,-109.,20.), labels=[0,0,0,1], fontsize=gridlinefont,color=border_c, linewidth=border_w)
    #define contour color and thickness
    contour_c = '0.1'
    contour_w = 0.7
    #create contour map
    contourm = map.contour(xi,yi,merrareduced,colors=contour_c,linewidths=contour_w)
    plt.clabel(contourm,fontsize=6,inline_spacing=1,colors='k',zorder=5,manual=False)
    plt.title(plottitle,fontweight="bold",fontsize=12)
        
    #add yuba shape
    #map.readshapefile(os.path.join(ws_directory,f'{watershed}'), watershed,linewidth=0.8,color='r')
    plt.scatter(-120.9,39.5,color='tomato',edgecolors='r',marker='*',linewidths=0.8,zorder=6)
    
    #cbar_ax = fig.add_axes([0.15,0.04,0.7,0.0216]) #bottom colorbar
    cbar = map.colorbar(colorm, location='right', pad="6%")
    #cbar = map.colorbar(colorm, location='bottom', pad="10%",ticks=np.arange(cbarstart[metvar],highlims[metvar]+1,cbarint[metvar]))
    cbar.ax.tick_params(labelsize=9)
    cbar.set_label(cbarlabs,fontsize=9)
    
    #SHOW MAP
    save_dir='I:\\Emma\\FIROWatersheds\\Figures\\SOMs'
    os.chdir(save_dir)
    plt.savefig(f'{metvar}_SOM_{i + 1}.png',dpi=300)
    plt.show()
#%%
   
# [bmus, qerrs] = som_bmus(sMap, sD); #% bmus is the SOM pattern time 
# # % series, and qerrs is the quantization error time series.
# # % Place the dates, best-matching pattern, and quantization errors 
# # % in one time series array
# clear sD
# timeseries = NaN(TotDays,4); % First two columns are year and day, 
# # % last two columns are best-matching pattern and rms error
# # % First add the best-matching pattern and quantization/rms error 
# # % values
# timeseries(:,3) = bmus;
# timeseries(:,4) = qerrs;
       
# # % Determine the frequency of occurrence of each pattern
# K = num_rows*num_cols; % Number of SOM patterns
        
# pat_freq = zeros(1,K);
# for p = 1:K
       
#     ind = find(timeseries(:,3) == p);
#     pat_freq(p) = length(ind)/num_obs; 
            
# end

# # % Get the SOM pattern matrix from the codebook
# pats = sMap.codebook;

# # Plotting
# # %identify the SOM assignments
# som=timeseries(:,3);
# #save your data

#%% plot training results
# plt.figure(figsize=(9, 9))
# plt.pcolor(som.distance_map().T, cmap='bone_r')  # plotting the distance map as background
# plt.colorbar()
# # Plotting the response for each pattern in the iris dataset
# for cnt, xx in enumerate(data):
#     w = som.winner(xx)  # getting the winner
#     plt.plot(w[0]+.5, w[1]+.5)
# plt.show()