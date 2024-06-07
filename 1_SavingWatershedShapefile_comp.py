# -*- coding: utf-8 -*-
"""
Created on Sun Jan 15 13:19:14 2023

@author: emrus2

Opens watershed dataset and restricts to desired single watershed shapefile

Saved File: 
    'I:\\Emma\\FIROWatersheds\\Data\\WatershedShapefiles\\California\\UpperYuba\\
        UpperYuba.shp''
"""

#IMPORT MODULES
import os
import geopandas as gpd


#IMPORT CA SHAPEFILE
#define CA Watersheds shapefile location
ws_filename = ('WBDHU8.shp')
ws_directory = 'I:\\Emma\\FIROWatersheds\\Data\\WatershedShapefiles\\California\\'
ws_path = os.path.join(ws_directory,ws_filename)
ws = gpd.read_file(ws_path, crs="epsg:4326")

#RESTRICT TO DESIRED WATERSHED
huc = "18020125"
#restrict dataframe to Upper Yuba watershed through HUC8 code
upperyuba = ws.loc[ws['huc8'] == huc]
upperyuba.to_file(os.path.join(ws_directory + 'UpperYuba','UpperYuba.shp'), driver='ESRI Shapefile')