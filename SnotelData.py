# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 14:27:33 2024

@author: emrus2
"""
import numpy as np
import os
import pandas as pd
from datetime import datetime

# import extreme precipitation days
extremedays = np.load('I://Emma//FIROWatersheds//Data//90Percentile_ExtremeDays.npy')
extremedays = [str(i) for i in extremedays]

# open csv file
savefolder = 'I:/Emma/Classes/SnowHydrology/'
df = pd.read_csv(os.path.join(savefolder,'SNOTEL_upperyuba.csv'))
print(df)
# remove excess information from dataset
dfred = df.iloc[58:,:-1]
# convert to numpy array
datanp = dfred.to_numpy()

# create empty list to store change in swe
deltaSWE = []
# convert datetimes to match format of extreme days
for i,arr in enumerate(datanp):
    # print(i,arr)
    datedt = datetime.strptime(arr[0],'%Y-%m-%d')
    # datanp[i,0] = datetime.strftime(datedt,'%Y%m%d')
    datestr = datetime.strftime(datedt,'%Y%m%d')
    # if an extreme day
    if datestr in extremedays:
        # calculate change in SWE across extreme day
        print(datanp[i-1,1],datestr,datanp[i+1,1])
        # print(datanp[i+1,0])
        SWE_prev = float(datanp[i-1,1])
        SWE_aft = float(datanp[i+1,1])
        SWE_diff = SWE_aft - SWE_prev
        deltaSWE.append(SWE_diff)

# save deltaSWE data
np.savetxt(os.path.join(savefolder,'deltaSWE.csv'), \
           deltaSWE,delimiter=",",fmt='%s')

