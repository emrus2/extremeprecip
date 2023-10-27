# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 14:18:57 2023

@author: emrus2
"""

import numpy as np
import matplotlib.pyplot as plt
import os
#%%
daysprior = 4
clusters = daysprior + 1

mat_dir='I:\\Emma\\FIROWatersheds\\Data\\SOMs\\SomOutput'
os.chdir(mat_dir)
sammon = np.load(f'Sammon_{clusters}d_test.npy')

#%%
fig, ax = plt.subplots(layout='constrained')

Y = -sammon[:,1]
X = sammon[:,0]
ax.plot(X,Y,'ko')
lw = 1
col = 'k'

leftcorners = (-800,-450)
rightcorners = (450,-450)
centers = ()
for i in range(len(sammon)):
    xannot = {0:X[i]+leftcorners[0],1:X[i]-50,2:X[i]+rightcorners[0],3:X[i]-950,4:X[i]+200,5:X[i]+450,6:X[i]+leftcorners[0],7:X[i]+250,8:X[i]+rightcorners[0]}
    yannot = {0:Y[i]+leftcorners[1],1:Y[i]-1700,2:Y[i]+rightcorners[1],3:Y[i]-250,4:Y[i]-1500,5:Y[i]-450,6:Y[i]+leftcorners[1],7:Y[i]+750,8:Y[i]+rightcorners[1]}
    ax.annotate(i+1,(xannot[i],yannot[i]))
    if i not in range(2,9,3):
        ax.plot(X[i:i+2],Y[i:i+2],color=col,linewidth=lw)
    if i < 6:
        x1,y1 = X[i],Y[i]
        x2,y2 = X[i+3],Y[i+3]
        Xlin = [x1,x2]
        Ylin = [y1,y2]
        ax.plot(Xlin,Ylin,color=col,linewidth=lw)
ax.set_ybound(-17500,14000)
ax.set_xticks([])
ax.set_yticks([])
for spine in ax.spines:
    ax.spines[spine].set_visible(False)

#%%
save_dir='I:\\Emma\\FIROWatersheds\\Figures\\SOMs'
os.chdir(save_dir)
plt.savefig(f'Sammon_{clusters}d.png',dpi=300)
plt.show()
