# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 14:18:57 2023

@author: emrus2
"""

import numpy as np
import matplotlib.pyplot as plt
import os
#%%
daysprior = 2
clusters = daysprior + 1

mat_dir='I:\\Emma\\FIROWatersheds\\Data\\SOMs\\SomOutput'
os.chdir(mat_dir)
sammon = np.load(f'Sammon_{clusters}d.npy')

#%%
fig, ax = plt.subplots(layout='constrained')

Y = -sammon[:,1]
X = sammon[:,0]
ax.plot(X,Y,'ko')
lw = 1
col = 'k'

for i in range(9):
    ax.annotate(i+1,(X[i]-400,Y[i]-1700))
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
