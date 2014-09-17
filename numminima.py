# Collect all the scripts to plot Action-levels in one file --URI

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

ntd = 0
D = 5
M = 1
PATH=50
dt=0.01
B=30
thesh = 500
action = np.zeros((B,PATH))
bad =0
for p in range(PATH):
#    data = np.loadtxt('path/D%d_M%d_PATH%d_Ntd%d_dt%d.dat'%(D,M,p,ntd,dt))
    data = np.loadtxt('path/D%d_M%d_PATH%d_Ntd%d_dt%e.dat'%(D,M,p,ntd,dt))
    if data.size==0:
        print "EMPTRY!"
        bad += 1
    else:
        action[:,p] = data[:,2]
        print "Incomplete Paths = ", sum(data[:,1]!=1)

numbelow = np.sum(action<thesh,1)
plt.plot(numbelow)
plt.show()
