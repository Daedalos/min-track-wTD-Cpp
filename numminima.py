# Count number of paths with Action below threshold

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

ntd = 3
taus = [10,20,30]
D = 5
M = 1
PATH=100
dt=0.01
Bstart = 20
B=30
thresh = 1
lowbnd = 0
action = np.zeros((B-Bstart,PATH))*np.nan
bad =0

filetemp = 'path/D%d_M%d_PATH%d_Ntd%d_'
for i in range(len(taus)):
    filetemp += '%d-' % taus[i]
filetemp += 'dt%e.dat'

for p in range(PATH):
#    data = np.loadtxt('path/D%d_M%d_PATH%d_Ntd%d_dt%d.dat'%(D,M,p,ntd,dt))
    try:
        data = np.loadtxt(filetemp%(D,M,p,ntd,dt))
    except: continue

    if data.size==0:
        print "EMPTRY!"
        bad += 1
    else:
        data = data[Bstart:B,:]
        action[:,p] = data[:,2]
        print "Incomplete Paths = ", sum(data[:,1]!=1)


action[np.isnan(action)] = -1

between = np.logical_and(action<thresh, action>lowbnd)

numbelow = np.sum( between ,1)
plt.plot(range(Bstart,B), numbelow)
plt.show()
