# Collect all the scripts to plot Action-levels in one file --URI

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

ntd = 1
D = 5
M = 1
PATH=50
B=30
minima = np.zeros((PATH*B,2))
bad =0
for p in range(PATH):
    data = np.loadtxt('path/D%d_M%d_PATH%d_Ntd%d.dat'%(D,M,p,ntd))
    if data.size==0:
        print "EMPTRY!"
        bad += 1
    else:
        minima[B*p:B*(p+1),:] = data[:,[0,2]]
        print "Failed Paths = ", sum(data[:,1]!=1)


np.savetxt('beta_minima.dat', minima, fmt='%e')

#matplotlib.use("Agg")
font = {'size'   : 18}
matplotlib.rc('font', **font)
minima = np.loadtxt('beta_minima.dat')
fig = plt.figure(figsize=(12,9))
ax1 = fig.add_subplot(111)
ax1.set_title('Lorenz96 SingleF Dim={0} Meas={1} NTau={2}'.format(D,M,ntd))
ax1.set_ylabel('Minimized Action')
ax1.set_xlabel('beta')
ax1.set_yscale('log')

ax1.scatter(minima[:,0],minima[:,1])
plt.savefig('Action_level_D{}_M{}_PATH{}_Ntd{}.png'.format(D,M,p,ntd),bbox_inches='tight')

plt.show()



