import numpy as np
D = 6
M = 2
PATH=100
B=30
minima = np.zeros((PATH*B,2))
for p in range(PATH):
    data = np.loadtxt('path/D%d_M%d_PATH%d_Ntd%d.dat'%(D,M,p,ntd))
    minima[B*p:B*(p+1),:] = data[:,[0,2]]
np.savetxt('beta_minima.dat', minima, fmt='%e')
