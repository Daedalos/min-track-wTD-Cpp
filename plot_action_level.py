import numpy as np
import matplotlib
import matplotlib.pyplot as plt
#matplotlib.use("Agg")
font = {'size'   : 18}
matplotlib.rc('font', **font)
minima = np.loadtxt('beta_minima.dat')
fig = plt.figure(figsize=(12,9))
ax1 = fig.add_subplot(111)
ax1.scatter(minima[:,0],minima[:,1])
ax1.set_title('Lorenz96 SingleF Dim=5 #Meas=2')
ax1.set_ylabel('Minimized Action')
ax1.set_xlabel('beta')
ax1.set_yscale('log')
plt.savefig('Action_level.png',bbox_inches='tight')
