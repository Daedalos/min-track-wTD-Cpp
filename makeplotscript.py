# Collect all the scripts to plot Action-levels in one file --URI

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

taus = []
ntd = len(taus)
D = 10
M = 1
PATH=100
dt=0.01
B=30
minima = np.zeros((PATH*B,2))
bad =0

filetemp = 'pbD{0}'
for i in range(len(taus)):
    filetemp += '-{}'.format(taus[i])
if ntd ==0:
    filetemp += '-0'

filetemp += '/path/D{0}_M{1}_PATH{2}_Ntd{3}_dt{4:e}.dat'

#filetemp = 'path/D%d_M%d_PATH%d_Ntd%d_'
#for i in range(len(taus)):
#    filetemp += '%d-' % taus[i]
#filetemp += 'dt%e.dat'

print filetemp
points = 0

for p in range(PATH):
    #    data = np.loadtxt('path/D%d_M%d_PATH%d_Ntd%d_dt%d.dat'%(D,M,p,ntd,dt)) 
    try:
        filename = filetemp.format(D,M,p,ntd,dt)
        data = np.loadtxt(filename)
        #data = np.loadtxt(filetemp%(D,M,p,ntd,dt))
    except: 
      
        continue

    if data.size==0:
        print "EMPTRY!"
        bad += 1
    else:
        if np.isnan(data[-1,0]): continue
        points +=1
        minima[B*p:B*(p+1),:] = data[:,[0,2]]
        print "Failed Paths = ", sum(data[:,1]!=1)


print 'points plotted = ', points
np.savetxt('beta_minima.dat', minima, fmt='%e')

#matplotlib.use("Agg")
font = {'size'   : 18}
matplotlib.rc('font', **font)
minima = np.loadtxt('beta_minima.dat')
fig = plt.figure(figsize=(12,9))
ax1 = fig.add_subplot(111)
ax1.set_title('Lorenz96 SingleF Dim={0} Meas={1} Tau={2}'.format(D,M,taus))
ax1.set_ylabel('Minimized Action')
ax1.set_xlabel('beta')
ax1.set_yscale('log')

#import ipdb; ipdb.set_trace()

ax1.scatter(minima[:,0],minima[:,1])
plt.savefig('Action_level_D{}_M{}_PATH{}_Ntd{}.png'.format(D,M,p,ntd),bbox_inches='tight')

plt.show()



