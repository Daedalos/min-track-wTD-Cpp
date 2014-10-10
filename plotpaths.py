
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

ntd = 0
D = 5
M = 1
PATH=100
dt= 0.1
thresh = 11000
lowerbnd = 0
B=29
N = 300


bad = 0
Ydata = np.loadtxt('dataN_dt01_noP.txt')

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
        if data[B, 2] < thresh and data[B,2]>lowerbnd:
            path = data[B,3:]
            path = np.reshape(path, (-1,D))
            if data.shape[1] != D*N+3:
                import ipdb; ipdb.set_trace()
                raise ValueError("Wrong Dims")

            plt.figure()
            plt.title("Path {0}, Action={1}, Beta={2}".format(p,data[B,2],B))
            plt.plot(path[:,1],'r',label='Estimate')
            plt.plot(Ydata[:N,1], 'b',label='Data')
        print "Incomplete Paths = ", sum(data[:,1]!=1)

plt.show()
