
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

ntd = 1
D = 5
M = 1
PATH=100
dt=0.01
thresh = 13000
lowerbnd = 11500
B=29
N = 300


bad = 0
Ydata = np.loadtxt('dataN_dt01_noP.txt')

for p in range(PATH):
#    data = np.loadtxt('path/D%d_M%d_PATH%d_Ntd%d_dt%d.dat'%(D,M,p,ntd,dt))
    data = np.loadtxt('path/D%d_M%d_PATH%d_Ntd%d_dt%e.dat'%(D,M,p,ntd,dt))
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
