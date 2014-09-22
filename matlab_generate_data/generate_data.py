from numpy import *
from scipy.integrate import odeint


def lorenz96(x, t):
    N = len(x)
    dxdt = zeros(N)
    for i in range(N):
        dxdt[i] = x[mod(i-1,N)]*(x[mod(i+1,N)]-x[mod(i-2,N)]) - x[i] + 8.17

    return dxdt

dt = 0.01

T_total = arange(0,100,dt)
#data_initial = randn(5+1,1)
#data_initial = [0.80, 0.95, 0.71, 0.24, 0.63,-1, 2, 1, -2.2, 0.3];
#data_initial[-1] = 8.17
data_initial = [0.80, 0.95, 0.71, 0.24, 0.63];

Y = odeint(lorenz96,data_initial,T_total)
Y = Y[1000:2501,:]
param = 8.17*ones(len(Y))
Y = column_stack((Y,param))

savetxt("data_D{}_dt{}.txt".format(len(data_initial),dt),Y)
noise = 0.5*random.standard_normal(Y.shape)
Ynoise = Y + noise
savetxt("dataN_D{}_dt{}.txt".format(len(data_initial),dt),Ynoise)

