import scipy.optimize as scp
import numpy as np
import matplotlib.pyplot as plt

sd = [
        [1, 20],
        [1, 20],
        [1, 0],
        [1, 0],
        [1, 0],
        [1, 0],
        [1, 0]
]

etau = 10000  # s
etaum = etau/86400

def fu(dx,A,sh):
    return [f(x,A,sh) for x in dx]

def f(x,A,sh):
    rx = round(x)
    K = sd[rx][0]
    if x-rx<2*etaum:
        return sh
    elif x-rx<3*etaum:
        return A*K+sh
    else:
        return sh

def s(dx):
    return [ sd[round(x)][1] for x in dx]

x = np.arange(0,1.99,0.001)
sig = x*0+0.5
sig[0:1000] = 0.1
ydat = fu(x,1,0) + np.random.rand(len(x))-0.5
ydat[0:1000] = ydat[0:1000]*2

popt, pcov = scp.curve_fit(fu, x, ydat, sigma=sig, absolute_sigma=True )
print('POPT:  ', popt)


plt.plot(x, ydat, x, fu(x,*popt), x, sig)
plt.show()
