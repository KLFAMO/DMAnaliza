
import numpy as np
import matplotlib.pyplot as plt

d1 = np.load('maxvs_c1.npy')
d2 = np.load('maxvs_c2.npy')
y = [min(x1,x2)*1e-18 for x1, x2 in zip(d1[:,1], d2[:,1])]
plt.plot(d1[:,0],y)
plt.yscale('log')
plt.grid()
plt.savefig('maxvs_c1_c2.png')