"""Plot results from maxvs_*.npy files from an_3D.py"""

import numpy as np
import matplotlib.pyplot as plt

camp = 'c2'
path = '../results/2023_08_no_ptb/'

data = np.load(f'{path}maxvs_{camp}.npy')
plt.plot(data[:,0],data[:,1]*1e-18)
plt.yscale('log')
plt.grid()
plt.savefig(f'{path}maxvs_{camp}.png')