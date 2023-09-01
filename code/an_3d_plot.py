import numpy as np
import matplotlib.pyplot as plt

camp = 'c2'

data = np.load('maxvs_'+camp+'.npy')
plt.plot(data[:,0],data[:,1]*1e-18)
plt.yscale('log')
plt.grid()
plt.savefig('maxvs_'+camp+'.png')