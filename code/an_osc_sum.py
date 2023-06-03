import matplotlib.pyplot as plt
import numpy as np

Ax = [0.1, 0.02,   0.01, 0.002,    0.001, ] #    
A = np.array([63, 124, 233 , 236,  224])*1e-18  # dm 0.1



plt.yscale('log')
plt.xscale('log')
plt.plot(Ax, A)
plt.show()

