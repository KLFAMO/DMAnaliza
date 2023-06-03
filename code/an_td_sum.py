import matplotlib.pyplot as plt
import numpy as np

Ax = [20, 50, 100, 150] #    
A = np.array([178,   214, 78 , 187 ])*1e-18  # dm 0.1

Bx = [20, 50, 100, 150] #    
B = np.array([879, 364, 305, 314 ])*1e-18  # dm 0.25

Cx = [20, 50, 100, 150] #    
C = np.array([1318, 403, 491, 441 ])*1e-18  # dm 0.5


plt.yscale('log')
plt.plot(Ax, A, Bx, B, Cx, C)
plt.show()

