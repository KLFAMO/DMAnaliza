import numpy as np
import matplotlib.pyplot as plt

d = np.load("out/out50abc_D13.npy")
mjd = d[:,0]
val = d[:,1]
plt.plot(mjd,val)
plt.show()