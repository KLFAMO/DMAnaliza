import matplotlib.pyplot as plt
import numpy as np

x = [-180, -135,-90,-45,0,45,90,135,180]
y50 = np.array([0.2773,0.301, 0.3624, 0.4417, 0.474, 0.4427, 0.3809, 0.305, 0.2773])
y20 = np.array( [0.621, 0.6323, 0.695, 0.871, 0.989, 0.846, 0.772, 0.668, 0.621]  )

plt.plot(x,y50/max(y50), x, y20/max(y20))
plt.show()