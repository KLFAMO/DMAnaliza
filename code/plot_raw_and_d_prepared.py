import numpy as np
import matplotlib.pyplot as plt
import pathlib as pa
progspath = pa.Path(__file__).absolute().parents[2]

print(progspath)
#d = np.load(str(progspath)+'/DMAnaliza/data/d_prepared/d_UMK2_c1.npy', allow_pickle=True)
r = np.load(str(progspath)+'/DMAnaliza/data/raw/raw_UMK2_c2.npy', allow_pickle=True)
#r = r[r[:,0].argsort()]
#print(r)
r[:,1] = r[:,1]*0.001
plt.plot(r[:,0],r[:,1])
#plt.plot(d[:,0],d[:,1])
plt.show()
np.save(str(progspath)+'/DMAnaliza/data/d_prepared/d_UMK2_c2.npy', r)