import sys
import pathlib as pa

progspath = pa.Path(__file__).absolute().parents[2]
sys.path.append(str(progspath / 'mytools'))

import matplotlib.pyplot as plt
import tools as tls
import numpy as np
import dm_tools as dmt

import time

#labs = ['UMK1','UMK2', 'NICT', 'NPLYb', 'NIST', 'NPLSr' ]
labs = ['UMK1', 'UMK2']
tcol = {'UMK1':'green',
        'UMK2':'red',
        'NIST':'blue',
        'NPLSr':'cyan',
        'NPLYb':'black',
        'NICT':'gray'}

d = dict()

start = time.time()

print("Start....")
a = dmt.DMpair(labs[0],labs[1],
                grid_s=10,
                te_s=20,
                fmjd=58665.6 ,tmjd=58665.605,
                grid_mode=0)  # 0 - steps (left val), 1 - interpol
a.plot()

out_tab = []
out_sh = np.arange(-30,30,1)
for sh in out_sh:
    out = a.calc_pair_sh(te_s=20, gshift_s=sh)
    out_tab.append(out)
print(out_sh,out_tab)
print('Time: ', time.time()-start )
plt.plot(out_sh,out_tab)
plt.show()