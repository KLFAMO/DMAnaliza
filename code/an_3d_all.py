import sys
import os
import pathlib as pa

progspath = pa.Path(__file__).absolute().parents[2]
sys.path.append(str(progspath / 'mytools'))

import matplotlib.pyplot as plt
import tools as tls
import numpy as np

import glob


for D in ['20','50','100','150']:
    for dm in [0.1, 0.25, 0.5, 1, 1.5]: 

        fname = 'D'+D+'_V_-1_-1_-1.npy'
        #fp = os.path.join(progspath,'DMAnaliza','code','out','out03',fname)
        fp = os.path.join(progspath,'DMAnaliza','code','out','out100c2',fname)
        data = np.load(fp)
        m = data[:,0]
        v = np.abs(data[:,1])

        last_mjd = m[0]
        mgap = (float(D)/86400)*0.5

        min_maxv = 1e6
        lenm = len(m)
        i=0
        while i < lenm-1:
            x = i
            start = m[i]
            last_mjd = start
            maxv = v[i]
            state = 1
            #print(m[i], min_maxv)
            while m[x]-start < dm:
                if m[x]-last_mjd > mgap:
                    state = 0
                    break
                if v[x] > maxv:
                    maxv = v[x]
                last_mjd = m[x]
                x = x+1
                if x>=lenm-1:
                    state = 0
                    break

            if min_maxv > maxv and state==1:
                min_maxv = maxv
            
            while m[i]-start < mgap:
                i=i+1
                if i >= lenm-1:
                    break
        print('D ',D,', dm ',dm, ', val ',min_maxv)
    



#di = os.path.join(progspath,'DMAnaliza','code','out','out03')
#flist = glob.glob(di+'/D150*.npy')
#for f in flist:
#    print(f)
#    data = np.load(f)
#    plt.plot(data[:,0],data[:,1])
#plt.show()