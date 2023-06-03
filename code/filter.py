import sys
import os
import pathlib as pa
import numpy as np
progspath = pa.Path(__file__).absolute().parents[2]
sys.path.append(str(progspath / 'mytools'))
import tools as tls

T = 200
f=1/T
x=np.arange(0,60*60, 1)
y=np.sin(2*np.pi*f*x)
a = tls.TSerie(mjd=x, val=y)
a.high_gauss_filter(stddev=350)
print(a)
a.plot()
