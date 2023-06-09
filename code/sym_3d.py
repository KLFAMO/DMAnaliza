import sys
import os
import pathlib as pa
progspath = pa.Path(__file__).absolute().parents[2]
sys.path.append(str(progspath / 'mytools'))
import matplotlib.pyplot as plt
import tools as tls
import numpy as np
import dm_tools as dmt
import scipy.optimize as scp
import time
import decimal as dec

etaum = 0

#labs = ['UMK1','UMK2', 'NIST', 'NPLSr', 'NPLYb', 'NICT', 'SYRTE']
labs = ['UMK1','UMK2', 'NIST', 'SYRTE', 'NPLSr', 'NPLYb', 'NICT']
#labs = ['NIST', 'SYRTE', 'NPLSr', 'NPLYb', 'NICT']
#labs = ['NIST', 'SYRTE', 'NPLSr', 'NPLYb', 'NICT']
#labs = ['NPLSr', 'NPLYb', 'NIST']
inf = { 'UMK1': {'col':'green', 'atom':'88Sr',
                 'X':3644273,  'Y':1226649,  'Z':5071736}, 
        'UMK2': {'col':'red',   'atom':'88Sr',
                 'X':3644273,  'Y':1226649,  'Z':5071736},
        'NIST': {'col':'blue',  'atom':'171Yb',
                 'X':-1288363, 'Y':-4721684, 'Z':4078659},
        'NPLSr':{'col':'cyan',  'atom':'87Sr',
                 'X':3985500,  'Y':-23625,   'Z':4962941},
        'NPLYb':{'col':'black', 'atom':'171Yb+',
                 'X':3985500,  'Y':-23625,   'Z':4962941},
        'NICT': {'col':'gray',  'atom':'87Sr',
                 'X':-3941931, 'Y':3368182,  'Z':3702068},
        'SYRTE':{'col':'brown', 'atom':'87Sr',
                 'X':4202777,  'Y':171368,   'Z':4778660}
}

lnum = {'UMK1':0, 'UMK2':1, 'NIST':2, 'NPLSr':3, 'NPLYb':4,
        'NICT':5, 'SYRTE':6}

# 0: K, 1: std
sd = [0,0,0,0,0,0,0]

def fu(dx,A,sh):
    return [f(x,A,sh) for x in dx]

def f(x,A,sh):
    ts = 8640       #1/ (10s / 86400 )
    global etaum
    rx = int(x)
    if x-rx<2*etaum:
        return sh
    elif x-rx<3*etaum:
        return A*(1-np.exp(-(x-rx-2*etaum)*ts))+sh   
    else:
        return ( A*(1-np.exp(-etaum*ts)) * np.exp(-(x-rx-3*etaum)*ts) ) + sh

def sigf(dx):
    return [ sd[int(x)] for x in dx]

def vmul(a,b):
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]

def calc_single(mjd, v, D, vec):
    """
    v - speed
    D - size
    vec - direction
    """
    global etaum
    global sd

    vec_abs = np.sqrt(vec[0]**2+vec[1]**2+vec[2]**2)
    if vec_abs!=0:
        vec = vec/vec_abs
    etau = D/v  # s
    etaum = etau/86400
    end_mjd = mjd+2*etaum + 7/8640
    durm = end_mjd-mjd
    datx=[]
    daty=[]

    cnt = 0
    clocks = 0
    # capture data for labs for given mjd, v, D, vec
    for lab in labs:
        #print('Lab: ',lab)
        sh = vmul(vec,[ inf[lab]['X'],inf[lab]['Y'],inf[lab]['Z'] ]) / v
        shmjd = sh/86400
        s = d[lab].getrange(mjd+shmjd,end_mjd+shmjd)
        
        
        if s != None and len(s.dtab)==1:
            if s.dtab[0].mjd_tab[-1]-s.dtab[0].mjd_tab[0] >= 0.95*durm:
                #s.rm_dc()
                s.rm_drift_each()
                #s.plot()
                datx.append(s.mjd_tab() - (mjd+shmjd) + lnum[lab])
                daty.append(s.val_tab())
                cnt = cnt+1
                clocks = clocks + (1 << lnum[lab])
                #sd[lnum[lab]]=s.std()

    if cnt>2:
            datx = np.concatenate(datx)
            daty = np.concatenate(daty)
            #sig = sigf(datx)
            try:
                popt, pcov = scp.curve_fit(fu, datx, daty )
                #print('POPT:  ', popt)
                #print('PCOV:  ', pcov)
                #plt.plot([x-int(x) for x in datx], daty, 
                #        [x-int(x) for x in datx], fu(datx,*popt), 
                        #[x-int(x) for x in datx], sig
                #        )
                #plt.plot(datx, daty, datx, fu(datx,*popt))
                #plt.show()
            except:
                return None
            return [popt[0], pcov[0,0]**0.5, pcov[1,1]**0.5, clocks]
    else:
            return None

#creating data file

v = 300000  # m/s
etaum = 20/86400
d = dict()
for lab in labs:
    d[lab] = tls.MTSerie(lab, color=inf[lab]['col'])
    mjd = np.arange(58658,58658.006,0.00001)
    sh = vmul([1,0,0],[ inf[lab]['X'],inf[lab]['Y'],inf[lab]['Z'] ]) / v
    shmjd = sh/86400
    val = fu(mjd-58658-0.001-shmjd,1,0)
    d[lab].add_mjdf_data(mjd, val)
    d[lab].plot(show=0)
    #d[lab].alphnorm(atom=inf[lab]['atom'])  #convert AOM freq to da/a
plt.show()

# loop prameters----------------------
v = 300000  # m/s

vecs = [   [-1,-1,0], [-1,0,0], [-1,1,0], 
           [0,-1,0], [0,1,0], [1,-1,0], [1,0,0], [1,1,0]  ]
#vecs = [[-1,-1,-1]]
vecs = vecs[::-1]
#Ds = [50*v]
#Ds = [100*v]
Ds = [ 20*v]
#Ds = [20*v, 50*v]
#Ds = [150*v, 100*v]

mjds = np.arange(58658,58659 ,0.00005)  #co ok 4s
# loop -----------------------------------
for D in Ds:
    for vec in vecs:
        start = time.time()
        #print(D, vec)
        mjd_t = []
        out_t = []
        err_t = []
        err2_t = []
        clocks_t = []
        #mjds = np.arange(58658,58670 ,0.00001*(D/v)/5.)
        for mjd in mjds:
            w = calc_single(mjd, v, D, vec)
            #print(w)
            if w!=None:
                mjd_t.append(mjd)
                out_t.append(w[0])
                err_t.append(w[1])
                err2_t.append(w[2])
                clocks_t.append(w[3])
                #print(int(D/v), vec, mjd, ' -> ', w[0], w[1], w[3])
        #plt.plot(mjd_t, out_t, mjd_t, err_t)
        #plt.show()
        fname = 'D'+str(int(D/v))+'_V_'+str(vec[0])+'_'+str(vec[1])+'_'+str(vec[2])+'.npy'
        outdat = np.column_stack((  np.array(mjd_t),np.array(out_t),
                                    np.array(err_t),np.array(err2_t),
                                    np.array(clocks_t)  ))
        print(vec, np.max(np.abs(out_t)))
        # print(outdat[1][1])
        np.save(os.path.join(progspath,'DMAnaliza','code',
                'out','out04sym',fname),outdat)
        print('time [min]: ',(time.time()-start)/60.)
# ----------------------------------------

