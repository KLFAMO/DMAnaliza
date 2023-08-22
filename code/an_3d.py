import sys
import os
import pathlib as pa
progspath = pa.Path(__file__).absolute().parents[2]
sys.path.append(str(progspath / 'mytools'))
#import matplotlib.pyplot as plt
import tools as tls
import numpy as np
import dm_tools as dmt
import scipy.optimize as scp
import time
import decimal as dec
import multiprocessing

from parameters import inf, labs, lnum

etaum = 0

# 0: K, 1: std
sd = [0,0,0,0,0,0,0,0,0]

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
    v - speed  [m/s]
    D - size [m]
    vec - direction
    """
    global etaum
    global sd

    vec_abs = np.sqrt(vec[0]**2+vec[1]**2+vec[2]**2)
    if vec_abs!=0:
        vec = vec/vec_abs
    etau = D/v  # defect duration [s]
    etaum = etau/86400 # defect duration [mjd]
    end_mjd = mjd+2*etaum + 7/8640
    durm = end_mjd-mjd  # analysis window duration [mjd]
    datx=[]
    daty=[]

    cnt = 0
    clocks = 0
    # capture data for labs for given mjd, v, D, vec
    for lab in labs:
        if lab in d.keys():
            #print('Lab: ',lab)
            sh = vmul(vec,[ inf[lab]['X'],inf[lab]['Y'],inf[lab]['Z'] ]) / v
            shmjd = sh/86400     #calculate mjd shift (delay) for given lab
            s = d[lab].getrange(mjd+shmjd,end_mjd+shmjd)  #get shifted data
            
            #if data exist, add and make inital calculations
            if s != None and len(s.dtab)==1:
                if s.dtab[0].mjd_tab[-1]-s.dtab[0].mjd_tab[0] >= 0.95*durm:
                    s.rm_drift_each()
                    #s.plot()
                    datx.append(s.mjd_tab() - (mjd+shmjd) + lnum[lab])
                    daty.append(s.val_tab())
                    cnt = cnt+1
                    clocks = clocks + (1 << lnum[lab])
                    sd[lnum[lab]]=s.std()

    #if data from at least 3 labs are captured, fit data
    if cnt>2:
            datx = np.concatenate(datx)
            daty = np.concatenate(daty)
            sig = sigf(datx)
            try:
                popt, pcov = scp.curve_fit(fu, datx, daty,  
                            sigma=sig, absolute_sigma=True )
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

#reading data from npy files "data/d_prepared/d_x.npy"
camp = 'c2'
d = dict()
for lab in labs:
    #print('\n'+lab)
    if os.path.isfile( str( progspath / (r'DMAnaliza/data/d_prepared/d_' +lab+'_'+camp+'.npy') ) ):
        d[lab] = tls.MTSerie(lab, color=inf[lab]['col'])
        d[lab].add_mjdf_from_file(
            str( progspath / (r'DMAnaliza/data/d_prepared/d_' +lab+'_'+camp+'.npy') )   )
        d[lab].split(min_gap=12)
        #d[lab].rm_dc_each()
        d[lab].rm_drift_each()
        d[lab].high_gauss_filter_each(stddev=350)
        d[lab].alphnorm(atom=inf[lab]['atom'])  #convert AOM freq to da/a
        #sd[lnum[lab]]=d[lab].std()

# loop prameters----------------------
v = 300000  # m/s

vecs = [   [-1,-1,-1],
        [-1,-1,0], [-1,-1,1], [-1,0,-1], [-1,0,0], [-1,0,1], [-1,1,-1], [-1,1,0], [-1,1,1], 
      [0,-1,-1], [0,-1,0], [0,-1,1], [0,0,-1], [0,0,0], [0,0,1], [0,1,-1], [0,1,0], [0,1,1],
      [1,-1,-1], [1,-1,0], [1,-1,1], [1,0,-1], [1,0,0], [1,0,1], [1,1,-1], [1,1,0], [1,1,1]  ]
#vecs = [[-1,-1,-1]]
#vecs = vecs[::-1]
#Ds = [50*v]
Ds = [ 50*v, 100*v, 150*v]

mjds_dict ={
    'c1' : np.arange(58658,58670 ,0.00005),  #co ok 4s
    'c2' : np.arange(58916,58935 ,0.00005),  #co ok 4s
    'c2' : np.arange(58917.8,58917.85 ,0.00005),  #co ok 4s  temporary
}
mjds = mjds_dict[camp]

def calc_for_single_mjd(p):
    w = calc_single(p['mjd'], p['v'], p['D'], p['vec'])
    if w!=None:
        return (p['mjd'], w[0], w[1], w[2], w[3])
    else:
        return None

time_all_start = time.time()
for D in Ds:
    for vec in vecs:
        start = time.time()
        params = [{'mjd':mjd, 'D':D, 'v':v, 'vec':vec} for mjd in mjds]
        with multiprocessing.Pool() as pool:
            out = pool.map(calc_for_single_mjd, params)
        #out = [calc_for_single_mjd(p) for p in params]
        out = [ x for x in out if x!=None]       
        fname = 'D'+str(int(D/v))+'_V_'+str(vec[0])+'_'+str(vec[1])+'_'+str(vec[2])+'.npy'
        outdat = np.array(out)
        np.save(os.path.join(progspath,'DMAnaliza',
                    'out','out50'+camp+'_'+fname),outdat)
        print('time [min]: ',(time.time()-start)/60.)


f = open(os.path.join(progspath,'DMAnaliza',
            'out','time.dat'), 'a')
f.write(f"\n{(time.time()-time_all_start)/60.} min")
f.close()
# ----------------------------------------
