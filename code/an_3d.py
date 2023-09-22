from itertools import chain
import sys
import os

from local_settings import progspath
sys.path.append(str(progspath / 'mytools'))

import tools as tls
import numpy as np
import scipy.optimize as scp
import time
import multiprocessing

import parameters as par
from earth_movement import earth_velocity_xyz
from input_data import InputData
import matplotlib.pyplot as plt

etaum = 0
default_inverse_ts = 1/(par.default_servo_time_s/86400)

# 0: K, 1: std
sd = [0]*len(par.labs)

def fu(dx,A,sh):
    return [f(x,A,sh) for x in dx]

def f(x,A,sh):
    inverse_ts = default_inverse_ts
    global etaum
    rx = int(x)
    if x-rx<2*etaum:
        return sh
    elif x-rx<3*etaum:
        return A*(1-np.exp(-(x-rx-2*etaum)*inverse_ts))+sh   
    else:
        return ( A*(1-np.exp(-etaum*inverse_ts)) * np.exp(-(x-rx-3*etaum)*inverse_ts) ) + sh

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
    for lab in par.labs:
        if lab in d.keys():
            sh = vmul(vec,[ par.inf[lab]['X'],par.inf[lab]['Y'],par.inf[lab]['Z'] ]) / v
            shmjd = sh/86400     #calculate mjd shift (delay) for given lab
            s = d[lab].getrange(mjd+shmjd,end_mjd+shmjd)  #get shifted data
            
            #if data exist, add and make inital calculations
            if s != None and len(s.dtab)==1:
                if s.dtab[0].mjd_tab[-1]-s.dtab[0].mjd_tab[0] >= 0.95*durm:
                    s.rm_drift_each()
                    #s.plot()
                    datx.append(s.mjd_tab() - (mjd+shmjd) + par.lnum[lab])
                    daty.append(s.val_tab())
                    cnt = cnt+1
                    clocks = clocks + (1 << par.lnum[lab])
                    sd[par.lnum[lab]]=s.std()

    #if data from at least 3 labs are captured, fit data
    if cnt>=par.min_required_clocks:
        datx = np.concatenate(datx)
        daty = np.concatenate(daty)
        sig = sigf(datx)
        try:
            popt, pcov = scp.curve_fit(fu, datx, daty,  
                        sigma=sig, absolute_sigma=True )
        except:
            return None
        return [popt[0], pcov[0,0]**0.5, pcov[1,1]**0.5, clocks]
    else:
        return None

#reading data from npy files
path = str( progspath / (r'DMAnaliza/data/d_prepared/') )
indat = InputData(campaigns=par.campaigns, labs=par.labs, inf=par.inf, path=path)
indat.load_data_from_raw_files()
indat.plot(file_name='indata1.png')
indat.split(min_gap=12)
indat.rm_dc_each()
indat.high_gauss_filter_each(stddev=350)
indat.alphnorm()
indat.plot(file_name='indat2.png')
d = indat.get_data_dictionary()


def calc_for_single_mjd(p):
    w = calc_single(p['mjd'], p['v'], p['D'], p['vec'])
    if w!=None:
        return (p['mjd'], w[0], w[1], w[2], w[3])
    else:
        return None

def calc_results_for_length(out, D, length_mjd):
    """
    calculate results from campaign for single set of parameters

    """
    outarr = np.array(out)
    
    m = outarr[:,0]
    v = np.abs(outarr[:,1])

    last_mjd = m[0]
    mgap = (float(D/par.v)/86400)*0.5

    min_maxv = 1e6
    lenm = len(m)
    i=0
    while i < lenm-1:
        x = i
        start = m[i]
        last_mjd = start
        maxv = v[i]
        state = 1
        while m[x]-start < length_mjd:
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
    maxvs.append([D/par.v, min_maxv, length_mjd])

maxvs = []
time_all_start = time.time()
mjd_ranges = list(par.mjds_dict.values())
mjds_chain = list(chain.from_iterable(mjd_ranges))
for D in par.Ds:
    print('event length [s]: ', D/par.v)
    for vec in par.vecs:
        # start = time.time()
        params = [{'mjd':mjd, 'D':D, 'v':par.v, 'vec':earth_velocity_xyz(mjd)} for mjd in mjds_chain]
        with multiprocessing.Pool() as pool:
            out = pool.map(calc_for_single_mjd, params)
        # out = [calc_for_single_mjd(p) for p in params]
        out = [ x for x in out if x!=None]
        if out:
            calc_results_for_length(out, D, par.expected_event_to_event_mjd)
        if par.save_mjd_calcs: 
            fname = 'D'+str(int(D/par.v))+'_V_'+str(vec[0])+'_'+str(vec[1])+'_'+str(vec[2])+'.npy'
            outdat = np.array(out)
            np.save(os.path.join(progspath,'DMAnaliza', 'out', 'out50_'+fname), outdat)

out_maxvs = np.array(maxvs)

f = open(os.path.join(progspath,'DMAnaliza',
            'out','time.dat'), 'a')
f.write(f"\n{(time.time()-time_all_start)/60.} min")
f.close()

plt.clf()
plt.plot(out_maxvs[:,0],out_maxvs[:,1]*1e-18)
plt.yscale('log')
plt.grid()
plt.savefig('maxvs.png')
# ----------------------------------------
