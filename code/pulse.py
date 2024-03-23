import logging
import numpy as np
import parameters as par
import timanda.tserie as tls

default_inverse_ts = 1/(par.default_servo_time_s/86400)

def vmul(a,b):
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]

def generate_mts_pulse(lab, mjd, amplitude, size, vec, speed):
    A=amplitude
    v=speed
    D=size
    vec_abs = np.sqrt(vec[0]**2+vec[1]**2+vec[2]**2)
    vec = vec/vec_abs
    etau = D/v
    etaum = etau/86400 # defect duration [mjd]
    end_mjd = mjd+2*etaum + 7/8640
    durm = end_mjd-mjd  # analysis window duration [mjd]
    sh = vmul(vec,[ par.inf[lab]['X'],par.inf[lab]['Y'],par.inf[lab]['Z'] ]) / v
    shmjd = sh/86400     #calculate mjd shift (delay) for given lab
    mjd_tab = np.arange(start=mjd+shmjd, stop=end_mjd+shmjd, step=0.00001)
    val_tab = np.array([servo_response(x-mjd-shmjd, A, etaum) for x in mjd_tab])
    mts = tls.MTSerie(TSerie=tls.TSerie(mjd=mjd_tab, val=val_tab))
    return mts

def servo_response(x, A, etaum):
    inverse_ts = default_inverse_ts
    if x < etaum:
        return A*(1-np.exp(-x*inverse_ts))   
    else:
        return ( A*(1-np.exp(-etaum*inverse_ts)) * np.exp(-(x-etaum)*inverse_ts) )