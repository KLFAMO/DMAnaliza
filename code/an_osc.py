
from itertools import chain
from local_settings import progspath
import matplotlib.pyplot as plt
import multiprocessing
import numpy as np
import os
import parameters as par
import scipy.optimize as scp
import sys
import time
import tools as tls
from input_data import InputData
sys.path.append(str(progspath / 'mytools'))
etau_mjd = 0

# 0: K, 1: std
sd = [0]*len(par.labs)

def fu(dx,A,sh):
    return [f(x,A,sh) for x in dx]

def f(x,A,sh):   
    global etau_mjd
    rx = int(x)
    return A*np.sin(etau_mjd*(x-rx)+sh)*ssf_osc(Om)

def sigf(dx):
    return [ sd[int(x)] for x in dx]

def ssf_osc(om_rad, Ts = 20):
    """
    Calculates servo sensitivity factor for oscillations

    :param om_rad: angular frequency in radians
    :param Ts: servo time constant in seconds
    :return: amplitude of servo sensitivity factor
    
    :Changes:
        2023-09-15 (Piotr Morzyński): First version
    """
    f = om_rad/(2*np.pi)  # Hz
    tosc = 1./f
    s = 1j*2*np.pi/tosc
    H = 1/(Ts*s + 1)
    return np.abs(H)

def vmul(a,b):
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]

def calc_single(mjd, om):
    """
    Fits sine to data starting at mjd

    :param mjd: modified julian date when data start
    :param om: angular frequency of fited sine
    :return: [amplitude of sine from fit, pcov, pcov, clocks participated in fit]
    
    :Changes:
        2023-09-15 (Piotr Morzyński): First version
    """
    global etau_mjd
    global sd

    etau_s = 10*2*np.pi/om  # [s]
    etau_mjd = etau_s/86400.
    end_mjd = mjd+etau_mjd
    datx=[]
    daty=[]
    cnt = 0
    clocks = 0
    # capture data for labs for given mjd, om
    for lab in par.labs:
        if lab in labs_data.keys():
            data_serie = labs_data[lab].getrange(mjd,end_mjd)
            
            if data_serie != None and len(data_serie.dtab)==1:
                if etau_mjd/data_serie.getTotalTimeWithoutGaps() > 0.5:
                    data_serie.rm_drift_each()
                    datx.append(data_serie.mjd_tab() - (mjd) + par.lnum[lab])
                    daty.append(data_serie.val_tab())
                    cnt = cnt+1
                    clocks = clocks + (1 << par.lnum[lab])
                    sd[par.lnum[lab]]=data_serie.std()

    if cnt>2:
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
indat.split(min_gap=12)
#indat.high_gauss_filter_each(stddev=350)
indat.alphnorm()
d = indat.get_data_dictionary()


def calc_for_single_mjd(p):
    w = calc_single(p['mjd'], p['Om'])
    if w!=None:
        return (p['mjd'], w[0], w[1], w[2], w[3])
    else:
        return None
    
#Oms = [0.1, 0.01, 0.001]
Oms = [ 0.02, 0.002]
Oms = np.arange(0.001, 0.4, 0.005)
outs = []

mjd_ranges = list(par.mjds_dict.values())
mjds_chain = list(chain.from_iterable(mjd_ranges))

for Om in Oms:
        start = time.time()
        params = [{'mjd':mjd, 'Om':Om} for mjd in mjds_chain]
        # w = calc_single(mjd, Om)
        with multiprocessing.Pool() as pool:
            out = pool.map(calc_for_single_mjd, params)
        out = [ x for x in out if x!=None]   
        out = np.array(out) 
        print(out)       
        # Ts = 10
        tosc = 1./Om
        outs.append([Om, max(np.abs(out[:,1])),  min(np.abs(out[:,1])), ])
        npout = np.array(outs)
        np.save('osc_'+par.camp+'.npy', npout)
        np.savetxt('osc_'+par.camp+'.txt', npout)
        
        plt.clf()
        plt.plot(npout[:,0], npout[:,1]*1e-18)
        plt.yscale('log')
        plt.grid()
        plt.savefig('osc_'+par.camp+'.png')

        print('time [min]: ',(time.time()-start)/60.)
# ----------------------------------------
