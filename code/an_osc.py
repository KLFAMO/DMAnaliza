import sys
import os
from local_settings import progspath
sys.path.append(str(progspath / 'mytools'))
import matplotlib.pyplot as plt
import tools as tls
import numpy as np
import scipy.optimize as scp
import time
import multiprocessing
import parameters as par

etaum = 0

# 0: K, 1: std
sd = [0]*len(par.labs)

def fu(dx,A,sh):
    return [f(x,A,sh) for x in dx]

def f(x,A,sh):   
    global etaum
    rx = int(x)
    return A*np.sin(etaum*(x-rx)+sh)*ssf_osc(Om)

def sigf(dx):
    return [ sd[int(x)] for x in dx]

def ssf_osc(om_rad, Ts = 10):
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
    global etaum
    global sd

    etau = om  # rad/s
    etaum = etau*86400
    end_mjd = mjd+0.05
    durm = end_mjd-mjd
    datx=[]
    daty=[]

    cnt = 0
    clocks = 0
    # capture data for labs for given mjd, om
    for lab in par.labs:
        if lab in d.keys():
            s = d[lab].getrange(mjd,end_mjd)
            
            if s != None and len(s.dtab)==1:
                if s.dtab[0].mjd_tab[-1]-s.dtab[0].mjd_tab[0] >= 0.95*durm:
                    #s.rm_dc()
                    s.rm_drift_each()
                    #s.plot()
                    datx.append(s.mjd_tab() - (mjd) + par.lnum[lab])
                    daty.append(s.val_tab())
                    cnt = cnt+1
                    clocks = clocks + (1 << par.lnum[lab])
                    sd[par.lnum[lab]]=s.std()

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

#reading data from npy files
d = dict()
for lab in par.labs:
    print('\n'+lab)
    print(os.path.isfile( str( progspath / (r'DMAnaliza/data/d_prepared/d_' +lab+'_'+par.camp+'.npy') ) ))
    if os.path.isfile( str( progspath / (r'DMAnaliza/data/d_prepared/d_' +lab+'_'+par.camp+'.npy') ) ):
        #print('\n'+lab)
        d[lab] = tls.MTSerie(lab, color=par.inf[lab]['col'])
        d[lab].add_mjdf_from_file(
            str( progspath / (r'DMAnaliza/data/d_prepared/d_' +lab+'_'+par.camp+'.npy') )   )
        d[lab].split(min_gap=12)
        #d[lab].rm_dc_each()
        #d[lab].high_gauss_filter_each(stddev=350)
        #d[lab].rm_drift_each()
        d[lab].alphnorm(atom=par.inf[lab]['atom'])  #convert AOM freq to da/a
        #sd[lnum[lab]]=d[lab].std()

def calc_for_single_mjd(p):
    w = calc_single(p['mjd'], p['Om'])
    if w!=None:
        return (p['mjd'], w[0], w[1], w[2], w[3])
    else:
        return None
    
#Oms = [0.1, 0.01, 0.001]
Oms = [ 0.02, 0.002]
Oms = np.arange(0.001, 0.2, 0.005)
outs = []

for Om in Oms:
        start = time.time()
        params = [{'mjd':mjd, 'Om':Om} for mjd in par.mjds]
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
        np.save('osc.npy', npout)
        np.savetxt('osc.txt', npout)
        
        plt.clf()
        plt.plot(npout[:,0], npout[:,1]*1e-18)
        plt.plot(npout[:,0], npout[:,2]*1e-18)
        plt.yscale('log')
        plt.grid()
        plt.savefig('osc.png')

        print('time [min]: ',(time.time()-start)/60.)
# ----------------------------------------
