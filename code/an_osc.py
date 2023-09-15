import sys
import os
from local_settings import progspath
sys.path.append(str(progspath / 'mytools'))
import matplotlib.pyplot as plt
import tools as tls
import numpy as np
import dm_tools as dmt
import scipy.optimize as scp
import time
import decimal as dec
import parameters as par

etaum = 0

# 0: K, 1: std
sd = [0]*len(par.labs)

def fu(dx,A,sh):
    return [f(x,A,sh) for x in dx]

def f(x,A,sh):   
    ts = 8640       #1/ (10s / 86400 )
    global etaum
    rx = int(x)
    return A*np.sin(etaum*(x-rx)+sh)

def sigf(dx):
    return [ sd[int(x)] for x in dx]

def vmul(a,b):
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]

def calc_single(mjd, om):
   
    global etaum
    global sd

    etau = 10/om  # rad/s
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


#Oms = [0.1, 0.01, 0.001]
Oms = [ 0.02, 0.002]
Oms = np.arange(0.001, 0.1, 0.01)
outs = []

mjds = np.arange(58658,58670 ,0.005)  #co ok 4s
for Om in Oms:
        start = time.time()
        mjd_t = []
        out_t = []
        err_t = []
        err2_t = []
        clocks_t = []
        for mjd in mjds:
            w = calc_single(mjd, Om)
            if w!=None:
                mjd_t.append(mjd)
                out_t.append(w[0])
                err_t.append(w[1])
                err2_t.append(w[2])
                clocks_t.append(w[3])
                # print(Om, mjd, ' -> ', w[0], w[1], w[3])
        # print(max(np.abs(out_t)), min(np.abs(out_t)) )
        
        Ts = 10
        tosc = 1./Om

        outs.append([Om, max(np.abs(out_t)),  min(np.abs(out_t)), ])
        npout = np.array(outs)
        np.save('osc.npy', npout)
        np.savetxt('osc.txt', npout)
        
        plt.plot(npout[:,0], npout[:,1]*1e-18)
        plt.yscale('log')
        plt.grid()
        plt.savefig('osc.png')

        print('time [min]: ',(time.time()-start)/60.)
# ----------------------------------------
