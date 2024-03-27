import sys
from local_settings import progspath
sys.path.append(str(progspath / 'mytools'))
import parameters as par
from input_data import InputData
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as scp

def ssf_osc(om_rad, Ts = par.default_servo_time_s):
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

def fu(dx,A,sh):
    return [f(x,A,sh) for x in dx]

def f(x,A,sh):   
    return A*np.sin(Om*x+sh)*ssf_osc(Om)

path = str( progspath / (r'DMAnaliza/data/d_prepared/') )
indat = InputData(campaigns=par.campaigns, labs=par.labs, inf=par.inf, path=path)
indat.load_data_from_raw_files()
indat.split(min_gap_s=20)
indat.rm_dc_each()
indat.rm_drift_each()
indat.high_gauss_filter_each(stddev=350)
indat.alphnorm()
indat.plot(savefig=True, file_name=f'indat_{par.campaigns[0]}.png')
labs_data = indat.get_data_dictionary()

mt = []
vt = []
st = []
for i, lab in enumerate(labs_data):
    ld = labs_data[lab]
    mt.append( ld.mjd_tab() )
    vt.append( ld.val_tab() )
    std = ld.std()
    sld = ld*0+std
    st.append( sld.val_tab() )

cmt = np.concatenate( mt )
cvt = np.concatenate( vt )
cst = np.concatenate( st )
cst = list(cst)

Om = 0.1
ot = []
at = []
for om in np.arange(0.01, 0.31, 0.01):
    Om = om
    popt, pcov = scp.curve_fit(fu, cmt, cvt, sigma=cst, absolute_sigma=True)
    ot.append(Om)
    at.append(np.abs(popt[0])*1e-18)
    plt.clf()
    plt.plot(ot, at)
    plt.yscale('log')
    plt.grid()
    plt.savefig(f'osc2_{par.campaigns[0]}.png')
    ff = open(f'dat_{par.campaigns[0]}.dat', 'a')
    ff.write(f'{om}\t{np.abs(popt[0])}\n')
    ff.close()
              
