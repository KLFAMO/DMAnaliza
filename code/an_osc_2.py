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

Om = 0.1

x = np.arange(100,200,0.1)
y = fu(x, 1, 0)

print(fu(x, 1, 0))


path = str( progspath / (r'DMAnaliza/data/d_prepared/') )
indat = InputData(campaigns=par.campaigns, labs=par.labs, inf=par.inf, path=path)
indat.load_data_from_raw_files()
# indat.get_mjd_range(58666, 58667)
indat.rm_dc()
indat.high_gauss_filter_each(stddev=50)
indat.split(min_gap=20)
indat.rm_dc_each()
indat.rm_drift_each()
indat.high_gauss_filter_each(stddev=50)
indat.rm_dc_each()
indat.rm_drift_each()
indat.alphnorm()
indat.plot(savefig=False)
labs_data = indat.get_data_dictionary()

mt = []
vt = []
st = []
for i, lab in enumerate(labs_data):
    ld = labs_data[lab]
    # ld.time_shift_each(i*40000)
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
for om in np.arange(0.001, 0.3, 0.03):
    Om = om
    popt, pcov = scp.curve_fit(fu, cmt, cvt, sigma=cst, absolute_sigma=True)
    ot.append(Om)
    at.append(np.abs(popt[0])*1e-18)
    plt.clf()
    plt.plot(ot, at)
    plt.yscale('log')
    plt.grid()
    plt.savefig('osc2.png')
              
