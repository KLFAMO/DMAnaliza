"""
This code is modification of an_3d.py.
It uses speed of the Earth in galaxy calculated from earth_movement_tests_pm.py.
In this simulations we use random generated data with added pulse, which simulates defect in halo.
"""
from itertools import chain
import os

from local_settings import project_path

import numpy as np
import scipy.optimize as scp
import time

import parameters as par
from earth_movement_tests_pm import earth_velocity_fast
from input_data import InputData
import matplotlib.pyplot as plt

etaum = 0
default_inverse_ts = 1/(par.default_servo_time_s/86400)
simulation_noise_std = 1.0
amplitude = 2

# 0: K, 1: std
sd = [0]*len(par.labs)
dd = {}

def get_d():
    """load and prepare data"""
    path = str( project_path / (r'data/d_prepared/') )
    indat = InputData(campaigns=par.campaigns, labs=par.labs, inf=par.inf, path=path)
    # indat.load_data_from_raw_files()
    indat.generate_random_data(from_mjd=58000, to_mjd=58000.006, dt_s=1, mean_val=0, std_val=simulation_noise_std)
    indat.add_pulse(mjd=58000.002, amplitude=amplitude, size=5e6, 
                    vec=earth_velocity_fast(58000.0005)/np.linalg.norm(earth_velocity_fast(58000.0005)), 
                    speed=np.linalg.norm(earth_velocity_fast(58000.005)))
    # indat.rm_dc_each()
    # indat.high_gauss_filter_each(stddev=350)
    # indat.alphnorm()
    # indat.plot(file_name='indat2.png')
    dd = indat.get_measurement_data()
    return indat.get_data_dictionary()

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

def calc_single(p):
    """
    v - speed  [m/s]
    D - size [m]
    vec - direction
    """

    mjd = p['mjd']
    v = p['v']
    D = p['D']
    vec = p['vec']
    d = p['data']

    global etaum
    global sd

    vec_abs = np.sqrt(vec[0]**2+vec[1]**2+vec[2]**2)
    if vec_abs!=0:
        vec = vec/vec_abs
    else:
        return None
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
                    # sd[par.lnum[lab]]=s.std()

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


def calc_for_single_mjd(p):
    global d
    w = calc_single(p)
    wm = calc_single({'mjd':p['mjd'], 'D':p['D'], 'v':p['v'], 'vec':-p['vec'], 'data':d})
    if w!=None and wm!=None:
        return (p['mjd'], w[0], w[1], wm[0], wm[1])
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


if __name__ == "__main__":
    d = get_d()
    for lab in par.labs:
        if lab in d:
            sd[par.lnum[lab]] = simulation_noise_std
    maxvs = []
    time_all_start = time.time()
    mjd_ranges = list(par.mjds_dict.values())
    mjds_chain = list(chain.from_iterable(mjd_ranges))
    
    for D in par.Ds:

        # start = time.time()
        params = [{
                'mjd':mjd,  # mjd of defect start
                'D':D,   # size of defect [m]
                'v':np.linalg.norm(earth_velocity_fast(mjd)),  # use speed of the Earth in galaxy
                'vec':earth_velocity_fast(mjd)/np.linalg.norm(earth_velocity_fast(mjd)),  # use direction of the Earth in galaxy
                'data':d,
            } for mjd in mjds_chain]
        
        out = [calc_for_single_mjd(p) for p in params]
        out = [ x for x in out if x!=None]

    
    # check if ../out exists, if not create it
    out_path = project_path / 'out'
    if not out_path.exists():
        out_path.mkdir()
    print(project_path)

    f = open(out_path / 'time.dat', 'a')
    f.write(f"\n{(time.time()-time_all_start)/60.} min")
    f.close()

    mjd = [i[0] for i in out]
    y = [i[1] for i in out]
    ym = [i[3] for i in out]

    

    out_arr = np.asarray(out, dtype=float)

    if out_arr.ndim != 2 or out_arr.shape[1] < 5:
        raise ValueError(
            "Nieprawidłowy format out. Oczekiwano kolumn: "
            "mjd, A_plus, sigma_plus, A_minus, sigma_minus."
        )

    # ========================================================
    # Odczytanie wyników analizy
    # ========================================================

    mjd_scan = out_arr[:, 0]

    A_plus = out_arr[:, 1]
    sigma_plus = out_arr[:, 2]

    A_minus = out_arr[:, 3]
    sigma_minus = out_arr[:, 4]

    # ========================================================
    # Obliczenie rho+ i rho-
    # ========================================================

    rho_plus = np.full_like(A_plus, np.nan)
    rho_minus = np.full_like(A_minus, np.nan)

    valid_plus = (
        np.isfinite(A_plus)
        & np.isfinite(sigma_plus)
        & (sigma_plus > 0)
    )

    valid_minus = (
        np.isfinite(A_minus)
        & np.isfinite(sigma_minus)
        & (sigma_minus > 0)
    )

    rho_plus[valid_plus] = (
        A_plus[valid_plus]
        / sigma_plus[valid_plus]
    )

    rho_minus[valid_minus] = (
        A_minus[valid_minus]
        / sigma_minus[valid_minus]
    )

    # ========================================================
    # Kontrast kierunkowy
    #
    # Cdir =  1: preferowany kierunek +v
    # Cdir = -1: preferowany kierunek -v
    # Cdir =  0: brak wyraźnej preferencji
    # ========================================================

    rho_plus_sq = rho_plus**2
    rho_minus_sq = rho_minus**2

    contrast_denominator = (
        rho_plus_sq + rho_minus_sq
    )

    Cdir = np.full_like(rho_plus, np.nan)

    valid_contrast = (
        np.isfinite(contrast_denominator)
        & (contrast_denominator > 0)
    )

    Cdir[valid_contrast] = (
        rho_plus_sq[valid_contrast]
        - rho_minus_sq[valid_contrast]
    ) / contrast_denominator[valid_contrast]

    # ========================================================
    # Korekta czasu analizy
    #
    # mjd_scan oznacza początek okna. Funkcja modelowa osiąga
    # maksimum po około 3 * etaum.
    # ========================================================

    speed_arr = np.asarray([
        np.linalg.norm(earth_velocity_fast(t))
        for t in mjd_scan
    ])

    etaum_arr = D / speed_arr / 86400

    mjd_result = (
        mjd_scan + 3 * etaum_arr
    )

    # ========================================================
    # Przygotowanie surowych danych
    # ========================================================

    raw_mjd_by_lab = {}

    for lab in par.labs:
        if lab in d:
            lab_mjd = np.asarray(
                d[lab].mjd_tab(),
                dtype=float,
            )

            if lab_mjd.size > 0:
                raw_mjd_by_lab[lab] = lab_mjd

    if not raw_mjd_by_lab:
        raise ValueError(
            "Brak surowych danych do narysowania."
        )

    # ========================================================
    # Wyznaczenie wspólnego zakresu czasu
    #
    # Pokazujemy tylko zakres, w którym dostępne są zarówno
    # surowe dane, jak i wyniki analizy.
    # ========================================================

    raw_mjd_min = min(
        np.min(lab_mjd)
        for lab_mjd in raw_mjd_by_lab.values()
    )

    raw_mjd_max = max(
        np.max(lab_mjd)
        for lab_mjd in raw_mjd_by_lab.values()
    )

    analysis_mjd_min = np.nanmin(mjd_result)
    analysis_mjd_max = np.nanmax(mjd_result)

    common_start_mjd = max(
        raw_mjd_min,
        analysis_mjd_min,
    )

    common_end_mjd = min(
        raw_mjd_max,
        analysis_mjd_max,
    )

    if common_end_mjd <= common_start_mjd:
        raise ValueError(
            "Surowe dane i wyniki analizy nie mają "
            "wspólnego zakresu czasu."
        )

    # Czas w sekundach od początku wspólnego zakresu.
    x_result_s = (
        mjd_result - common_start_mjd
    ) * 86400

    duration_s = (
        common_end_mjd - common_start_mjd
    ) * 86400

    # ========================================================
    # Rysowanie trzech paneli
    # ========================================================

    fig, (ax1, ax2, ax3) = plt.subplots(
        3,
        1,
        sharex=True,
        figsize=(14, 10),
        gridspec_kw={
            'height_ratios': [1.4, 1.0, 1.0],
        },
    )

    # --------------------------------------------------------
    # Panel 1: surowe dane z laboratoriów
    # --------------------------------------------------------

    for lab, lab_mjd in raw_mjd_by_lab.items():
        lab_time_s = (
            lab_mjd - common_start_mjd
        ) * 86400

        lab_value = np.asarray(
            d[lab].val_tab()
        )

        ax1.scatter(
            lab_time_s,
            lab_value,
            s=2,
            label=lab,
        )

    ax1.set_ylabel('Sygnał sensora')
    ax1.set_title(
        'Surowe, symulowane dane pomiarowe'
    )

    ax1.legend(
        ncol=5,
        fontsize=8,
        markerscale=3,
    )

    ax1.grid(alpha=0.25)

    # --------------------------------------------------------
    # Panel 2: rho+ i rho-
    # --------------------------------------------------------

    ax2.plot(
        x_result_s,
        rho_plus,
        color='tab:blue',
        label=r'$\rho_{+}$: kierunek $\mathbf{v}$',
    )

    ax2.plot(
        x_result_s,
        rho_minus,
        color='tab:orange',
        label=r'$\rho_{-}$: kierunek $-\mathbf{v}$',
    )

    ax2.axhline(
        0,
        color='black',
        linewidth=0.7,
    )

    ax2.set_ylabel(
        r'$\rho=A/\sigma_A$'
    )

    ax2.set_title(
        'Znormalizowane amplitudy dopasowania'
    )

    ax2.legend()
    ax2.grid(alpha=0.25)

    # --------------------------------------------------------
    # Panel 3: kontrast kierunkowy
    # --------------------------------------------------------

    ax3.plot(
        x_result_s,
        Cdir,
        color='tab:purple',
        linewidth=1.2,
        label=r'$C_{\mathrm{dir}}$',
    )

    ax3.axhline(
        0,
        color='black',
        linewidth=0.7,
    )

    ax3.set_ylim(-1.05, 1.05)
    ax3.set_xlim(0, duration_s)

    ax3.set_xlabel(
        'Czas [s]'
    )

    ax3.set_ylabel(
        r'$C_{\mathrm{dir}}$'
    )

    ax3.set_title(
        'Kontrast kierunkowy'
    )

    ax3.legend()
    ax3.grid(alpha=0.25)

    # Usunięcie dodatkowych marginesów poziomych.
    for ax in (ax1, ax2, ax3):
        ax.margins(x=0)

    plt.tight_layout(
        rect=[0, 0, 1, 0.95]
    )

    plt.show()