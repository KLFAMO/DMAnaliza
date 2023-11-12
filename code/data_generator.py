#%%

import numpy as np 
import random 
import matplotlib.pyplot as plt
from scipy import signal

def data_simulator(amplitude, duration, mjd):

    """ 
    Generates random numbers between -10 and 10 every 4s for 4min
    
    inputs: 
    amplitude: standard deviation of the stimulated data 
    duration: in seconds
    mjd: starting MJD in MJD
    
    outputs:
    c1: time array with a point every 4s in mjd
    nps: stimulated data 

    """

    # generate mjd every 4s for the duration in seconds 
    duration_mjd = duration * 0.000675 #change the duration from seconds to MJD
    c1 = np.arange(mjd, mjd+duration_mjd, 0.00005) 
    c1_list = list(c1)

    S = []
    for i in range(len(c1_list)):
        # random numbers between -10 and 10
        s = np.random.normal(0, amplitude)
        S.append(s)
    nps = np.array(S)

    return c1, nps

correction = data_simulator(10, 4*60, 58658)[1]
mjd = 58658
c1 = data_simulator(10, 4*60, 58658)[0]

# plot the generated data 

plt.plot(c1,correction)
plt.title("Generated data for 4min")
plt.xlabel("MJD")
plt.ylabel("Frequency correction (au)")
plt.show


def pulse(c1, correction):

    """ generates 20s pulses """

    #generating a 20s pulse
    pulse_duration = 0.00005 * 5  # 5*4 seconds
    duty_cycle = 0.5  # Adjust this value to control the pulse width
    pulse_signal = signal.square(2 * np.pi * pulse_duration * c1, duty=duty_cycle)

    print(pulse_signal)

    # Combine with your existing data
    existing_data = correction 

    print(existing_data)

    combined_data = existing_data + pulse_signal

    plt.figure()
    plt.plot(c1, combined_data)
    plt.title("Generated data with 20s pulse")
    plt.xlabel("MJD")
    plt.ylabel("Frequency correction (au)")
    plt.show

    return combined_data

JLAB1 = data_simulator(10, 4*60, 58658)
np.save('d_JLAB1_c1.npy', JLAB1)

JLAB2 = data_simulator(10, 4*60, 58658)
np.save('d_JLAB2_c1.npy', JLAB2)
# %%
