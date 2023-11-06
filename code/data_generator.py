#%%

import numpy as np 
import random 
import matplotlib.pyplot as plt
from scipy import signal

def correction_simulator(amplitude, duration, mjd):

    """ generates random numbers between -10 and 10 every 4s for 4min"""

    # generate mjd every 4s for the duration in seconds 
    duration_mjd = duration * 0.000675
    c1 = np.arange(mjd, mjd+duration_mjd, 0.00005) 
    c1_list = list(c1)

    S = []
    for i in range(len(c1_list)):
        # random numbers between -10 and 10
        s = np.random.normal(0, amplitude)
        S.append(s)
    nps = np.array(S)

    return c1, nps

correction = correction_simulator(10, 4*60, 58658)[1]
print(correction)
mjd = 58658
c1 = correction_simulator(10, 4*60, 58658)[0]
print(c1)

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

    return combined_data

print(pulse(c1, correction))


JLAB1 = correction_simulator(c1)
np.save('d_JLAB1_c1.npy', JLAB1)

JLAB2 = correction_simulator(c1)
np.save('d_JLAB2_c1.npy', JLAB2)
# %%
