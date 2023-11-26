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

#c1 = list(c1)

def nearest_index(array, value):
    # calculate the difference array
    difference_array = np.absolute(array-value)
 
    # find the index of minimum element from the array
    index = difference_array.argmin()
    return index

def pulse(c1, correction, pulse_duration, pulse_amplitude, starting):

    """ 
    Generates a pulse
    
    input: 
    c1: time array with a point every 4s in mjd
    correction: stimulated data for (amplitude, duration, mjd)=(10, 4*60, 58658)
    pulse_duration: length of the pulse in seconds
    pulse_amplitude: amplitude of the pulse
    starting: how long after the data begins should the pulse be, in seconds

    outputs:

       """

    # converting to MJD

    pulse_duration = 0.000675 * pulse_duration 
    starting = 0.000675 * starting 

    start = c1[0] + starting
    end = start + pulse_duration

    # generating a pulse
    
    print(c1[0], c1[1])
    print(start, pulse_duration, start+pulse_duration)

    pulse_begin = nearest_index(c1, start)
    pulse_end = nearest_index(c1, end)

    print(pulse_begin, pulse_end)

    existing_data = correction 

    print(existing_data)

    combined_data = existing_data

    i = 0

    while i < len(c1):
        if i in range(pulse_begin, pulse_end):
            combined_data[i] += pulse_amplitude
        else:
            combined_data[i] = existing_data[i]
        i += 1

    plt.figure()
    plt.plot(c1, combined_data)
    plt.title("Generated data with 20s pulse")
    plt.xlabel("MJD")
    plt.ylabel("Frequency correction (au)")
    plt.show

    return combined_data

results = pulse(c1, correction, 30, 30, 80)
print(results)

JLAB1 = data_simulator(10, 4*60, 58658)
np.save('d_JLAB1_c1.npy', JLAB1)

JLAB2 = data_simulator(10, 4*60, 58658)
np.save('d_JLAB2_c1.npy', JLAB2)
# %%
