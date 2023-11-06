#%%

import numpy as np 
import random 
import matplotlib.pyplot as plt
from scipy import signal

# generate mjd every 4s for 4min 

c1 = np.arange(58658, 58658.0027, 0.00005) 
c1 = list(c1)

# generate random frequency correction 

correction = random.triangular(-10, 10)
print("correction is:", correction)

def correction_simulator(c1):

    """ generates random numbers between -10 and 10 every 4s for 4min"""

    S = []
    for i in range(len(c1)):
        # random numbers between -10 and 10
        s = random.triangular(-10, 10)
        S.append(s)
    nps = np.array(S)

    return nps

correction = correction_simulator(c1)
c1 = np.array(c1)

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


umk1 = correction_simulator(c1)
np.save('d_UMK1_c1.npy', umk1)

umk2 = correction_simulator(c1)
np.save('d_UMK2_c1.npy', umk2)

umk3 = correction_simulator(c1)
np.save('d_NIST _c1.npy', umk3)
# %%
