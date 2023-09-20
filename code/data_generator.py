import numpy as np 
import random 

# generate mjd every 4s for 4min 

c1 = np.arange(58658, 58658.0027, 0.00005)
c1 = list(c1)

# generate random frequency correction 

correction = random.triangular(-10, 10)
print("correction is:", correction)

def correction_simulator(c1):
    S = []
    for i in range(len(c1)):
        s = [c1[i], random.triangular(-10, 10)]
        S.append(s)
    nps = np.array(S)

    return nps

umk1 = correction_simulator(c1)
np.save('d_UMK1_c1.npy', umk1)

umk2 = correction_simulator(c1)
np.save('d_UMK2_c1.npy', umk2)

umk3 = correction_simulator(c1)
np.save('d_NIST_c1.npy', umk3)