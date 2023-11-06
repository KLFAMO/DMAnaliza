import numpy as np

labs = ['UMK1','UMK2', 'NIST']

inf = { 'JLAB1': {'col':'green', 'atom':'88Sr', 'X':3144873,  'Y':1928649,  'Z':3679736, 'servo_time_s': 20}, 
        'JLAB2': {'col':'red',   'atom':'88Sr', 'X':3949271,  'Y':5224649,  'Z':7071336, 'servo_time_s': 20},
}

lnum = {'JLAB1':0, 'JLAB2':1}

# calc loop parameters
v = 300000  # m/s   - speed of the Earth in space

vecs = [ [1,1,1], ]

Ds = [ 10*v, 20*v, 30*v]
#Ds = [ x*v for x in range(10,201)]

campaigns = ['c1', 'c2', 'c3']
campaigns = ['c1']

mjds_dict_fast = {
    'c1' : np.arange(58658,58659 ,0.03),  # for fast tests
    'c2' : np.arange(58917.8,58919 ,0.05),  # for fast tests
}

mjds_dict ={
    'c1' : np.arange(58658,58670 ,0.00005),  #every 4s
    'c2' : np.arange(58916,58935 ,0.00005),  #every 4s
    'c3' : np.arange(59638, 59670, 0.00005),
}

mjds_dict_osc ={
    'c1' : np.arange(58658,58670 ,0.005),  #every 400s
    'c2' : np.arange(58916,58935 ,0.005),  #every 400s
}

# mjds_dict = mjds_dict_fast

save_mjd_calcs = False

default_servo_time_s = 20
min_required_clocks = 2
expected_event_to_event_mjd = 0.1
# at the moment earth's speed vector is fixed to [1,1,1] because of some issues with function calculating it