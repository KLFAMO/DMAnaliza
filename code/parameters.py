import numpy as np

labs = ['UMK1','UMK2', 'NIST', 'SYRTE', 'NPLSr', 'NPLYb', 'NICT','NMIJ', 'KRISS', 'PTB']

inf = { 'UMK1': {'col':'green', 'atom':'88Sr', 'X':3644273,  'Y':1226649,  'Z':5071736, 'servo_time_s': 20}, 
        'UMK2': {'col':'red',   'atom':'88Sr', 'X':3644273,  'Y':1226649,  'Z':5071736, 'servo_time_s': 20},
        'NIST': {'col':'blue',  'atom':'171Yb', 'X':-1288363, 'Y':-4721684, 'Z':4078659, 'servo_time_s': 20},
        'NPLSr':{'col':'cyan',  'atom':'87Sr', 'X':3985500,  'Y':-23625,   'Z':4962941, 'servo_time_s': 20},
        'NPLYb':{'col':'black', 'atom':'171Yb+', 'X':3985500,  'Y':-23625,   'Z':4962941, 'servo_time_s': 20},
        'NICT': {'col':'gray',  'atom':'87Sr', 'X':-3941931, 'Y':3368182,  'Z':3702068, 'servo_time_s': 20},
        'SYRTE':{'col':'brown', 'atom':'87Sr', 'X':4202777,  'Y':171368,   'Z':4778660, 'servo_time_s': 20},
        'NMIJ' :{'col':'yellow', 'atom':'87Sr', 'X':-3953004,  'Y':3305232,   'Z':3758967, 'servo_time_s': 20},
        'KRISS':{'col':'brown', 'atom':'171Yb', 'X':-3116663,  'Y':4080538,   'Z':3783681, 'servo_time_s': 20},
        'PTB':{'col':'brown', 'atom':'87Sr', 'X':3836112,  'Y':708145,   'Z':5046077, 'servo_time_s': 20},
}

lnum = {'UMK1':0, 'UMK2':1, 'NIST':2, 'NPLSr':3, 'NPLYb':4,
        'NICT':5, 'SYRTE':6, 'NMIJ':7, 'KRISS':8, 'PTB':9}

# calc loop parameters
v = 230000  # m/s   - speed of the Earth in space

vecs = [ [1,1,1], ]

#Ds = [ 20*v, 50*v, 100*v, 150*v]
Ds = [ x*v for x in range(10,51)]

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

mjds_dict = mjds_dict_fast

save_mjd_calcs = False

default_servo_time_s = 20
min_required_clocks = 2
expected_event_to_event_mjd = 0.1
