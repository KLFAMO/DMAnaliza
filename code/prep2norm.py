import sys
import pathlib as pa

progspath = pa.Path(__file__).absolute().parents[2]
sys.path.append(str(progspath / 'mytools'))

import matplotlib.pyplot as plt
import tools as tls

#labs = ['UMK1','UMK2', 'NICT', 'SYRTE', 'NPLYb', 'NIST', 'NPLSr' ]
#labs = ['UMK1', 'UMK2']
labs = ['NPLYb', 'NIST', 'NPLSr' ]
tcol = {'UMK1':'green', 
        'UMK2':'red',
        'NIST':'blue',
        'NPLSr':'cyan',
        'NPLYb':'black',
        'NICT':'gray',
        'SYRTE':'brown'}

inf = { 'UMK1': {'col':'green', 'atom':'88Sr',
                 'X':3644273,  'Y':1226649,  'Z':5071736}, 
        'UMK2': {'col':'red',   'atom':'88Sr',
                 'X':3644273,  'Y':1226649,  'Z':5071736},
        'NIST': {'col':'blue',  'atom':'171Yb',
                 'X':-1288363, 'Y':-4721684, 'Z':4078659},
        'NPLSr':{'col':'cyan',  'atom':'87Sr',
                 'X':3985500,  'Y':-23625,   'Z':4962941},
        'NPLYb':{'col':'black', 'atom':'171Yb+',
                 'X':3985500,  'Y':-23625,   'Z':4962941},
        'NICT': {'col':'gray',  'atom':'87Sr',
                 'X':-3941931, 'Y':3368182,  'Z':3702068},
        'SYRTE':{'col':'brown', 'atom':'87Sr',
                 'X':4202777,  'Y':171368,   'Z':4778660}
}

d = dict()

for lab in labs:
    print('\n'+lab)
    d[lab] = tls.MTSerie(lab, color=tcol[lab])
    d[lab].add_mjdf_from_file(
           str( progspath / (r'DMAnaliza/data/d_prepared/d_' +lab+'.npy') )   )
    d[lab].rmrange(58669.42,60000)
    d[lab].split(min_gap=12)
    d[lab].rm_dc_each()
    d[lab].rm_drift_each()
    d[lab].high_gauss_filter_each(stddev=50)
    #d[lab]+=100
    d[lab].alphnorm(atom=inf[lab]['atom'])  #convert AOM freq to da/a
    d[lab].plot(show=0)
    print(d[lab].std())
plt.show()
#a = tls.MGserie(d['UMK1'], d['UMK2'], grid_s=1)
