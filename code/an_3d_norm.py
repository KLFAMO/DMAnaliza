import sys
import pathlib as pa

progspath = pa.Path(__file__).absolute().parents[2]
sys.path.append(str(progspath / 'mytools'))

import matplotlib.pyplot as plt
import tools as tls

#labs = ['UMK1','UMK2', 'NICT', 'NPLYb', 'NIST', 'NPLSr', 'SYRTE' ]
labs = ['UMK1','UMK2', 'NICT', 'NPLYb', 'NIST', 'NPLSr' ]
#labs = ['UMK1', 'UMK2']
#labs = ['SYRTE']
inf = { 'UMK1': {'col':'green', 'atom':'88Sr'}, 
        'UMK2': {'col':'red',   'atom':'88Sr'},
        'NIST': {'col':'blue',  'atom':'171Yb'},
        'NPLSr':{'col':'cyan',  'atom':'87Sr'},
        'NPLYb':{'col':'black', 'atom':'171Yb+'}, #sprawdzic
        'NICT': {'col':'gray',  'atom':'87Sr'},
        'SYRTE':{'col':'brown', 'atom':'87Sr'}
}

d = dict()

for lab in labs:
    print('\n'+lab)
    d[lab] = tls.MTSerie(lab, color=inf[lab]['col'])
    d[lab].add_mjdf_from_file(
           str( progspath / (r'DMAnaliza/data/d_prepared/d_' +lab+'.npy') )   )
    d[lab].split(min_gap=12)
    d[lab].rm_dc_each()
    d[lab].high_gauss_filter_each(stddev=100)
    #d[lab]+=100
    d[lab].alphnorm(atom=inf[lab]['atom'])
    d[lab].plot(show=0)
    

plt.show()
#a = tls.MGserie(d['UMK1'], d['UMK2'], grid_s=1)
