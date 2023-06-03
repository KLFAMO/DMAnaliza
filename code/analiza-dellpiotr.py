import sys
sys.path.insert(1, '/home/piotr/progs/mytools/')

import tools as tls

labs = ['UMK1','UMK2', 'NICT', 'NPLYb', 'NIST', 'NPLSr' ]
#labs = ['UMK1']
tcol = {'UMK1':'green', 
        'UMK2':'red',
        'NIST':'blue',
        'NPLSr':'cyan',
        'NPLYb':'black',
        'NICT':'gray'}

d = dict()

for lab in labs:
    print('\n'+lab)
    d[lab] = tls.MTSerie(lab, color=tcol[lab])
    d[lab].add_mjdf_from_file('d_'+lab+'.npy')
    d[lab].split(min_gap=12)
    d[lab].rm_dc_each()
    d[lab].high_gauss_filter_each()
    d[lab].plot()
