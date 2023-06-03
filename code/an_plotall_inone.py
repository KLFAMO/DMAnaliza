import sys
import pathlib as pa
import numpy as np
from os.path import exists

progspath = pa.Path(__file__).absolute().parents[2]
sys.path.append(str(progspath / "mytools"))

import matplotlib.pyplot as plt
import tools as tls

labs = [
    "UMK1",
    "UMK2",
    "NICT",
   "SYRTE",
    "NMIJ",
    "KRISS",
    "NPLSr",
    "NPLYb",
    "NIST",
    " ",
]
#labs = ['UMK2', 'UMK1']
#labs = [ "SYRTE"]

tcol = {
    "UMK1": "green",
    "UMK2": "red",
    "NIST": "blue",
    "NPLSr": "cyan",
    "NPLYb": "orange",
    "NICT": "gray",
    "SYRTE": "brown",
    "NMIJ": "purple",
    "KRISS": "green",
    " ": "black",
}

cx = ["c1","c2","c3"]

fig, tax = plt.subplots(1,3, gridspec_kw={'width_ratios': [12,18,30]})
#fig, tax = plt.subplots(1,2, gridspec_kw={'width_ratios': [12,30]})

for j in range(0, len(cx)):

    camp = cx[j]

    d = dict()
    i = 0
    for lab in labs:
        i += 1
        print("\n" + lab)
        d[lab] = tls.MTSerie(lab, color=tcol[lab])
        file_name = str(
            progspath / (r"DMAnaliza/data/d_prepared/d_" + lab + "_" + camp + ".npy")
        )
        if exists(file_name):
            d[lab].add_mjdf_from_file(file_name)
            d[lab].split(min_gap=12)
            #d[lab].rmrange(58670, 58915)
            # d[lab].rmrange(58934, 59915)
            d[lab].rm_dc_each()
            d[lab].high_gauss_filter_each()
            d[lab].rm_drift_each()
            #d[lab] += i*100
            d[lab].plot(show=0, ax=tax[j])

    #tax[j].set_yticks(np.arange(0, 3, 1))
    #ax = plt.gca()
    #ax = tax[j]
    #labels = [item.get_text() for item in tax[j].get_yticklabels()]
    #print(labels)
    #for i in range(0, len(labs)):
    #    print(i, labs[i])
    #    labels[i+1 ] = labs[i]
    #    if j>0:
    #        labels[i ] = ""
    #tax[j].set_yticklabels(labels)
    #ax.set_xlabel("time / MJD")

fig.tight_layout()
plt.show()
