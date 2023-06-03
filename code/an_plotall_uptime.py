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
    "NPLYb",
    "NIST",
    "NPLSr",
    "SYRTE",
    "NMIJ",
    "KRISS",
    " ",
]
# labs = ['UMK1', 'UMK2']
# labs = ["UMK1","NICT", "SYRTE"]

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

camp = "c3"


d = dict()
i = 0
for lab in labs:
    i += 1
    print("\n" + lab)
    d[lab] = tls.MTSerie(lab, color=tcol[lab])
    file_name = str(
        progspath / (r"DMAnaliza/data/raw/raw_" + lab + "_" + camp + ".npy")
    )
    if exists(file_name):
        d[lab].add_mjdf_from_file(file_name)
        d[lab].split(min_gap=12)
        d[lab].rmrange(58670, 58915)
        # d[lab].rmrange(58934, 59915)
        d[lab].rm_dc_each()
        d[lab].high_gauss_filter_each()
        d[lab] *= 0
        d[lab] += i
        d[lab].plot(show=0)

plt.yticks(np.arange(0, 11, 1))
ax = plt.gca()
labels = [item.get_text() for item in ax.get_yticklabels()]
for i, x in enumerate(labs):
    labels[i + 1] = labs[i]
# labels[1] = "FAMO Sr1"
# labels[2] = "FAMO Sr2"
# labels[3] = "NIST Sr"
# labels[4] = "NPL Sr"
# labels[5] = "NPL Yb+"
# labels[6] = "NICT Sr"
# labels[7] = "SYRTE Sr"
ax.set_yticklabels(labels)
ax.set_xlabel("time / MJD")
plt.show()
# a = tls.MGserie(d['UMK1'], d['UMK2'], grid_s=1)
