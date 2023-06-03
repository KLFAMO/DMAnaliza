import sys
import pathlib as pa

progspath = pa.Path(__file__).absolute().parents[2]
sys.path.append(str(progspath / "mytools"))

import matplotlib.pyplot as plt
import tools as tls

labs = ["UMK1", "UMK2", "NICT", "SYRTE", "NPLYb", "NIST", "NPLSr"]
#labs = ['UMK1', 'UMK2']
#labs = ['SYRTE' ]
tcol = {
    "UMK1": "green",
    "UMK2": "red",
    "NIST": "blue",
    "NPLSr": "cyan",
    "NPLYb": "black",
    "NICT": "gray",
    "SYRTE": "brown",
}

d = dict()

for lab in labs:
    print("\n" + lab)
    d[lab] = tls.MTSerie(lab, color=tcol[lab])
    d[lab].add_mjdf_from_file(
        str(progspath / (r"DMAnaliza/data/d_prepared/d_" + lab + "_c1.npy"))
    )
    d[lab].rmrange(58669.42, 60000)
    d[lab].split(min_gap=12)
    d[lab].rm_dc_each()

    d[lab].high_gauss_filter_each()
    d[lab].rm_drift_each()
    d[lab] += 100
    d[lab].plot(show=0)
    print(d[lab].std())
plt.show()
# a = tls.MGserie(d['UMK1'], d['UMK2'], grid_s=1)
