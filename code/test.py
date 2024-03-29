from input_data import InputData

import parameters as par
import pulse
from local_settings import progspath

v = 230000 # m/s - speed of the Earth in space
D = 10*v

path = str( progspath / (r'DMAnaliza/data/d_prepared/') )
print(path)
indat = InputData(campaigns=par.campaigns, labs=par.labs, inf=par.inf, path=path)
indat.load_data_from_raw_files()
indat.split(min_gap_s=12)
indat.rm_dc_each()
indat.high_gauss_filter_each(stddev=350)
indat.alphnorm()
indat.add_pulse(
    mjd=58666.045,
    amplitude=50000,
    size=D,
    vec=[1, 0, 0],
    speed=v,
)
indat.plot(savefig=False)
