from input_data import InputData
import parameters as par
from local_settings import progspath

path = str( progspath / (r'DMAnaliza/data/d_prepared/') )
print(path)
indat = InputData(campaigns=par.campaigns, labs=par.labs, inf=par.inf, path=path)
indat.load_data_from_raw_files()
indat.split(min_gap_s=12)
indat.rm_dc_each()
indat.high_gauss_filter_each(stddev=350)
indat.alphnorm()
T_s = 100
freq = 1/T_s
omega=2*3.14*freq
indat.add_sin(
    amplitude=100000,
    omega=omega,
)
indat.plot(savefig=False)
