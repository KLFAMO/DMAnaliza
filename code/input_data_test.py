from input_data import InputData
import parameters as par

indat = InputData(campaigns=par.campaigns, labs=par.labs, inf=par.inf)
indat.generate_random_data(from_mjd=58000, to_mjd=58000.1, dt_s=2, mean_val=0, std_val=1)
indat.print_info()