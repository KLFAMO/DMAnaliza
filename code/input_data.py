import logging
import numpy as np
import os
import pulse
# import tools as tls
import timanda.tserie as tls
import matplotlib.pyplot as plt

logging.basicConfig(
    level=logging.INFO, format='%(levelname)s - %(message)s'
)

class InputData:

    def __init__(self, campaigns=[], labs=[], inf=[], path=''):
        self.campaigns = campaigns
        self.labs = labs
        self.inf = inf
        self.path = path
        self.d = dict()
        self.loaded_labs = []
    
    def load_data_from_raw_files(self):
        for lab in self.labs:
            logging.info("-------------------")
            logging.info(f"lab: {lab}")
            for campaign in self.campaigns:
                lab_path = self.path+'/d_'+lab+'_'+campaign+'.npy'
                logging.info(f"lab_path: {lab_path}")
                is_lab_file = os.path.isfile(lab_path)
                logging.info(f"is_lab_file: {is_lab_file}")
                if is_lab_file:
                    if lab not in self.loaded_labs:
                        self.loaded_labs.append(lab)
                        self.d[lab] = tls.MTSerie(lab, color=self.inf[lab]['col'])
                    self.d[lab].add_mjdf_from_file( lab_path )
    
    def split(self, min_gap_s=12):
        for lab in self.loaded_labs:
            self.d[lab].split(min_gap_s=min_gap_s)
    
    def rm_dc(self):
        for lab in self.loaded_labs:
            self.d[lab].rm_dc()

    def rm_dc_each(self):
        for lab in self.loaded_labs:
            self.d[lab].rm_dc_each()
    
    def rm_drift_each(self):
        for lab in self.loaded_labs:
            self.d[lab].rm_drift_each()

    def high_gauss_filter_each(self, stddev=350):
        for lab in self.loaded_labs:
            self.d[lab].high_gauss_filter_each(stddev=stddev)
    
    def rmoutlayers(self):
        for lab in self.loaded_labs:
            self.d[lab].rmoutlayers()

    def alphnorm(self):
        """
        convert AOM freq to da/a
        """
        for lab in self.loaded_labs:
            self.d[lab].alphnorm(atom=self.inf[lab]['atom'])

    def get_data_dictionary(self):
        return self.d
    
    def get_mjd_range(self, from_mjd, to_mjd):
        for lab in self.loaded_labs:
            self.d[lab].getrange_on_self(from_mjd, to_mjd)
    
    def plot(self, file_name='indata.png', savefig=True):
        for lab in self.loaded_labs:
            # self.d[lab].rm_dc_each()
            self.d[lab].plot(show=0)
        if savefig:
            plt.savefig(file_name)
        else:
            plt.show()

    def add_pulse(self, mjd, amplitude, size, vec, speed):
        """
        Add artificial pulse to existing data.

        params:
            mjd - event mjd
            amplitude - amplitude of the pulse
            direction - numpy 3d vector of the defect speed
            size - size of the defect
        """
        # direction_ampl = np.linalg.norm(direction)
        # defect_duration = size/speed
        mjd_tab = [58666.0001, 58666.2, 58666.4, 58666.6, 58666.8 ,58667.1]
        val_tab = [220, 120, 220, 120, -220, 120]
        off_mts = tls.MTSerie(TSerie=tls.TSerie(mjd=mjd_tab, val=val_tab))
        for lab in self.loaded_labs:
            off_mts = pulse.generate_mts_pulse(
                mjd=mjd,
                lab=lab,
                amplitude=amplitude,
                size=size,
                vec=vec,
                speed=speed,
            )
            self.d[lab].add_val_offset_from_mts(off_mts)
    
    def add_sin(self, ref_mjd=0, amplitude=1, omega=0):
        """
        Add artificial sinusoidal signal to existing data.

        params:
            amplitude - amplitude A in A*sin(omega*t)
        """
        for lab in self.loaded_labs:
            self.d[lab].add_sin(ref_mjd=0, amplitude=amplitude, omega=omega)
