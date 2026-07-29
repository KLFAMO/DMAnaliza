import logging
import numpy as np
import os
import pulse
# import tools as tls
import timanda.tserie as tls
from timanda.mtserie import MTSerie
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
                    self.loaded_labs.append(lab)
                    self.d[lab] = MTSerie(lab, color=self.inf[lab]['col'])
                    self.d[lab].add_mjdf_from_file( lab_path )
    
    def generate_random_data(self, from_mjd, to_mjd, dt_s=1, mean_val=0, std_val=1):
        """
        Generate random data for testing purposes.

        Parameters:
        from_mjd (float): Starting MJD
        to_mjd (float): Ending MJD
        dt_s (float): Sampling period in seconds (default 1)
        mean_val (float): Mean value of the random data (default 0)
        std_val (float): Standard deviation of the random data (default 1)
        """
        for lab in self.labs:
            logging.info("-------------------")
            logging.info(f"Generating random data for lab: {lab}")
            for campaign in self.campaigns:
                self.loaded_labs.append(lab)
                self.d[lab] = MTSerie.generate_random(from_mjd, to_mjd, dt_s, mean_val, std_val)

    
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
    
    # def plot(self, file_name='indata.png', savefig=True, split_horizontal=0):
    #     for i, lab in enumerate(self.loaded_labs):
    #         print(f"lab: {lab}, {i}")
    #         self.d[lab].plot(show=0, val_offset=split_horizontal*i)
    #     if savefig:
    #         plt.savefig(file_name)
    #     else:
    #         plt.show()

    def plot(self, file_name='indata.png', savefig=True, split_horizontal=0):
        """Plot all labs data in separated subplots."""
        num_labs = len(self.loaded_labs)
        fig, axes = plt.subplots(num_labs, 1, figsize=(10, 5 * num_labs), sharex=True)
        for i, lab in enumerate(self.loaded_labs):
            print(f"lab: {lab}, {i}")
            self.d[lab].plot(ax=axes[i], show=0)
            axes[i].set_title(f"Lab: {lab}")
        plt.xlabel("MJD")
        plt.tight_layout()
        if savefig:
            plt.savefig(file_name)
        else:
            plt.show()
       
    
    def print_info(self):
        for lab in self.loaded_labs:
            print(lab)
            print(self.d[lab])

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
        # mjd_tab = [58666.0001, 58666.2, 58666.4, 58666.6, 58666.8 ,58667.1]
        # val_tab = [220, 120, 220, 120, -220, 120]
        # off_mts = tls.MTSerie(TSerie=tls.TSerie(mjd=mjd_tab, val=val_tab))
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
    
    def add_sin(self, amplitude=1, omega=0):
        """
        Add artificial sinusoidal signal to existing data.

        params:
            amplitude - amplitude A in A*sin(omega*t)
        """
        for lab in self.loaded_labs:
            self.d[lab].add_sin(amplitude=amplitude, omega=omega)


    def get_measurement_data(self, labs=None):
        """
        Return measurement data for selected labs as NumPy arrays.

        Parameters:
        labs (iterable[str] or str, optional): Labs to return. If omitted,
            data for all loaded labs are returned.

        Returns:
        dict: ``{lab: {'mjd': np.ndarray, 'value': np.ndarray}}``.
            If an ``MTSerie`` contains several ``TSerie`` objects, their data
            are concatenated and sorted by MJD.

        Example:
            data = indat.get_measurement_data(['UMK', 'PTB'])
            for lab, series in data.items():
                plt.plot(series['mjd'], series['value'], label=lab)
            plt.legend()
        """
        if labs is None:
            labs = list(dict.fromkeys(self.loaded_labs))
        elif isinstance(labs, str):
            labs = [labs]
        else:
            labs = list(labs)

        unknown_labs = [lab for lab in labs if lab not in self.d]
        if unknown_labs:
            raise KeyError(
                f"No data loaded for labs: {unknown_labs}. "
                f"Available labs: {list(self.d)}"
            )

        result = {}
        for lab in labs:
            mts = self.d[lab]

            result[lab] = {'mjd': mts.mjd_tab(), 'value': mts.val_tab()}

        return result
