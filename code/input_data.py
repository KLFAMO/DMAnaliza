import os
import tools as tls
import matplotlib.pyplot as plt

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
            print('\nlab: ', lab)
            for campaign in self.campaigns:
                lab_path = self.path+'/d_'+lab+'_'+campaign+'.npy'
                print('lab_path: ', lab_path)
                is_lab_file = os.path.isfile(lab_path)
                print('is_lab_file: ', is_lab_file )
                if is_lab_file:
                    if lab not in self.loaded_labs:
                        self.loaded_labs.append(lab)
                        self.d[lab] = tls.MTSerie(lab, color=self.inf[lab]['col'])
                    self.d[lab].add_mjdf_from_file( lab_path )
    
    def split(self, min_gap=12):
        for lab in self.loaded_labs:
            self.d[lab].split(min_gap=min_gap)
    
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
            self.d[lab].rm_dc_each()
            self.d[lab].plot(show=0)
        if savefig:
            plt.savefig(file_name)
        else:
            plt.show()