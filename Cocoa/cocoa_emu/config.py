import yaml
import numpy as np
# from .sampling import get_starting_pos

class Config:
    def __init__(self, configfile):
        with open(configfile, "r") as stream:
            config_args = yaml.safe_load(stream)
        
        self.config_args_emu = config_args['emulator'] 
        self.params          = config_args['params'] 
        config_args_lkl = config_args['likelihood']
        
        self.savedir         = self.config_args_emu['io']['savedir']
        try:
            self.save_train_data = self.config_args_emu['io']['save_train_data'].lower()=='true'
        except:
            self.save_train_data = False
        try:
            self.save_intermediate_model = self.config_args_emu['io']['save_intermediate_model'].lower()=='true'
        except:
            self.save_intermediate_model = False
        try:
            self.n_pcas_baryon = self.config_args_emu['baryons']['n_pcas_baryon']
        except:
            self.n_pcas_baryon = 0

        self.dv_fid_path   = self.config_args_emu['training']['dv_fid']
        self.n_lhs         = int(self.config_args_emu['training']['n_lhs'])
        self.n_train_iter  = int(self.config_args_emu['training']['n_train_iter'])
        self.n_resample    = int(self.config_args_emu['training']['n_resample'])
        self.emu_type      = self.config_args_emu['training']['emu_type']
        self.batch_size    = int(self.config_args_emu['training']['batch_size'])
        self.n_epochs      = int(self.config_args_emu['training']['n_epochs'])
        
        self.config_data(config_args_lkl)
        
        self.n_emcee_walkers = int(self.config_args_emu['sampling']['n_emcee_walkers'])
        self.n_mcmc          = int(self.config_args_emu['sampling']['n_mcmc'])
        self.n_burn_in       = int(self.config_args_emu['sampling']['n_burn_in'])
        self.n_thin          = int(self.config_args_emu['sampling']['n_thin'])
        self.temper0          = float(self.config_args_emu['sampling']['temper0'])
        self.temper_increment = float(self.config_args_emu['sampling']['temper_increment'])
        
        self.lhs_minmax    = self.get_lhs_minmax()
        self.n_dim         = len(self.lhs_minmax)
        
        self.param_labels = list(self.lhs_minmax.keys())
    
    def config_data(self, config_args_lkl):
        self.likelihood      = list(config_args_lkl.keys())[0]
        self.config_args_lkl = config_args_lkl[self.likelihood]
        self.likelihood_path = self.config_args_lkl['path']
        self.dataset_path    = self.likelihood_path + '/' + self.config_args_lkl['data_file']

        with open(self.dataset_path, 'r') as f:
            for line in f.readlines():
                split_line = line.split()
                if(split_line[0]=='data_file'):
                    self.dv_obs_path = self.likelihood_path + '/' + split_line[-1]
                if(split_line[0]=='cov_file'):
                    cov_file        = self.likelihood_path + '/' + split_line[-1]
                if(split_line[0]=='mask_file'):
                    mask_file       = self.likelihood_path + '/' + split_line[-1]
                if(split_line[0]=='baryon_pca_file'):
                    baryon_pca_file = self.likelihood_path + '/' + split_line[-1]
                if(split_line[0]=='source_ntomo'):
                    self.source_ntomo = int(split_line[-1])
                if(split_line[0]=='lens_ntomo'):
                    self.lens_ntomo = int(split_line[-1])
        
        self.baryon_pcas = np.loadtxt(baryon_pca_file)
        self.mask        = np.loadtxt(self.config_args_emu['sampling']['scalecut_mask'])[:,1].astype(bool)
        self.dv_fid      = np.loadtxt(self.dv_fid_path)[:,1]
        self.dv_obs      = np.loadtxt(self.dv_obs_path)[:,1]
        self.output_dims = len(self.dv_obs)
        self.shear_calib_mask = np.load(self.config_args_emu['shear_calib']['mask'])
        
        assert len(self.dv_obs)==len(self.dv_fid),"Observed data vector is of different size compared to the fiducial data vector."
        self.cov         = self.get_full_cov(cov_file)
        self.dv_std      = np.sqrt(np.diagonal(self.cov))
        
    def get_lhs_minmax(self):
        lh_minmax = {}
        for x in self.params:
            if('prior' in self.params[x]):
                prior = self.params[x]['prior']
                if('dist' in prior):
                    loc   = prior['loc']
                    scale = prior['scale']
                    lh_min = loc - 4. * scale
                    lh_max = loc + 4. * scale
                else:
                    lh_min = prior['min']
                    lh_max = prior['max']
                lh_minmax[x] = {'min': lh_min, 'max': lh_max}
        return lh_minmax
    
    def get_full_cov(self, cov_file):
        print("Getting covariance...")
        full_cov = np.loadtxt(cov_file)
        cov = np.zeros((self.output_dims, self.output_dims))
        for line in full_cov:
            i = int(line[0])
            j = int(line[1])

            cov_g_block  = line[-2]
            cov_ng_block = line[-1]

            cov_ij = cov_g_block + cov_ng_block

            cov[i,j] = cov_ij
            cov[j,i] = cov_ij

        return cov
    
    