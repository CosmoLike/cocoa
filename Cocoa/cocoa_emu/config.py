import yaml
import numpy as np

class Config:
    def __init__(self, configfile):
        with open(configfile, "r") as stream:
            config_args = yaml.safe_load(stream)
        
        self.config_args_emu = config_args['emulator'] 
        self.params          = config_args['params'] 
        config_args_lkl = config_args['likelihood']
        self.config_lkl(config_args_lkl)

        self.N_lhs        = int(self.config_args_emu['N_lhs'])
        self.savedir      = self.config_args_emu['savedir']
        self.N_train_iter = int(self.config_args_emu['N_train_iter'])
        self.dv_fid_path  = self.config_args_emu['dv_fid']
        
        self.emu_type     = self.config_args_emu['emu_type']
        self.batch_size   = int(self.config_args_emu['batch_size'])
        self.n_epochs     = int(self.config_args_emu['n_epochs'])

        self.dv_fid       = np.loadtxt(self.dv_fid_path)[:,1]
        self.output_dims  = len(self.dv_fid)
        self.cov          = self.get_full_cov()
        self.dv_std       = np.sqrt(np.diagonal(self.cov))
        
        self.priors       = self.get_priors()
        self.N_dim        = len(self.priors)
        
        self.param_labels = list(self.priors.keys())
    
    def config_lkl(self, config_args_lkl):
        self.likelihood      = list(config_args_lkl.keys())[0]
        self.config_args_lkl = config_args_lkl[self.likelihood]
        self.likelihood_path = self.config_args_lkl['path']
        self.dataset_path    = self.likelihood_path + '/' + self.config_args_lkl['data_file']

        with open(self.dataset_path, 'r') as f:
            for line in f.readlines():
                split_line = line.split()
                if(split_line[0]=='cov_file'):
                    self.cov_file        = self.likelihood_path + '/' + split_line[-1]
                if(split_line[0]=='mask_file'):
                    self.mask_file       = self.likelihood_path + '/' + split_line[-1]
                if(split_line[0]=='baryon_pca_file'):
                    self.baryon_pca_file = self.likelihood_path + '/' + split_line[-1]
            
    def get_priors(self):
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
    
    def get_full_cov(self):
        print("Getting covariance...")
        full_cov = np.loadtxt(self.cov_file)
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