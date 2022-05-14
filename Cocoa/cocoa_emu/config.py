import yaml

class Config:
    def __init__(self, configfile):
        with open(configfile, "r") as stream:
            config_args = yaml.safe_load(stream)
        
        self.config_args_emu = config_args['emulator'] 
        self.params          = config_args['params'] 
        self.likelihood      = list(config_args['likelihood'].keys())[0]

        self.N_lhs        = int(self.config_args_emu['N_lhs'])
        self.savedir      = self.config_args_emu['savedir']
        self.N_train_iter = int(self.config_args_emu['N_train_iter'])
        
        self.priors = self.get_priors()
        self.N_dim = len(self.priors)
        
        self.param_labels = list(self.priors.keys())
    
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