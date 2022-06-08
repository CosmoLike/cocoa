import yaml
import numpy as np
import torch

def hard_prior(theta, params_prior):
    is_lower_than_min = bool(np.sum(theta < params_prior[:,0]))
    is_higher_than_max = bool(np.sum(theta > params_prior[:,1]))
    if is_lower_than_min or is_higher_than_max:
        return -np.inf
    else:
        return 0.
    
def gaussian_prior(theta, params_prior):
    mu  = params_prior[:,0]
    std = params_prior[:,1]
    y = (theta - mu) / std
    return -0.5 * np.sum(y * y)

def split_with_comma(configline):
    configline_split = configline.split(',')
    configline_list = []
    for obj in configline_split:
        configline_list.append(float(obj))
    return np.array(configline_list)
    
class EmuSampler:
    def __init__(self, emu, config):
        self.emu               = emu
        self.params            = config.params
        
        self.probe             = config.probe
        
        self.n_walkers         = config.n_emcee_walkers
        self.n_pcas_baryon     = config.n_pcas_baryon
        self.baryon_pcas       = config.baryon_pcas
        
        self.emu_type          = config.emu_type
        self.mask              = config.mask
        self.cov               = config.cov
        self.masked_inv_cov    = np.linalg.inv(self.cov[self.mask][:,self.mask])
        self.dv_obs            = config.dv_obs
        self.shear_calib_mask  = config.shear_calib_mask
        
        self.source_ntomo      = config.source_ntomo
        if config.probe != 'cosmic_shear':
            self.lens_ntomo = config.lens_ntomo
        
        self.n_fast_pars       = config.n_fast_pars
        
        self.m_shear_prior_std = split_with_comma(config.config_args_emu['shear_calib']['prior_std'])
        self.config_args_baryons = config.config_args_emu['baryons']
        
        if self.probe!='cosmic_shear':
            self.config_args_bias  = config.config_args_emu['galaxy_bias']
            self.bias_fid          = split_with_comma(self.config_args_bias['bias_fid'])
            self.galaxy_bias_mask  = config.galaxy_bias_mask
        
        self.get_priors()
        
        if self.probe!='cosmic_shear':
            self.n_sample_dims    = config.n_dim + self.lens_ntomo + self.source_ntomo + config.n_pcas_baryon
        else:
            self.n_sample_dims    = config.n_dim + self.source_ntomo + config.n_pcas_baryon
        
    def get_priors(self):
        gaussian_prior_indices = []
        gaussian_prior_parameters = []

        flat_prior_indices    = []
        flat_prior_parameters = []

        ind = 0
        for x in self.params:
            if 'prior' in self.params[x]:
                prior = self.params[x]['prior']
                if('dist' in prior):
                    dist = prior['dist']
                    assert dist=='norm', "Got unexpected value"
                    gaussian_prior_indices.append(ind)
                    gaussian_prior_parameters.append([prior['loc'], prior['scale']])
                else:
                    flat_prior_indices.append(ind)
                    flat_prior_parameters.append([prior['min'], prior['max']])
                ind += 1

        self.flat_prior_indices     = flat_prior_indices
        self.gaussian_prior_indices = gaussian_prior_indices
        
        self.gaussian_prior_parameters = np.array(gaussian_prior_parameters)
        self.flat_prior_parameters     = np.array(flat_prior_parameters)
        
        if self.probe!='cosmic_shear':
            self.bias_prior        = split_with_comma(self.config_args_bias['bias_prior'])
            self.bias_prior        = np.tile(self.bias_prior[np.newaxis], (self.lens_ntomo, 1))
        
        self.m_shear_prior_parameters = np.array([np.zeros(self.source_ntomo), self.m_shear_prior_std]).T
        
        if(self.n_pcas_baryon > 0):
            baryon_priors = []
            for i in range(self.n_pcas_baryon):
                baryon_prior_i = self.config_args_baryons['prior_Q%d'%(i+1)].split(',')
                baryon_priors.append([float(baryon_prior_i[0]), float(baryon_prior_i[-1])])
            self.baryon_priors = np.array(baryon_priors)

    def get_starting_pos(self):
        p0 = []
        for x in self.params:
            if('prior' in self.params[x]):
                loc   = float(self.params[x]['ref']['loc'])
                scale = float(self.params[x]['ref']['scale'])
                p0_i = loc + scale * np.random.normal(size=self.n_walkers)
                p0.append(p0_i)        
        p0 = np.array(p0).T
        if self.probe!='cosmic_shear':
            bias_pars_std = 0.1 * (self.bias_prior[:,1] - self.bias_prior[:,0]) * np.ones(self.lens_ntomo)
            fast_pars_std = np.hstack([bias_pars_std, self.m_shear_prior_std])
        else:
            fast_pars_std = self.m_shear_prior_std
        if(self.n_pcas_baryon > 0):
            baryon_std = np.hstack([0.1 * np.ones(self.n_pcas_baryon)])
            fast_pars_std = np.hstack([fast_pars_std, baryon_std])
        p0_fast = fast_pars_std * np.random.normal(size=(self.n_walkers, self.n_fast_pars))
        if self.probe!='cosmic_shear':
            p0_fast[:,:self.lens_ntomo] = p0_fast[:,:self.lens_ntomo] + self.bias_fid
        p0 = np.hstack([p0, p0_fast])
        return p0
            
    def compute_datavector(self, theta):        
        if(self.emu_type=='nn'):
            theta = torch.Tensor(theta)
        elif(self.emu_type=='gp'):
            theta = theta[np.newaxis]
        datavector = self.emu.predict(theta)[0]        
        return datavector
    
    def add_bias(self, bias_theta, datavector):
        for i in range(self.lens_ntomo):
            factor = (bias_theta[i] / self.bias_fid[i])**self.galaxy_bias_mask[i]
            datavector = factor * datavector
        return datavector

    def add_baryon_q(self, Q, datavector):
        for i in range(self.n_pcas_baryon):
            datavector = datavector + Q[i] * self.baryon_pcas[:,i]
        return datavector

    def add_shear_calib(self, m, datavector):
        for i in range(self.source_ntomo):
            factor = (1 + m[i])**self.shear_calib_mask[i]
            datavector = factor * datavector
        return datavector

    def get_data_vector_emu(self, theta):
        theta_emu     = theta[:-self.n_fast_pars]
        datavector = self.compute_datavector(theta_emu)
        if self.probe!='cosmic_shear':
            bias_theta = theta[self.n_sample_dims-(self.n_pcas_baryon + self.source_ntomo + self.lens_ntomo):
                                  self.n_sample_dims-(self.n_pcas_baryon + self.source_ntomo)]
            datavector = self.add_bias(bias_theta, datavector)
        m_shear_theta = theta[self.n_sample_dims-(self.n_pcas_baryon + self.source_ntomo):
                              self.n_sample_dims-self.n_pcas_baryon]
        datavector = self.add_shear_calib(m_shear_theta, datavector)
        if(self.n_pcas_baryon > 0):
            baryon_q   = theta[-self.n_pcas_baryon:]
            datavector = self.add_baryon_q(baryon_q, datavector)
        return datavector

    def ln_prior(self, theta):
        flat_prior_theta     = theta[self.flat_prior_indices]
        gaussian_prior_theta = theta[self.gaussian_prior_indices]
        if self.probe!='cosmic_shear':
            bias_theta = theta[self.n_sample_dims-(self.n_pcas_baryon + self.source_ntomo + self.lens_ntomo):
                                  self.n_sample_dims-(self.n_pcas_baryon + self.source_ntomo)]
            prior_galaxy_bias = hard_prior(bias_theta, self.bias_prior)
        else:
            prior_galaxy_bias = 0.
        m_shear_theta        = theta[self.n_sample_dims-(self.n_pcas_baryon + self.source_ntomo):
                                     self.n_sample_dims-self.n_pcas_baryon]
        
        prior_flat    = hard_prior(flat_prior_theta, self.flat_prior_parameters)
        prior_gauss   = gaussian_prior(gaussian_prior_theta, self.gaussian_prior_parameters)
        prior_m_shear = gaussian_prior(m_shear_theta, self.m_shear_prior_parameters)
        if(self.n_pcas_baryon > 0):
            baryon_q   = theta[-self.n_pcas_baryon:]
            prior_baryons = hard_prior(flat_prior_theta, self.flat_prior_parameters)
        else:
            prior_baryons = 0.
        
        return prior_flat + prior_gauss + prior_galaxy_bias + prior_m_shear + prior_baryons
    
    def ln_lkl(self, theta):
        model_datavector = self.get_data_vector_emu(theta)
        delta_dv = (model_datavector - self.dv_obs)[self.mask]
        return -0.5 * delta_dv @ self.masked_inv_cov @ delta_dv        

    def ln_prob(self, theta, temper=1.):
        return self.ln_prior(theta) + temper * self.ln_lkl(theta)
