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
    
class EmuSampler:
    def __init__(self, emu, config):
        self.emu            = emu
        self.params         = config.params
        
        self.get_priors()
        
        self.n_walkers      = config.n_emcee_walkers
        self.n_pcas_baryon  = config.n_pcas_baryon
        self.baryon_pcas    = config.baryon_pcas
        
        self.emu_type       = config.emu_type
        self.mask           = config.mask
        self.cov            = config.cov
        self.masked_inv_cov = np.linalg.inv(self.cov[self.mask][:,self.mask])
        self.dv_obs         = config.dv_obs
        
        self.n_sample_dims = config.n_dim + config.n_pcas_baryon

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

    def get_starting_pos(self):
        p0 = []
        for x in self.params:
            if('prior' in self.params[x]):
                loc   = float(self.params[x]['ref']['loc'])
                scale = float(self.params[x]['ref']['scale'])
                p0_i = loc + scale * np.random.normal(size=self.n_walkers)
                p0.append(p0_i)        
        p0 = np.array(p0).T
        if self.n_pcas_baryon!=0:
            p0_baryons = 0.1 * np.random.normal(size=(self.n_walkers, self.n_pcas_baryon))
            p0 = np.hstack([p0, p0_baryons])
        return p0

    def ln_prior(self, theta):
        flat_prior_theta     = theta[self.flat_prior_indices]
        gaussian_prior_theta = theta[self.gaussian_prior_indices]

        prior_flat  = hard_prior(flat_prior_theta, self.flat_prior_parameters)
        prior_gauss = gaussian_prior(gaussian_prior_theta, self.gaussian_prior_parameters)
        
        return prior_flat + prior_gauss
            
    def compute_datavector(self, theta):
        if(self.emu_type=='nn'):
            theta = torch.Tensor(theta)
        elif(self.emu_type=='gp'):
            theta = theta[np.newaxis]
        datavector = self.emu.predict(theta)[0]        
        return datavector
    
    def add_baryon_q(self, Q, datavector):
        for i in range(self.n_pcas_baryon):
            datavector = datavector + Q[i] * self.baryon_pcas[:,i]
        return datavector

    def get_data_vector_emu(self, theta):
        theta_emu  = theta[:-self.n_pcas_baryon]
        baryon_q   = theta[-self.n_pcas_baryon:]
        datavector = self.compute_datavector(theta_emu)
        datavector = self.add_baryon_q(baryon_q, datavector)
        return datavector

    def ln_lkl(self, theta):
        model_datavector = self.get_data_vector_emu(theta)
        delta_dv = (model_datavector - self.dv_obs)[self.mask]
        return -0.5 * delta_dv @ self.masked_inv_cov @ delta_dv        

    def ln_prob(self, theta, temper=1.):
        return self.ln_prior(theta) + temper * self.ln_lkl(theta)