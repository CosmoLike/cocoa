import sys
import numpy as np
from cocoa_emu import Config, CocoaModel
from cocoa_emu.emulator import NNEmulator, GPEmulator
from cocoa_emu.sampling import EmuSampler
import emcee

from multiprocessing import Pool

def ln_prob(theta, temper_val):
     return emu_sampler.ln_prob(theta, temper_val)
    
configfile = sys.argv[1]
config = Config(configfile)

try:
    emu_iter = int(sys.argv[2])
except:
    emu_iter = config.n_train_iter - 1

cocoa_model = CocoaModel(configfile, config.likelihood)
emu = NNEmulator(config.n_dim, config.output_dims, config.dv_fid, config.dv_std)
emu.load(config.savedir + '/model_%d'%(emu_iter))

emu_sampler = EmuSampler(emu, config)
pos0 = emu_sampler.get_starting_pos()

import os
os.environ["OMP_NUM_THREADS"] = "1"

with Pool() as pool:
    sampler = emcee.EnsembleSampler(config.n_emcee_walkers, emu_sampler.n_sample_dims, 
                                        ln_prob, args=(1.,), pool=pool)
    sampler.run_mcmc(pos0, config.n_mcmc, progress=True)
    
samples = sampler.chain[:,config.n_burn_in::config.n_thin].reshape((-1, emu_sampler.n_sample_dims))    

np.save(config.savedir + '/emu_chain_%d.npy'%(emu_iter), samples)
