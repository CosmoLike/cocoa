import sys
from mpi4py import MPI
import numpy as np
import torch
from cocoa_emu import Config, get_lhs_params_list, CocoaModel
from cocoa_emu.emulator import NNEmulator, GPEmulator
from cocoa_emu.sampling import EmuSampler
import emcee

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

configfile = sys.argv[1]
config = Config(configfile)
    
# ============= LHS samples =================
from pyDOE import lhs

def get_lhs_samples(N_dim, N_lhs, lhs_minmax):
    unit_lhs_samples = lhs(N_dim, N_lhs)
    lhs_params = get_lhs_params_list(unit_lhs_samples, lhs_minmax)
    return lhs_params

if(rank==0):
    lhs_params = get_lhs_samples(config.n_dim, config.n_lhs, config.lhs_minmax)
else:
    lhs_params = None

lhs_params = comm.bcast(lhs_params, root=0)

# ================== Calculate data vectors ==========================

cocoa_model = CocoaModel(configfile, config.likelihood)

def get_local_data_vector_list(params_list, rank):
    train_params_list      = []
    train_data_vector_list = []
    N_samples = len(params_list)
    N_local   = N_samples // size    
    for i in range(rank * N_local, (rank + 1) * N_local):
        params_arr  = np.array(list(params_list[i].values()))
        data_vector = cocoa_model.calculate_data_vector(params_list[i])
        train_params_list.append(params_arr)
        train_data_vector_list.append(data_vector)
    return train_params_list, train_data_vector_list

def get_data_vectors(params_list, comm, rank):
    local_params_list, local_data_vector_list = get_local_data_vector_list(params_list, rank)
    if rank!=0:
        comm.send([local_params_list, local_data_vector_list], dest=0)
        train_params       = None
        train_data_vectors = None
    else:
        data_vector_list = local_data_vector_list
        params_list      = local_params_list
        for source in range(1,size):
            new_params_list, new_data_vector_list = comm.recv(source=source)
            data_vector_list = data_vector_list + new_data_vector_list
            params_list      = params_list + new_params_list
        train_params       = np.vstack(params_list)    
        train_data_vectors = np.vstack(data_vector_list)        
    return train_params, train_data_vectors

params_list = lhs_params
train_params, train_data_vectors = get_data_vectors(params_list, comm, rank)    

# ================== Train emulator ==========================

if(rank==0):
    print("Done computing the data vectors...")
    np.save(config.savedir + '/train_data_vectors.npy', train_data_vectors)
    np.save(config.savedir + '/train_params.npy', train_params)
    if(config.emu_type=='nn'):
        emu = NNEmulator(config.n_dim, config.output_dims, config.dv_fid, config.dv_std)
        print("Training emulator...")
        emu.train(torch.Tensor(train_params), torch.Tensor(train_data_vectors), 
                  batch_size=config.batch_size, n_epochs=config.n_epochs)
    elif(config.emu_type=='gp'):
        emu = GPEmulator(config.n_dim, config.output_dims, config.dv_fid, config.dv_std)
        emu.train(train_params, train_data_vectors)
# ============= Sample from the posterior ======================
if(rank==0):
    emu_sampler = EmuSampler(emu, config)
    
    pos0 = emu_sampler.get_starting_pos()

#     print("pos0.shape: "+str(pos0.shape))
#     print("config.n_emcee_walkers: "+str(config.n_emcee_walkers))
#     print("emu_sampler.n_sample_dims: "+str(emu_sampler.n_sample_dims))
    
    sampler = emcee.EnsembleSampler(config.n_emcee_walkers, emu_sampler.n_sample_dims, emu_sampler.ln_prob)
    sampler.run_mcmc(pos0, config.n_mcmc, progress=True)
    
MPI.Finalize