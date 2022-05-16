import sys
from mpi4py import MPI
import numpy as np
from cocoa_emu import Config, get_lhs_params_list, CocoaModel

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

configfile = sys.argv[1]
config = Config(configfile)
    
# ============= LHS samples =================
from pyDOE import lhs

def get_lhs_samples(N_dim, N_lhs, priors):
    unit_lhs_samples = lhs(N_dim, N_lhs)
    lhs_params = get_lhs_params_list(unit_lhs_samples, priors)
    return lhs_params

if(rank==0):
    lhs_params = get_lhs_samples(config.N_dim, config.N_lhs, config.priors)
else:
    lhs_params = None

lhs_params = comm.bcast(lhs_params, root=0)
# ============================================

cocoa_model = CocoaModel(configfile, config.likelihood)

def get_local_data_vector_list(params_list, rank):
    train_params_list      = []
    train_data_vector_list = []
    N_samples = len(params_list)
    N_local   = N_samples // size    
    for i in range(rank * N_local, (rank + 1) * N_local):
        print("i: %d"%(i))
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
            print("source: %d"%(source))
            new_params_list, new_data_vector_list = comm.recv(source=source)
            data_vector_list = data_vector_list + new_data_vector_list
            params_list      = params_list + new_params_list
        train_params       = np.vstack(params_list)    
        train_data_vectors = np.vstack(data_vector_list)        
    return train_params, train_data_vectors

params_list = lhs_params
train_params, train_data_vectors = get_data_vectors(params_list, comm, rank)    
if(rank==0):
    print("Done computing the data vectors...")
    np.save(config.savedir + '/train_data_vectors.npy', train_data_vectors)
    np.save(config.savedir + '/train_params.npy', train_params)

MPI.Finalize
# ============================================
