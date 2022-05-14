import sys
from cocoa_emu import Config, get_lhs_params_list, CocoaModel

configfile = sys.argv[1]
config = Config(configfile)
    
# ============= LHS samples =================
from pyDOE import lhs

def get_lhs_samples(N_dim, N_lhs, priors):
    unit_lhs_samples = lhs(N_dim, N_lhs)
    lhs_params = get_lhs_params_list(unit_lhs_samples, priors)
    return lhs_params

lhs_params = get_lhs_samples(config.N_dim, config.N_lhs, config.priors)
# ============================================

cocoa_model = CocoaModel(configfile, config.likelihood)

def get_data_vector_list(params_list):
    train_data_vector_list = []
    for i in range(0, len(params_list)):
        print("i: %d"%(i))
        data_vector = cocoa_model.calculate_data_vector(params_list[i])
        train_data_vector_list.append(data_vector)
    return train_data_vector_list

get_data_vector_list(lhs_params)
# ============================================