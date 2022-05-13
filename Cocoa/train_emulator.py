import sys
from cocoa_emu import Config, get_lhs_params_list

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
