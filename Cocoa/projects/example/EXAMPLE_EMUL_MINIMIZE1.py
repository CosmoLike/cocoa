import functools, ipyparallel, scipy, iminuit, copy, argparse, random
import emcee, itertools
import numpy as np
from cobaya.theories.emulcmb.emulcmb2 import emulcmb
from cobaya.theories.emultheta.emultheta2 import emultheta
from getdist import IniFile
from mflike import TTTEEE, BandpowerForeground
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# initial guess
x = np.array([
    3.045845885,         # logA
    9.652308970e-01,     # ns
    1.0410562,           # thetastar
    2.246801442e-02,     # omegabh2
    1.198257361e-01 ,    # omegach2
    5.433339482e-02,     # tau
    3.299661491,         # a_tSZ
    1.603283271,         # a_kSZ
    6.909072608,         # a_p
    2.081592320,         # beta_p
    4.903593522,         # a_c
    2.201183246,         # beta_c
    3.100233500,         # a_s
    2.843185734,         # a_gtt
    1.910150877e-01 ,    # xi
    1.028875217e+01,     # T_d
    0.10,                # a_gte
    5.316779467e-04,     # a_pste
    9.891127676e-02,     # a_gee
    0.0,          # a_psee
], dtype='float64')

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

MKL = TTTEEE(
    {
        "data_folder": "simons_observatory/v0.8",
        "input_file": "LAT_simu_sacc_00044.fits",
        "cov_Bbl_file": "data_sacc_w_covar_and_Bbl.fits",
        "defaults": {
            "polarizations": ['TT', 'TE', 'ET', 'EE'],
            "scales": {
                "TT": [50, 3002],
                "TE": [50, 3002],
                "ET": [50, 3002],
                "EE": [50, 3002],
            },
            "symmetrize": False,
        },
        "lmax_theory": 5000
    }
)
fg  = BandpowerForeground(MKL.get_fg_requirements())

ecmb = emulcmb(extra_args={  
    'eval': [True, True, True, False],
    'file': ['external_modules/data/emultrf/CMB_TRF/chiTTAstautrf1dot2milnewlhcevansqrtrescalec16.pt',
             'external_modules/data/emultrf/CMB_TRF/chiTEAstautrf1dot2milnewlhcevansqrtrescalec16.pt',
             'external_modules/data/emultrf/CMB_TRF/chiEEAstautrf1dot2milnewlhcevansqrtrescalec16.pt',
             None],
    'extra': ['external_modules/data/emultrf/CMB_TRF/extrainfo_lhs_tt_96.npy',
              'external_modules/data/emultrf/CMB_TRF/extrainfo_lhs_te_96.npy',
              'external_modules/data/emultrf/CMB_TRF/extrainfo_lhs_ee_96.npy',
              None],
    'ord': [['omegabh2','omegach2','H0','tau','ns','logA','mnu','w','wa'],
            ['omegabh2','omegach2','H0','tau','ns','logA','mnu','w','wa'],
            ['omegabh2','omegach2','H0','tau','ns','logA','mnu','w','wa'], 
            None],
    'extrapar': [{'ellmax' : 5000, 'MLA': 'TRF', 'NCTRF': 16, 'INTDIM': 4, 'INTTRF': 5120},
                 {'ellmax' : 5000, 'MLA': 'TRF', 'NCTRF': 16, 'INTDIM': 4, 'INTTRF': 5120},
                 {'ellmax' : 5000, 'MLA': 'TRF', 'NCTRF': 16, 'INTDIM': 4, 'INTTRF': 5120},
                 None]})

etheta = emultheta(extra_args={ 
    'file': ['external_modules/data/emultrf/CMB_TRF/thetaH0GP.joblib'],
    'extra':['external_modules/data/emultrf/CMB_TRF/extrainfotheta.npy'],
    'ord':  [['omegabh2','omegach2','thetastar']]})

def chi2(p):
    NP = {
        "T_effd": 19.6,
        "beta_d": 1.5,
        "beta_s": -2.5,
        "alpha_s": 1,
        "bandint_shift_LAT_93": 0,
        "bandint_shift_LAT_145": 0,
        "bandint_shift_LAT_225": 0,
        "cal_LAT_93": 1,
        "cal_LAT_145": 1,
        "cal_LAT_225": 1,
        "calG_all": 1,
        "alpha_LAT_93": 0,
        "alpha_LAT_145": 0,
        "alpha_LAT_225": 0,
        "a_tSZ": p[6],
        "a_kSZ": p[7],
        "a_p": p[8],
        "beta_p": p[9],
        "a_c": p[10],
        "beta_c": p[11],
        "a_s": p[12],
        "T_d": p[15],
        "a_gtt": p[13],
        "xi": p[14],
        "alpha_dT": -0.6,
        "alpha_p": 1.0,
        "alpha_tSZ": 0.0,
        "calT_LAT_93": 1.0,
        "calT_LAT_145": 1.0,
        "calT_LAT_225": 1.0,
        "a_gte": p[16],
        "a_pste": p[17],
        "alpha_dE": -0.4,
        "a_gee": p[18],
        "a_psee": p[19],
        "alpha_dE": -0.4,
        "calE_LAT_93": 1,
        "calE_LAT_145": 1,
        "calE_LAT_225": 1
    }
    param = {'logA': p[0],
             'ns': p[1],
             'thetastar': p[2],
             'omegabh2': p[3],
             'omegach2': p[4],
             'tau': p[5],
             'mnu': 0.06,
             'w': -1.0,
             'wa': 0.0,
             'omegamh2': p[4]+p[3]+(0.06*(3.046/3)**0.75)/94.0708
             }
    param = param | etheta.get_params(param)
    cl    = ecmb.get_Cl(params=param, ell_factor=True)
    fmt   = fg.get_foreground_model_totals(**NP)
    return -2*MKL.loglike(cl, fmt, **NP)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

parser = argparse.ArgumentParser(prog='EXAMPLE_EMUL_MINIMIZE1')

parser.add_argument("--tol",
                    dest="tolerance",
                    help="Minimizer Tolerance",
                    type=float,
                    nargs='?',
                    const=1,
                    default=0.05)

parser.add_argument("--maxfeval",
                    dest="maxfeval",
                    help="Minimizer: maximum number of likelihood evaluations",
                    type=int,
                    nargs='?',
                    const=1,
                    default=5000)

parser.add_argument("--maxiter",
                    dest="maxiter",
                    help="Minimizer: maximum number of minimizer iterations",
                    type=int,
                    nargs='?',
                    const=1,
                    default=3)

parser.add_argument("--minmethod",
                    dest="minmethod",
                    help="Minimizer: minimizer method",
                    type=int,
                    nargs='?',
                    const=1,
                    default=5)

parser.add_argument("--root",
                    dest="root",
                    help="Name of the Output File",
                    nargs='?',
                    const=1,
                    default="./projects/example/")

parser.add_argument("--outroot",
                    dest="outroot",
                    help="Name of the Output File",
                    nargs='?',
                    const=1,
                    default="test.dat")

parser.add_argument("--numpts",
                    dest="numpts",
                    help="Number of Points to Compute Minimum",
                    type=int,
                    nargs='?',
                    const=1,
                    default=20)

parser.add_argument("--factor",
                    dest="factor",
                    help="Factor that set the bounds (multiple of cov matrix)",
                    type=int,
                    nargs='?',
                    const=1,
                    default=3)

parser.add_argument("--ref",
                    dest="ref",
                    help="Refinement Factor (min method 2)",
                    type=int,
                    nargs='?',
                    const=1,
                    default=1) # zero or one

# need to use parse_known_args because of mpifuture 
args, unknown = parser.parse_known_args()
tol         = args.tolerance
maxfeval    = args.maxfeval
maxiter     = args.maxiter
min_method  = args.minmethod
oroot       = args.root + "chains/" + args.outroot
numpts      = args.numpts
ref         = args.ref

cov_file = args.root + 'EXAMPLE_EMUL_MCMC1.covmat'
cov      = np.loadtxt(cov_file)[0:len(x),0:len(x)]
sigma    = np.sqrt(np.diag(cov))
bounds   = np.c_[x - args.factor*sigma, x + args.factor*sigma]

class GaussianStep:
    def __init__(self, stepsize=0.2):
       self.cov = stepsize*cov
    def __call__(self, x):
        print(len(x), len(self.cov))
        return np.random.multivariate_normal(x, self.cov, size=1)

def min_chi2(x0, 
             bounds, 
             min_method, 
             tol=0.01, 
             maxfeval=3000, 
             maxiter=10, 
             cov=cov):

    def mychi2(params, *args):
        T = args
        return chi2(p=params)/T

    args = (1.0)

    def log_prior(params):
        return 1.0
    
    def logprob(params, *args):
        lp = log_prior(params)
        if not np.isfinite(lp):
            return -np.inf
        else:
            return -0.5*mychi2(params, *args) + lp

    if min_method == 1:    

        x = copy.deepcopy(x0)
        partial_samples = [x]
        partial = [mychi2(x, *args)]
        
        for i in range(maxiter):
            tmp = iminuit.minimize(fun=mychi2, 
                                   x0=x, 
                                   args=args, 
                                   bounds=bounds, 
                                   method="migrad", 
                                   tol=tol,
                                   options = {'stra'  : 2, 
                                              'maxfun': maxfeval})
            partial_samples.append(tmp.x)
            partial.append(tmp.fun)
            j = np.argmin(np.array(partial))
            x = GaussianStep(stepsize=0.005)(partial_samples[j])[0,:]
            print(f"MN: i = {i}, chi2 = {partial[j]}, params = {partial_samples[j]}")
        
        j = np.argmin(np.array(partial))
        result = [partial_samples[j], partial[j]]
        print(f"All (MN): chi2 = {partial[j]}, params = {partial_samples[j]}")

    elif min_method == 2: # adapted from PROCOLI       
    
        ndim        = int(x0.shape[0])
        nwalkers    = int(2*x0.shape[0])
        nsteps      = maxfeval
        temperature = np.array([1.0, 0.25, 0.1, 0.005, 0.001], dtype='float64')
        stepsz      = temperature/4.0
        partial_samples = []
        partial = []

        for i in range(5):
            x = [] # Initial point
            for j in range(nwalkers):
                x.append(GaussianStep(stepsize=stepsz[i])(x0)[0,:])
            x = np.array(x,dtype='float64')

            GScov  = copy.deepcopy(cov)
            GScov *= temperature[i]*stepsz[i] 
      
            sampler = emcee.EnsembleSampler(nwalkers, 
                                            ndim, 
                                            logprob, 
                                            args=(args[0], args[1], temperature[i]),
                                            moves=[(emcee.moves.GaussianMove(cov=GScov),1.)])
            
            sampler.run_mcmc(x, nsteps, skip_initial_state_check=False)
            samples = sampler.get_chain(flat=True, thin=1, discard=0)

            j = np.argmin(-1.0*np.array(sampler.get_log_prob(flat=True)))
            
            partial_samples.append(samples[j])
            partial.append(mychi2(samples[j], *args))
            j = np.argmin(np.array(partial))
            x0 = copy.deepcopy(partial_samples[j])
            sampler.reset()
            print(f"emcee: i = {i}, chi2 = {partial[j]}, params = {partial_samples[j]}")
        
        j = np.argmin(np.array(partial))
        result = [partial_samples[j], partial[j]]
        print(f"All (emcee): chi2 = {partial[j]}, params = {partial_samples[j]}")
    
    return result

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

def minimizer(x0, index, min_method, maxiter, maxfeval):
    res =  min_chi2(x0=np.array(x0[index,:], dtype='float64'), 
                    bounds=np.array(bounds, dtype='float64'), 
                    min_method=min_method,
                    tol=tol, 
                    maxfeval=maxfeval, 
                    maxiter=maxiter)
    return res

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

from mpi4py.futures import MPIPoolExecutor

if __name__ == '__main__':
    
    print(f"min_method={min_method}, maxiter={maxiter}, maxfeval={maxfeval}, tol={tol}")

    executor = MPIPoolExecutor()
    
    x0 = np.zeros((numpts,len(x)), dtype='float64')
    for i in range(numpts):
        x0[i,:] = GaussianStep(stepsize=1.0)(x0[:,i])[0,:]
    
    res = np.array(list(executor.map(functools.partial(prf,
                                                       index=index, 
                                                       min_method=min_method,
                                                       maxiter=maxiter,
                                                       maxfeval=maxfeval),
                                     x0)))

    j = np.argmin(np.array(res[:,1]))

    rnd = random.randint(0,1000)
    out = oroot + "_" + str(rnd) + "_method_" + str(min_method)
    print("Output file = ", out +".txt")
    print(res[j,1])
    print(res[j,0])
    np.savetxt(out + ".txt", np.c_[res[j,1], res[j,0]]) # + here = concatenate

    if ref > 0:
        for i in range(numpts):
            print(res[i,0])
            x0[i,:] = GaussianStep(stepsize=1.0)(res[i,0])[0,:]

        res = np.array(list(executor.map(functools.partial(prf,
                                                           index=index, 
                                                           min_method=1,
                                                           maxiter=maxiter,
                                                           maxfeval=maxfeval/10),
                                         x0)))
        print("Output file = ", out + "_ref" + ".txt")
        np.savetxt(out + "_ref" + ".txt", np.c_[res[j,1], res[j,0]])
    executor.shutdown()
