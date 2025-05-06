import sys, platform, os
from cobaya.theories.emulcmbtrf.emuprofile import EmulProfile
import functools
import numpy as np
import ipyparallel
import sys, platform, os
import math
import scipy
from getdist import IniFile
import itertools
import iminuit
import copy
import argparse
import random
import emcee
import itertools
from mflike import TTTEEE, BandpowerForeground

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

x = np.array([
                3.04,         # logA
                0.96,         # ns
                67.0,         # H0
                0.02,         # omegabh2
                .12,          # omegach2
                0.05,         # tau
                3.29,         # a_tSZ
                1.60,         # a_kSZ
                6.90,         # a_p
                2.08,         # beta_p
                4.90,         # a_c
                2.20,         # beta_c
                3.10,         # a_s
                2.84,         # a_gtt
                1.9e-01,      # xi
                10.2,         # T_d
                0.10,         # a_gte
                5.31e-04,     # a_pste
                9.8e-02,      # a_gee
                0.0,          # a_psee
              ], dtype='float64')

bounds = np.array([
                    [1.61, 3.91],   # logAs
                    [0.92, 1.05],   # ns 
                    [60.0, 90.0],   # H0
                    [0.01, 0.04],   # omegabh2
                    [0.06, 0.20],   # omegach2
                    [0.02, 0.09],   # tau
                    [3.0, 3.6],     # a_tSZ
                    [1.4, 1.8],     # a_kSZ
                    [6.2, 7.6],     # a_p
                    [1.8, 2.2],     # beta_p
                    [4.4, 5.4],     # a_c
                    [2.0, 2.4],     # beta_c
                    [2.8, 3.4],     # a_s
                    [1.5, 4.0],     # a_gtt
                    [0.0, 0.2],     # xi
                    [8.60, 10.60],  # T_d
                    [0.20, 0.50],   # a_gte
                    [-1.00, 1.00],  # a_pste
                    [-0.05, 0.30],  # a_gee
                    [0.00, 0.5],    # a_psee
                  ], dtype='float64')

start = np.array([ 
                     2.974,      # logAs
                     0.9439,     # ns 
                     64.66,      # H0
                     0.0216,     # omegabh2
                     0.114,      # omegach2
                     0.0252,     # tau
                     3.0,        # a_tSZ
                     1.4,        # a_kSZ
                     6.2,        # a_p
                     1.8,        # beta_p                  
                     4.4,        # a_c
                     2.0,        # beta_c
                     2.8,        # a_s
                     0.54,       # a_gtt
                     0.0,        # xi
                     8.60,       # T_d
                     0.16,       # a_gte
                    -0.5,        # a_pste
                    -0.0,        # a_gee
                    -0.20        # a_psee
                ], dtype='float64')

stop  = np.array([ 
                    3.114,      # logAs
                    0.9859,     # ns
                    70.06,      # H0
                    0.02312,    # omegabh2
                    0.126,      # omegach2
                    0.0836,     # tau
                    3.6,        # a_tSZ
                    1.8,        # a_kSZ
                    7.6,        # a_p
                    2.2,        # beta_p                  
                    5.4,        # a_c
                    2.4,        # beta_c
                    3.4,        # a_s
                    5.04,       # a_gtt
                    0.2,        # xi
                    10.60,      # T_d
                    0.56,       # a_gte
                    -0.5,       # a_pste
                    0.3,        # a_gee
                    0.20        # a_psee
                 ], dtype='float64')

name  = [ 
            "logAs",     # As
            "ns",        # ns
            "H0",        # H0
            "omegabh2",  # omegabh2
            "omegach2",  # omegach2
            "tau",       # tau
            "a_tSZ",     # a_tSZ
            "a_kSZ",     # a_kSZ
            "a_p",       # a_p
            "beta_p",    # beta_p          
            "a_c",       # a_c
            "beta_c",    # beta_c
            "a_s",       # a_s
            "a_gtt",     # a_gtt
            "xi"         # xi
            "T_d"        # T_d
            "a_gte",     # a_gte
            "a_pste"     # a_pste
            "a_gee",     # a_gee
            "a_psee"     # a_psee
        ]
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

MKL = TTTEEE(
    {
        "data_folder": "simons_observatory/v0.8",
        #"packages_path": packages_path,
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

tmp='external_modules/data/emultrf/CMB_TRF/'
emu=EmulProfile(ROOTDIR=os.environ.get("ROOTDIR") , 
                ttfilename=tmp+'/chiTTAstautrf1dot2milnewlhcevansqrtrescalec16', 
                tefilename=tmp+'/chiTEAstautrf1dot2milnewlhcevansqrtrescalec16', 
                eefilename=tmp+'/chiEEAstautrf1dot2milnewlhcevansqrtrescalec16', 
                ttextraname=tmp+'/extrainfo_lhs_tt_96.npy', 
                teextraname=tmp+'/extrainfo_lhs_te_96.npy', 
                eeextraname=tmp+'/extrainfo_lhs_ee_96.npy')

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
    } | {
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
        "alpha_p": 1,
        "alpha_tSZ": 0.,
        "calT_LAT_93": 1,
        "calT_LAT_145": 1,
        "calT_LAT_225": 1,
    } | {
        "a_gte": p[16],
        "a_pste": p[17],
        "alpha_dE": -0.4
    } | {
        "a_gee": p[18],
        "a_psee": p[19],
        "alpha_dE": -0.4,
        "calE_LAT_93": 1,
        "calE_LAT_145": 1,
        "calE_LAT_225": 1
    }
    #[ombh2, omch2, H0, tau, ns, logAs]
    return -2 * (MKL.loglike(emu.get_Cl(
        cmb_params=[p[3],p[4],p[2],p[5],p[1],p[0]]),
        fg.get_foreground_model_totals(**NP),**NP)-MKL.logp_const)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

parser = argparse.ArgumentParser(prog='EXAMPLE_PROFILE1')

parser.add_argument("--tol",
                    dest="tolerance",
                    help="Minimizer Tolerance",
                    type=float,
                    nargs='?',
                    const=1,
                    default=0.02)

parser.add_argument("--maxfeval",
                    dest="maxfeval",
                    help="Minimizer: maximum number of likelihood evaluations",
                    type=int,
                    nargs='?',
                    const=1,
                    default=75000)

parser.add_argument("--maxiter",
                    dest="maxiter",
                    help="Minimizer: maximum number of minimizer iterations",
                    type=int,
                    nargs='?',
                    const=1,
                    default=5)

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

parser.add_argument("--profile",
                    dest="profile",
                    help="Which Parameter to Profile",
                    type=int,
                    nargs='?',
                    const=1,
                    default=1)

parser.add_argument("--mpi",
                    dest="mpi",
                    help="Number of MPI cores",
                    type=int,
                    nargs='?',
                    const=1)

# need to use parse_known_args because of mpifuture 
args, unknown = parser.parse_known_args()

tol             = args.tolerance
maxfeval        = args.maxfeval
maxiter         = args.maxiter
min_method      = args.minmethod
oroot           = args.root + "chains/" + args.outroot
index           = args.profile
nummpi          = args.mpi

cov_file = args.root + 'EXAMPLE_EVALUATE20.covmat'
cov      = np.loadtxt(cov_file)

print(cov.shape)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

def min_chi2(x0, 
             bounds, 
             min_method, 
             fixed=-1, 
             tol=0.01, 
             maxfeval=3000, 
             maxiter=10, 
             cov=cov):

    def mychi2(params, *args):
        z, fixed, T = args
        params = np.array(params)
        if fixed > -1:
            params = np.insert(params, fixed, z)
        res = chi2(params=params)/T
        return res

    if fixed > -1:
        z      = x0[fixed]
        x0     = np.delete(x0, (fixed))
        bounds = np.delete(bounds, (fixed), axis=0)
        args = (z, fixed, 1.0)
    else:
        args = (0.0, -2.0, 1.0)

    def log_prior(params):
        params = np.array(copy.deepcopy(params))
        prior = 1.0
        for i in range(len(params)):
            if (params[i] < bounds[i][0]) or (params[i] > bounds[i][1]):
                return -np.inf
            else:
                prior *= np.log(1.0/(bounds[i][1]-bounds[i][0]))
        return prior
    
    def logprob(params, *args):
        lp = log_prior(params)
        if not np.isfinite(lp):
            return -np.inf
        else:
            return -0.5*mychi2(params, *args) + lp
    
    cov = np.delete(cov, (fixed), axis=0)
    cov = np.delete(cov, (fixed), axis=1)

    class GaussianStep:
       def __init__(self, stepsize=0.2):
           self.cov = stepsize*cov
       def __call__(self, x):
           print(self.cov.shape)
           print(x.shape) 
           return np.random.multivariate_normal(x, self.cov, size=1)

    if min_method == 1:
        x = copy.deepcopy(x0)
        tmp = iminuit.minimize(fun=mychi2, 
                               x0=x, 
                               args=args, 
                               bounds=bounds, 
                               method="migrad", 
                               tol=tol,
                               options = {'stra'  : 2, 
                                          'maxfun': int(maxfeval/3.0)}) 
        
        x = copy.deepcopy(tmp.x)  
        GS = GaussianStep()
        tmp = scipy.optimize.basinhopping(func=mychi2, 
                                          x0=x, 
                                          T=0.1, 
                                          target_accept_rate=0.25, 
                                          niter=maxiter, 
                                          stepsize=0.02,
                                          interval=100,
                                          niter_success=3,
                                          take_step=GS,
                                          minimizer_kwargs={"method": 'Nelder-Mead', 
                                                            "args": args, 
                                                            "bounds": bounds, 
                                                            "options": {'adaptive' : True, 
                                                                        'xatol'    : tol,
                                                                        'fatol'    : tol, 
                                                                        'maxfev'   : maxfeval,
                                                                        'maxiter'  : maxiter}})
        x = copy.deepcopy(tmp.x)  
        tmp = iminuit.minimize(fun=mychi2, 
                               x0=x, 
                               args=args, 
                               bounds=bounds, 
                               method="migrad", 
                               tol=tol,
                               options = {'stra'  : 2, 
                                          'maxfun': int(maxfeval/3.0)})
        result = [tmp.x, tmp.fun]
    
    elif min_method == 2:
        x = copy.deepcopy(x0)
        tmp = iminuit.minimize(fun=mychi2, 
                               x0=x, 
                               args=args, 
                               bounds=bounds, 
                               method="migrad", 
                               tol=tol,
                               options = {'stra'  : 2, 
                                          'maxfun': int(maxfeval/3.0)})
        
        x = copy.deepcopy(tmp.x)
        # https://stats.stackexchange.com/a/456073
        tmp = scipy.optimize.dual_annealing(func=mychi2, 
                                            x0=x, 
                                            args=args, 
                                            bounds=bounds, 
                                            maxfun=maxfeval,
                                            no_local_search=True, 
                                            maxiter=maxiter, 
                                            visit=1.01, 
                                            accept=1, 
                                            initial_temp=5230.0, 
                                            restart_temp_ratio=0.0002)
        x = copy.deepcopy(tmp.x)  
        tmp = iminuit.minimize(fun=mychi2, 
                               x0=x, 
                               args=args, 
                               bounds=bounds, 
                               method="migrad", 
                               tol=tol,
                               options = {'stra'  : 2, 
                                          'maxfun': int(maxfeval/3.0)})
        result = [tmp.x, tmp.fun]
    
    elif min_method == 3:
        tmp = scipy.optimize.shgo(func=mychi2, 
                                  args=args, 
                                  bounds=bounds,
                                  n=50,
                                  iters=maxiter,
                                  options={'f_tol'    : tol, 
                                           'maxfev'   : maxfeval,
                                           'maxiter'  : maxiter},
                                  minimizer_kwargs={"method": 'Nelder-Mead', 
                                                              "args": args, 
                                                              "bounds": bounds, 
                                                              "options": {'adaptive' : True, 
                                                                          'xatol'    : tol,
                                                                          'fatol'    : tol, 
                                                                          'maxfev'   : maxfeval,
                                                                          'maxiter'  : maxiter}})
        result = [tmp.x, tmp.fun]
    elif min_method == 4:
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
            print(f"MN: i = {i}, chi2 = {partial[j]}, param = {args[0]}")
        j = np.argmin(np.array(partial))
        result = [partial_samples[j], partial[j]]   
    
    elif min_method == 5: # adapted from PROCOLI       
        ndim        = int(x0.shape[0])
        nwalkers    = int(2*x0.shape[0])
        nsteps      = maxfeval
        #temperature = np.array([1.0, 0.3, 0.25, 0.2, 0.1, 0.005, 0.001], dtype='float64')
        temperature = np.array([1.0, 0.25, 0.1, 0.005, 0.001], dtype='float64')
        stepsz      = temperature/4.0

        partial_samples = []
        partial = []

        for i in range(5):
            print(f"emcee: i = {i} begins")

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
            sampler.run_mcmc(x, nsteps, skip_initial_state_check=True)
            
            samples = sampler.get_chain(flat=True, thin=1, discard=0)
            j = np.argmin(-1.0*np.array(sampler.get_log_prob(flat=True)))
            
            partial_samples.append(samples[j])
            partial.append(mychi2(samples[j], *args))
            
            j = np.argmin(np.array(partial))
            x0 = copy.deepcopy(partial_samples[j])
            
            sampler.reset()

            print(f"emcee: i = {i}, chi2 = {partial[j]}, param = {args[0]}")
        
        tmp = iminuit.minimize(fun=mychi2, 
                               x0=x0, 
                               args=args, 
                               bounds=bounds, 
                               method="migrad", 
                               tol=tol,
                               options = {'stra'  : 2, 
                                          'maxfun': int(maxfeval/3)})
        
        partial_samples.append(tmp.x)
        partial.append(tmp.fun)  
        
        j = np.argmin(np.array(partial))
        result = [partial_samples[j], partial[j]]
        print(f"MN: chi2 = {result}, param = {args[0]}")
        print(f"All: chi2 = {partial}, param = {args[0]}")
    return result

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

def prf(x0, index, min_method, maxiter, maxfeval):
    res =  min_chi2(x0=np.array(x0, dtype='float64'), 
                    bounds=np.array(bounds, dtype='float64'), 
                    min_method=min_method, 
                    fixed=index, 
                    tol=tol, 
                    maxfeval=maxfeval, 
                    maxiter=maxiter)
    return res


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

from mpi4py.futures import MPIPoolExecutor

if __name__ == '__main__':
    
    print(f"min_method={min_method}, maxiter={maxiter}, maxfeval={maxfeval}, tol={tol}, param={index}")

#    executor = MPIPoolExecutor()
    
    param = np.linspace(start=start[index], stop=stop[index], num=nummpi)
    print(f"profile param values = {param}")

    x0 = []
    for i in range(nummpi):
        tmp = copy.deepcopy(x)
        tmp[index] = param[i]
        x0.append(tmp)
    x0 = np.array(x0, dtype=object)

#    res = np.array(list(executor.map(functools.partial(prf,
    res = np.array(list(map(functools.partial(prf,
                                                       index=index, 
                                                       min_method=min_method,
                                                       maxiter=maxiter,
                                                       maxfeval=maxfeval),
                                     x0)))

    rnd = random.randint(0,1000)
    out = oroot + "_" + str(rnd) + "_method_" + str(min_method) + "_" + name[index] + ".txt"
    print("Output file = ", out)
    np.savetxt(out, np.c_[param, res[:,1]])

    if min_method != 4:
        x0 = []
        for i in range(nummpi):
            x0.append(np.insert(res[i,0], index, param[i]))
        x0 = np.array(x0, dtype=object)
        res = np.array(list(executor.map(functools.partial(prf,
                                                           index=index, 
                                                           min_method=4,
                                                           maxiter=maxiter,
                                                           maxfeval=int(2500)),
                                         x0)))
        rnd = random.randint(0,1000)
        out = oroot + "_" + str(rnd) + "_method_" + str(4) + "_" + name[index] + ".txt"
        print("Output file = ", out)
        np.savetxt(out, np.c_[param, res[:,1]])

#    executor.shutdown()


#mpirun -n 10 python -m mpi4py.futures ./projects/example/EXAMPLE_PROFILE1.py --tol 0.02 --profile 1 --maxiter 5 --maxfeval 100 --mpi 10 --outroot "test1" --minmethod 5 

#mpirun -n ${SLURM_NTASKS} --oversubscribe --mca btl vader,tcp,self \
#  --bind-to core:overload-allowed --map-by numa:pe=${OMP_NUM_THREADS} \
#  python -m mpi4py.futures EXAMPLE_PROFILE1.py --AB 1.0 --tol 0.02 --profile 1 \
#  --maxiter 5 --maxfeval 10000 --mpi ${tmp} --outroot "monday" --minmethod 1