import sys, platform, os
sys.path.insert(0, os.environ['ROOTDIR'] + 
                   '/external_modules/code/CAMB/build/lib.linux-x86_64-'
                   +os.environ['PYTHON_VERSION'])
import functools
import numpy as np
import ipyparallel
import sys, platform, os
import math
import euclidemu2
import scipy
from getdist import IniFile
import itertools
import iminuit
import camb
import cosmolike_lsst_y1_interface as ci
import copy
import argparse
import random
import emcee
import itertools

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

CLprobe="xi"
path= os.environ['ROOTDIR'] + "/external_modules/data/lsst_y1"
data_file="lsst_y1_M1_GGL0.05.dataset"

IA_model = 0
IA_redshift_evolution = 3

# Init Cosmolike & Read LSST-Y1 data file
ini = IniFile(os.path.normpath(os.path.join(path, data_file)))
data_vector_file = ini.relativeFileName('data_file')
cov_file = ini.relativeFileName('cov_file')
mask_file = ini.relativeFileName('mask_file')
ntheta = ini.int("n_theta")
theta_min_arcmin = ini.float("theta_min_arcmin")
theta_max_arcmin = ini.float("theta_max_arcmin")

lens_file = ini.relativeFileName('nz_lens_file')

source_file = ini.relativeFileName('nz_source_file')

lens_ntomo = ini.int("lens_ntomo")

source_ntomo = ini.int("source_ntomo")

ci.initial_setup()

ci.init_cosmo_runmode(is_linear=False)

ci.init_source_sample(filename=source_file, ntomo_bins=int(source_ntomo))

ci.init_lens_sample(filename=lens_file, ntomo_bins=int(lens_ntomo))

ci.init_IA(ia_model=int(IA_model), ia_redshift_evolution=int(IA_redshift_evolution))

# Init Cosmolike
ci.init_probes(possible_probes=CLprobe)
ci.init_binning(int(ntheta), theta_min_arcmin, theta_max_arcmin)
ci.init_data_real(cov_file, mask_file, data_vector_file)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

def get_camb_cosmology(omegam, omegab, H0, ns, As_1e9 , w, w0pwa, mnu,
                       AccuracyBoost=1.0, kmax=10, k_per_logint=20, 
                       CAMBAccuracyBoost=1.1, non_linear_emul=2):
    As = lambda As_1e9: 1e-9 * As_1e9
    wa = lambda w0pwa, w: w0pwa - w
    omegabh2 = lambda omegab, H0: omegab*(H0/100)**2
    omegach2 = lambda omegam, omegab, mnu, H0: (omegam-omegab)*(H0/100)**2-(mnu*(3.046/3)**0.75)/94.0708
    omegamh2 = lambda omegam, H0: omegam*(H0/100)**2

    CAMBAccuracyBoost = CAMBAccuracyBoost*AccuracyBoost
    
    kmax = max(kmax/2.0, kmax*(1.0 + 3*(AccuracyBoost-1)))
    
    k_per_logint = max(k_per_logint/2.0, int(k_per_logint) + int(3*(AccuracyBoost-1)))
    
    extrap_kmax = max(max(2.5e2, 3*kmax), max(2.5e2, 3*kmax) * AccuracyBoost)

    z_interp_1D = np.concatenate( (np.concatenate((np.linspace(0,2.0,1000),
                                                   np.linspace(2.0,10.1,200)),
                                                   axis=0),
                                   np.linspace(1080,2000,20)),
                                   axis=0)
    
    z_interp_2D = np.concatenate(( np.linspace(0, 2.0, 95), 
                                   np.linspace(2.25, 10, 5)),  
                                 axis=0)

    log10k_interp_2D = np.linspace(-4.2, 2.0, 1200)

    pars = camb.set_params(H0=H0, 
                           ombh2=omegabh2(omegab, H0), 
                           omch2=omegach2(omegam, omegab, mnu, H0), 
                           mnu=mnu, 
                           omk=0, 
                           tau=0.06,  
                           As=As(As_1e9), 
                           ns=ns, 
                           halofit_version='takahashi', 
                           lmax=10,
                           AccuracyBoost=CAMBAccuracyBoost,
                           lens_potential_accuracy=1.0,
                           num_massive_neutrinos=1,
                           nnu=3.046,
                           accurate_massive_neutrino_transfers=False,
                           k_per_logint=k_per_logint,
                           kmax=kmax);
    
    pars.set_dark_energy(w=w, wa=wa(w0pwa, w), dark_energy_model='ppf');    
    
    pars.NonLinear = camb.model.NonLinear_both
    
    pars.set_matter_power(redshifts = z_interp_2D, kmax = kmax, silent = True);
    
    results = camb.get_results(pars)
    
    PKL  = results.get_matter_power_interpolator(var1="delta_tot", 
                                                 var2="delta_tot", 
                                                 nonlinear=False, 
                                                 extrap_kmax=extrap_kmax, 
                                                 hubble_units=False, 
                                                 k_hunit=False);
    
    PKNL = results.get_matter_power_interpolator(var1="delta_tot", 
                                                 var2="delta_tot",  
                                                 nonlinear=True, 
                                                 extrap_kmax=extrap_kmax, 
                                                 hubble_units=False, 
                                                 k_hunit=False);
    
    lnPL = np.empty(len(log10k_interp_2D)*len(z_interp_2D))
    for i in range(len(z_interp_2D)):
        lnPL[i::len(z_interp_2D)] = np.log(PKL.P(z_interp_2D[i], 
                                                 np.power(10.0,log10k_interp_2D)))
    lnPL  += np.log(((H0/100.)**3)) 
    
    lnPNL  = np.empty(len(log10k_interp_2D)*len(z_interp_2D))
    if non_linear_emul == 1:
        params = { 'Omm'  : omegam, 
                   'As'   : As(As_1e9), 
                   'Omb'  : omegab,
                   'ns'   : ns, 
                   'h'    : H0/100., 
                   'mnu'  : mnu,  
                   'w'    : w, 
                   'wa'   : wa(w0pwa, w)
                 }
        kbt, bt = euclidemu2.get_boost(params, 
                                       z_interp_2D, 
                                       np.power(10.0, 
                                                np.linspace(-2.0589, 
                                                            0.973, 
                                                            len(log10k_interp_2D)))
                                      )
        log10k_interp_2D = log10k_interp_2D - np.log10(H0/100.)
        
        for i in range(len(z_interp_2D)):    
            lnbt = scipy.interpolate.interp1d(np.log10(kbt), 
                                              np.log(bt[i]), 
                                              kind = 'linear', 
                                              fill_value = 'extrapolate', 
                                              assume_sorted = True)(log10k_interp_2D)
            lnbt[np.power(10,log10k_interp_2D) < 8.73e-3] = 0.0
            lnPNL[i::len(z_interp_2D)] = lnPL[i::len(z_interp_2D)] + lnbt
    elif non_linear_emul == 2:
        for i in range(len(z_interp_2D)):
            lnPNL[i::len(z_interp_2D)] = np.log(PKNL.P(z_interp_2D[i], 
                                                       np.power(10.0,log10k_interp_2D)))            
        log10k_interp_2D = log10k_interp_2D - np.log10(H0/100.)
        lnPNL += np.log(((H0/100.)**3))

    G_growth = np.sqrt(PKL.P(z_interp_2D,0.0005)/PKL.P(0,0.0005))
    G_growth = G_growth*(1 + z_interp_2D)/G_growth[len(G_growth)-1]

    chi = results.comoving_radial_distance(z_interp_1D, tol=1e-4) * (H0/100.)

    return (log10k_interp_2D, z_interp_2D, lnPL, lnPNL, G_growth, z_interp_1D, chi)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

def chi2(params, AccuracyBoost=1.0, non_linear_emul=2):

    As_1e9 = params[0]
    ns     = params[1]
    H0     = params[2]
    omegab = params[3]
    omegam = params[4]
    
    DZ_S1  = params[5]
    DZ_S2  = params[6]
    DZ_S3  = params[7]
    DZ_S4  = params[8]
    DZ_S5  = params[9]    
    
    A1_1   = params[10] 
    A1_2   = params[11]
    
    M1     = params[12]
    M2     = params[13]
    M3     = params[14]
    M4     = params[15]
    M5     = params[16]

    w      = -0.9 
    w0pwa  = -0.9
    mnu    =  0.06

    (log10k_interp_2D, z_interp_2D, lnPL, lnPNL, 
        G_growth, z_interp_1D, chi) = get_camb_cosmology(omegam=omegam, 
                                                         omegab=omegab, 
                                                         H0=H0, 
                                                         ns=ns, 
                                                         As_1e9=As_1e9,
                                                         w=w, 
                                                         w0pwa=w0pwa, 
                                                         mnu=mnu,
                                                         AccuracyBoost=AccuracyBoost,
                                                         non_linear_emul=non_linear_emul)
    
    ci.init_accuracy_boost(AccuracyBoost, 
                           AccuracyBoost, 
                           int(1+5*(AccuracyBoost-1)))
    
    ci.set_cosmology(omegam=omegam,
                     H0=H0, 
                     log10k_2D=log10k_interp_2D, 
                     z_2D=z_interp_2D, 
                     lnP_linear=lnPL,
                     lnP_nonlinear=lnPNL,
                     G=G_growth,
                     z_1D=z_interp_1D,
                     chi=chi)
    
    ci.set_nuisance_shear_calib(M=[M1, M2, M3, M4, M5])
    
    ci.set_nuisance_shear_photoz(bias=[DZ_S1, DZ_S2, DZ_S3, DZ_S4, DZ_S5])
    
    ci.set_nuisance_ia(A1   = [A1_1, A1_2, 0, 0, 0], 
                       A2   = [0, 0, 0, 0, 0], 
                       B_TA = [0, 0, 0, 0, 0])

    datavector = np.array(ci.compute_data_vector_masked())
    return ci.compute_chi2(datavector)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

x = np.array([
                2.1,          # As
                0.96,         # ns
                67.0,         # H0
                0.04,         # omegab
                0.30,         # omegam
                0.04,         # S1
                0.0016,       # S2
                0.03,         # S3
                -0.08,        # S4
                -8.67127e-05, # S5
                0.7,          # A11
                -1.5,         # A12
                0.001,        # M1
                0.002,        # M2
                0.003,        # M3
                0.004,        # M4
                0.001         # M5
              ], dtype='float64')

bounds = np.array([
                    [1.95, 2.4],   # As
                    [0.89, 1.05],  # ns 
                    [55.0, 87.0],  # H0
                    [0.027, 0.065], # omegab
                    [0.22, 0.38],   # omegam
                    [-0.12, 0.12], # S1
                    [-0.12, 0.12], # S2
                    [-0.12, 0.12], # S3
                    [-0.12, 0.12], # S4
                    [-0.12, 0.12], # S5
                    [-4.5, 4.5],   # A11
                    [-4.5, 4.5],   # A12
                    [-0.12, 0.12], # M1
                    [-0.12, 0.12], # M2
                    [-0.12, 0.12], # M3
                    [-0.12, 0.12], # M4
                    [-0.12, 0.12]  # M5
                  ], dtype='float64')

start = np.array([ 
                    1.96,         # As
                    0.90,         # ns 
                    61.0,         # H0
                    0.028,        # omegab
                    0.25,         # omegam
                    -0.20,        # S1
                    -0.20,        # S2
                    -0.20,        # S3
                    -0.20,        # S4
                    -0.20,        # S5                   
                    -3.00,        # A11
                    -3.00,        # A12
                    -0.20,        # M1
                    -0.20,        # M2
                    -0.20,        # M3
                    -0.20,        # M4
                    -0.20         # M5
                ], dtype='float64')

stop  = np.array([ 
                    2.3,         # As
                    1.04,         # ns
                    78.0,         # H0
                    0.05,         # omegab
                    0.35,         # omegam
                    0.20,         # S1
                    0.20,         # S2
                    0.20,         # S3
                    0.20,         # S4                    
                    0.20,         # S5                    
                    3.00,         # A11
                    3.00,         # A12
                    0.20,         # M1
                    0.20,         # M2
                    0.20,         # M3
                    0.20,         # M4
                    0.20          # M5
                 ], dtype='float64')

name  = [ 
            "As",       # As
            "ns",       # ns
            "H0",       # H0
            "omegab",   # omegab
            "omm",      # omegam
            "S1",       # S1
            "S2",       # S2
            "S3",       # S3
            "S4",       # S4
            "S5",       # S5          
            "A11",      # A11
            "A12",      # A12
            "M1",       # M1
            "M2",       # M2
            "M3",       # M3
            "M4",       # M4
            "M5"        # M5
        ]

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

parser = argparse.ArgumentParser(prog='EXAMPLE_PROFILE1')

parser.add_argument("--AB",
                    dest="AccuracyBoost",
                    help="Accuracy Boost of CAMB/Cosmolike calculation",
                    type=float,
                    nargs='?',
                    const=1,
                    default=1.1)

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
                    default=50000)

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
                    default=1)

parser.add_argument("--outroot",
                    dest="outroot",
                    help="Name of the Output File",
                    nargs='?',
                    const=1)

parser.add_argument("--profile",
                    dest="profile",
                    help="Which Parameter to Profile",
                    type=int,
                    nargs='?',
                    const=1)

parser.add_argument("--mpi",
                    dest="mpi",
                    help="Number of MPI cores",
                    type=int,
                    nargs='?',
                    const=1)

# need to use parse_known_args because of mpifuture 
args, unknown = parser.parse_known_args()

non_linear_emul = 2
AccuracyBoost   = args.AccuracyBoost
tol             = args.tolerance
maxfeval        = args.maxfeval
maxiter         = args.maxiter
min_method      = args.minmethod
oroot           = "chains/" + args.outroot
index           = args.profile
nummpi          = args.mpi

cov_file = 'EXAMPLE_MCMC1.covmat'
cov      = np.loadtxt(cov_file)
# need to delete the w line
cov = np.delete(cov, (5), axis=0)
cov = np.delete(cov, (5), axis=1)
cov = np.array(cov, dtype='float64')
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

def min_chi2(x0, bounds, min_method, fixed=-1, AccuracyBoost=1.0, 
             tol=0.01, maxfeval=3000, non_linear_emul=2, maxiter=10, cov=cov):

    def mychi2(params, *args):
        z, fixed, T = args
        params = np.array(params)
        if fixed > -1:
            params = np.insert(params, fixed, z)
        res = chi2(params=params, 
                   AccuracyBoost=AccuracyBoost, 
                   non_linear_emul=non_linear_emul)/T
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
            x = [] # Initial point
            for j in range(nwalkers):
                x.append(GaussianStep(stepsize=stepsz[i])(x0)[0,:])
            x = np.array(x, dtype='float64')

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
                    AccuracyBoost=AccuracyBoost, 
                    tol=tol, 
                    maxfeval=maxfeval, 
                    non_linear_emul=non_linear_emul, 
                    maxiter=maxiter)
    return res


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

from mpi4py.futures import MPIPoolExecutor

if __name__ == '__main__':
    
    print(f"min_method={min_method}, maxiter={maxiter}, maxfeval={maxfeval}, tol={tol}, param={index}")

    executor = MPIPoolExecutor()

    param = np.linspace(start=start[index], stop=stop[index], num=nummpi)
    
    print(f"profile param values = {param}")

    x0 = []
    for i in range(nummpi):
        tmp = copy.deepcopy(x)
        tmp[index] = param[i]
        x0.append(tmp)
    x0 = np.array(x0, dtype=object)

    res = np.array(list(executor.map(functools.partial(prf,
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
    
    executor.shutdown()
