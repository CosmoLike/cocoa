import warnings
import os, psutil
from sklearn.exceptions import InconsistentVersionWarning
warnings.filterwarnings("ignore", category=InconsistentVersionWarning)
warnings.filterwarnings(
    "ignore",
    message=".*column is deprecated.*",
    module=r"sacc\.sacc"
)
warnings.filterwarnings(
    "ignore",
    category=RuntimeWarning,
    message=r".*invalid value encountered.*"
)
warnings.filterwarnings(
    "ignore",
    category=RuntimeWarning,
    message=r".*overflow encountered*"
)
warnings.filterwarnings(
    "ignore",
    category=UserWarning,
    message=r".*Hartlap correction*"
)
import functools, iminuit, copy, argparse, random, time 
import emcee, itertools
import numpy as np
from cobaya.yaml import yaml_load
from cobaya.model import get_model
from getdist import IniFile
from schwimmbad import MPIPool
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
_affinity_set = False

def enforce_affinity():
    # schwimmbad.MPIPool uses one Python process per MPI rank
    # So each rank can directly control its own CPU affinity via psutil
    # No reliance on mpirun or OpenMPI doing the right thing
    rank = int(os.environ.get("OMPI_COMM_WORLD_RANK", 0))
    omp_threads = int(os.environ.get("OMP_NUM_THREADS", 1))
    first_core = rank * omp_threads
    last_core  = first_core + omp_threads - 1
    try:
        psutil.Process().cpu_affinity(list(range(first_core, last_core + 1)))
    except Exception as e:
        print(f"[Rank {rank}] Failed to set affinity: {e}")
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
parser = argparse.ArgumentParser(prog='EXAMPLE_MINIMIZE_LCDM_CAMB1')
parser.add_argument("--nstw",
                    dest="nstw",
                    help="Number of likelihood evaluations (steps) per temperature per walker",
                    type=int,
                    nargs='?',
                    const=1,
                    default=200)
parser.add_argument("--root",
                    dest="root",
                    help="Name of the Output File",
                    nargs='?',
                    const=1,
                    default="./projects/lsst_y1/")
parser.add_argument("--outroot",
                    dest="outroot",
                    help="Name of the Output File",
                    nargs='?',
                    const=1,
                    default="example_min1")
# need to use parse_known_args because of mpifuture 
args, unknown = parser.parse_known_args() 
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
yaml_string = r"""
likelihood:
  planck_2018_highl_plik.TTTEEE_lite: 
    path: ./external_modules/
    clik_file: plc_3.0/hi_l/plik_lite/plik_lite_v22_TTTEEE.clik
  planck_2018_lowl.TT: 
    path: ./external_modules
  planck_2018_lowl.EE:
    path: ./external_modules
  sn.desy5: 
    path: ./external_modules/data/sn_data
  bao.desi_dr2.desi_bao_all:
    path: ./external_modules/data/ 
  act_dr6_lenslike.ACTDR6LensLike:
    lens_only: True
    variant: actplanck_baseline
params:
  logA:
    prior:
      min: 1.61
      max: 3.91
    ref:
      dist: norm
      loc: 3.0448
      scale: 0.05
    proposal: 0.05
    latex: '\log(10^{10} A_\mathrm{s})'
  ns:
    prior:
      min: 0.92
      max: 1.05
    ref:
      dist: norm
      loc: 0.96605
      scale: 0.005
    proposal: 0.005
    latex: 'n_\mathrm{s}'
  theta_MC_100:
    prior:
      min: 1
      max: 1.2
    ref:
      dist: norm
      loc: 1.04109
      scale: 0.0004
    proposal: 0.0002
    latex: 100\theta_\mathrm{MC}
    drop: true
    renames: theta
  cosmomc_theta:
    value: 'lambda theta_MC_100: 1.e-2*theta_MC_100'
    derived: false
  omegabh2:
    prior:
      min: 0.01
      max: 0.04
    ref:
      dist: norm
      loc: 0.022383
      scale: 0.005
    proposal: 0.005
    latex: '\Omega_\mathrm{b} h^2'
  omegach2:
    prior:
      min: 0.06
      max: 0.2
    ref:
      dist: norm
      loc: 0.12011
      scale: 0.03
    proposal: 0.03
    latex: '\Omega_\mathrm{c} h^2'
  tau:
    prior:
      min: 0.04
      max: 0.09
    ref:
      dist: norm
      loc: 0.055
      scale: 0.01
    proposal: 0.005
    latex: \tau_\mathrm{reio}
  mnu:
    value: 0.06
  As:
    value: 'lambda logA: 1e-10*np.exp(logA)'
    latex: A_\mathrm{s}
  omegab:
    derived: 'lambda omegabh2, H0: omegabh2/((H0/100)**2)'
    latex: \Omega_\mathrm{b}
  omegac:
    derived: 'lambda omegach2, H0: omegach2/((H0/100)**2)'
    latex: \Omega_\mathrm{c}
  H0:
    derived: true
    latex: H_0
  omegam:
    derived: true
    latex: \Omega_\mathrm{m}
  rdrag:
    derived: true
    latex: r_\mathrm{drag}
  omegamh2:
    derived: 'lambda omegach2, omegabh2, mnu: omegach2+omegabh2+(mnu*(3.046/3)**0.75)/94.0708'
    latex: \Omega_\mathrm{m} h^2
  thetastar:
   derived: true
   latex: \Theta_\star
theory:
  camb:
    path: ./external_modules/code/CAMB
    extra_args:
      AccuracyBoost: 1.10
      lens_potential_accuracy: 2
      num_massive_neutrinos: 1
      nnu: 3.046
      nonlinear: true
      halofit_version: takahashi
      dark_energy_model: ppf
      kmax: 5
"""
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
model = get_model(yaml_load(yaml_string))
def chi2(p):
    p = [float(v) for v in p.values()] if isinstance(p, dict) else p
    if np.any(np.isinf(p)) or  np.any(np.isnan(p)):
      raise ValueError(f"At least one parameter value was infinite (CoCoa) param = {p}")
    point = dict(zip(model.parameterization.sampled_params(), p))
    res1 = model.logprior(point,
                          make_finite=False)
    if np.isinf(res1) or np.any(np.isnan(res1)):
      return 1.e20
    res2 = model.loglike(point,
                         make_finite=False,
                         cached=False,
                         return_derived=False)
    if np.isinf(res2) or np.isnan(res2):
      return 1.e20
    return -2.0*(res1+res2)
def chi2v2(p):
    p = [float(v) for v in p.values()] if isinstance(p, dict) else p
    point = dict(zip(model.parameterization.sampled_params(), p))
    logposterior = model.logposterior(point, as_dict=True)
    chi2likes=-2*np.array(list(logposterior["loglikes"].values()))
    chi2prior=-2*np.atleast_1d(model.logprior(point,make_finite=False))
    return np.concatenate((chi2likes, chi2prior))
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
def min_chi2(x0, 
             cov, 
             fixed=-1, 
             nstw=200,
             nwalkers=5,
             pool=None):
    def mychi2(params, *args):
        z, fixed, T = args
        params = np.array(params, dtype='float64')
        if fixed > -1:
            params = np.insert(params, fixed, z)
        return chi2(p=params)/T

    if fixed > -1:
        z      = x0[fixed]
        x0     = np.delete(x0, (fixed))
        args = (z, fixed, 1.0)
        
        cov = np.delete(cov, (fixed), axis=0)
        cov = np.delete(cov, (fixed), axis=1)
    else:
        args = (0.0, -2.0, 1.0)

    def logprob(params, *args):
        global _affinity_set
        if not _affinity_set:
          enforce_affinity()  # enforce per-rank affinity on pool workers!
          _affinity_set = True
          start_time = time.time()
          res = mychi2(params, *args)
          etime = time.time() - start_time
          rank = int(os.environ.get("OMPI_COMM_WORLD_RANK",0))
          print(f"Emcee: Like Eval Time: {etime:.4f} secs and MPI Rank: {rank}")
        else:
          res = mychi2(params, *args)
        if (res > 1.e19 or np.isinf(res) or  np.isnan(res)):
          return -np.inf
        else:
          return -0.5*res
    
    class GaussianStep:
       def __init__(self, stepsize=0.2):
           self.cov = stepsize*cov
       def __call__(self, x):
           return np.random.multivariate_normal(x, self.cov, size=1)   
    
    ndim        = int(x0.shape[0])
    nwalkers    = int(nwalkers)
    nstw        = int(nstw)
    #temperature = np.array([1.0, 0.25, 0.1, 0.005, 0.001], dtype='float64')
    temperature = np.array([1.0], dtype='float64')
    ntemp       = len(temperature)
    stepsz      = temperature/4.0

    start_time = time.time()
    mychi2(GaussianStep(stepsize=0.001)(x0)[0,:], *args)
    elapsed_time = time.time() - start_time
    print(f"nTemp = {len(temperature)}, "
          f"feval (per Temp) = {nstw}, "
          f"feval = {nstw*len(temperature)}")
    print(f"Emcee: Like Eval Time: {elapsed_time:.4f} secs, "
          f"Eval Time = {elapsed_time*nstw*ntemp/3600.:.4f} hours.")

    partial_samples = []
    partial = []
    for i in range(len(temperature)):
        x = [] # Initial point
        for j in range(nwalkers):
            x.append(GaussianStep(stepsize=stepsz[i])(x0)[0,:])  
        sampler = emcee.EnsembleSampler(nwalkers, 
                                        ndim, 
                                        logprob, 
                                        args=(args[0], args[1], temperature[i]),
                                        moves=[(emcee.moves.DEMove(), 0.8),
                                               (emcee.moves.DESnookerMove(), 0.2)],
                                        pool=pool)    
        sampler.run_mcmc(np.array(x, dtype='float64'), 
                         nstw, 
                         skip_initial_state_check=True)
        samples = sampler.get_chain(flat=True, discard=0)
        j = np.argmin(-1.0*np.array(sampler.get_log_prob(flat=True)))
        partial_samples.append(samples[j])
        partial.append(mychi2(samples[j], *args))
        x0 = copy.deepcopy(samples[j])
        sampler.reset()    
        j = np.argmin(np.array(partial))
        print(f"Partial ({i+1}/{len(temperature)}): "
              f"params = {partial_samples[j]}, and chi2 = {partial[j]}")
    # min chi2 from the entire emcee runs
    j = np.argmin(np.array(partial))
    result = [partial_samples[j], partial[j]]
    return result
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
def prf(x0, nstw, cov, fixed=-1, nwalkers=5, pool=None):
    res =  min_chi2(x0=np.array(x0, dtype='float64'), 
                    fixed=fixed,
                    cov=cov, 
                    nstw=nstw, 
                    nwalkers=nwalkers,
                    pool=pool)
    return res
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
from cobaya.theories.emultheta.emultheta2 import emultheta
etheta = emultheta(extra_args={ 
    'file': ['external_modules/data/emultrf/CMB_TRF/emul_lcdm_thetaH0_GP.joblib'],
    'extra': ['external_modules/data/emultrf/CMB_TRF/extra_lcdm_thetaH0.npy'],
    'ord': [['omegabh2','omegach2','thetastar']],
    'extrapar': [{'MLA' : "GP"}]})
from cobaya.theories.emulrdrag.emulrdrag2 import emulrdrag
erd = emulrdrag(extra_args={ 
    'file': ['external_modules/data/emultrf/BAO_SN_RES/emul_lcdm_rdrag_GP.joblib'],
    'extra':['external_modules/data/emultrf/BAO_SN_RES/extra_lcdm_rdrag.npy'],
    'ord':  [['omegabh2','omegach2']]})
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if __name__ == '__main__':
    with MPIPool() as pool:
        if not pool.is_master():
            pool.wait()
            sys.exit(0)
        enforce_affinity() # enforce affinity (so Hybrid MPI-OpenMP works)!
        
        dim = model.prior.d()  
        nwalkers = max(3*dim, pool.comm.Get_size())
        nstw = args.nstw
        (x0, results) = model.get_valid_point(max_tries=50, 
                                              ignore_fixed_ref=False,
                                              logposterior_as_dict=True)
        # 1st: Get covariance --------------------------------------------------
        cov = model.prior.covmat(ignore_external=False) # cov from prior
        
        # 2nd: Run Procoli -----------------------------------------------------
        res = np.array(list(prf(np.array(x0, dtype='float64'), 
                               fixed=-1, 
                               nstw=nstw,
                               nwalkers=nwalkers,
                               pool=pool,
                               cov=cov)), dtype="object")
        xf = np.array([res[0]],dtype='float64')
        
        # 3rd Append derived parameters ----------------------------------------
        H0 = []
        omm = []
        rdrag = []
        for d in xf:
            tmp = model.logposterior(d, 
                                     as_dict=True, 
                                     make_finite=True, 
                                     return_derived=True)
            H0.append(tmp['derived']['H0'])
            omm.append(tmp['derived']['omegam'])
            rdrag.append(tmp['derived']['rdrag'])
        xf = np.column_stack((xf, 
                              np.array(H0,dtype='float64'), 
                              np.array(omm,dtype='float64'),
                              np.array(rdrag,dtype='float64'), 
                              np.array([chi2v2(d) for d in xf], dtype='float64'),
                              res[1]))
  
        # 4th Save output file -------------------------------------------------
        os.makedirs(os.path.dirname(f"{args.root}chains/"), exist_ok=True)
        names = list(model.parameterization.sampled_params().keys()) # Cobaya Call
        hd = names + ['H0', 'omegam', 'rdrag']
        hd = hd + list(model.info()['likelihood'].keys()) + ["prior"] + ["chi2"]
        np.savetxt(f"{args.root}chains/{args.outroot}.txt", 
                   xf,
                   fmt="%.12e",
                   header=f"nswt (evals/Temp/walker)={nstw}\n"+' '.join(hd),
                   comments="# ")
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------