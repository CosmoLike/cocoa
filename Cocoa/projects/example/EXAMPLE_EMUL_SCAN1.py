r"""Scan chi2 along one parameter over its whole prior, one minimization per process.

A scan is a profile (EXAMPLE_EMUL_PROFILE1.py defines it): one parameter is fixed
at each value of a grid, and chi2 = -2 log posterior = -2 (log prior + log
likelihood) is minimized over all the other parameters. Here the grid spans the
whole prior range of the profiled parameter, and the MPI work is split
differently: each worker process runs a complete minimization for one grid point,
evaluating its walkers itself, instead of all workers sharing the walkers of one
minimization. Every grid point starts from its own random point, so the scan
suits parameters whose chi2 can have several separated minima, such as those of
oscillatory (monodromic) dark-energy models.

As written, the first minimization stops with a TypeError (min_chi2 explains
why), so the script produces no output; the text below describes the code
around that call.

The model is the YAML in yaml_string, the same as in EXAMPLE_EMUL_MINIMIZE1.py:
Planck 2018 high-l TT, TE, EE (plik-lite), low-l TT and low-l EE (SRoll2), DES-Y5
supernovae, DESI DR2 BAO and ACT DR6 lensing, in LCDM with the neutrino mass
fixed at 0.06 eV. Emulators (neural networks and Gaussian processes trained on
Boltzmann-code outputs) replace CAMB. Seven sampled parameters, in sampled order:
logA, ns, thetastar = 100 theta_*, omegabh2, omegach2, tau, and the Planck
calibration A_planck (Gaussian prior, mean 1, standard deviation 0.0025).

Each minimization is annealed emcee sampling: a set of points called walkers
samples the posterior raised to the power 1/T while the temperature T falls
along the ladder 1.0, 0.25, 0.1, 0.005, 0.001, so the walkers crowd around the
minimum (min_chi2 gives the details; the module docstring of
external_modules/code/cosmolike_core/cocoa_hybrid_sampling.py explains the same
annealing scheme). Settings: 3 * n_sampled walkers, the prior covariance, and
starting clouds of covariance cov*T/4.

The grid has one point per MPI process, comm.Get_size(), evenly spaced over the
central 0.999999 interval of the profiled parameter's prior (its full range for a
flat prior). mpi4py.futures runs the main block on rank 0 and the minimizations
on the other ranks, so N processes give N grid points for N - 1 workers.

Inputs (command line; paths are relative to the Cocoa folder):
  --root     folder, ending in "/", that receives chains/ (default
             ./projects/example/).
  --outroot  basename of the output file (default test.dat).
  --profile  zero-based index, in sampled order, of the scanned parameter
             (default 1, ns).
  --nstw     emcee steps per walker at each temperature (default 200).

Output: <root>chains/<outroot>.<name>.txt, with <name> the scanned parameter,
one row per grid point: the grid value, the minimized chi2, the sampled
parameters of the minimization, the chi2 of each likelihood, and -2 log prior.
The header records maxfeval = nstw * 5 * nwalkers, the chi2 evaluations of one
grid point.

Run from the Cocoa folder after "source start_cocoa.sh"; python -m
mpi4py.futures starts the worker processes:

    mpirun -n 12 --oversubscribe python -m mpi4py.futures \
        ./projects/example/EXAMPLE_EMUL_SCAN1.py --root ./projects/example/ \
        --outroot EXAMPLE_EMUL_SCAN1 --nstw 200 --profile 1
"""

import warnings
import os
from sklearn.exceptions import InconsistentVersionWarning
# The filters below hide messages; they change no computed value. They hide the
# scikit-learn notice raised when a Gaussian-process emulator file (joblib) was
# saved by another scikit-learn version, a deprecation notice of the sacc data
# library, numpy's "overflow" and "invalid value" warnings at extreme parameter
# values (chi2 turns those points into rejections), a "Function not smooth or
# differentiable" notice, and the ACT DR6 lensing message that prints its Hartlap
# factor (the correction for an inverse covariance estimated from a finite
# number of simulations).
warnings.filterwarnings("ignore", category=InconsistentVersionWarning)
warnings.filterwarnings(
    "ignore",
    message=".*column is deprecated.*",
    module=r"sacc\.sacc"
)
warnings.filterwarnings(
    "ignore",
    category=RuntimeWarning,
    message=r".*invalid value encountered*"
)
warnings.filterwarnings(
    "ignore",
    category=RuntimeWarning,
    message=r".*overflow encountered*"
)
warnings.filterwarnings(
    "ignore",
    category=UserWarning,
    message=r".*Function not smooth or differentiabl*"
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
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Command-line options; the module docstring gives their meaning. Every option
# uses nargs='?' with const=1: written without a value, an option takes the value
# 1 (for --root that is the integer 1, not a path), so always give a value.
parser = argparse.ArgumentParser(prog='EXAMPLE_EMUL_SCAN1')
parser.add_argument("--nstw",
                    dest="nstw",
                    help="Number of likelihood evaluations per temperature per walker",
                    type=int,
                    nargs='?',
                    const=1,
                    default=200)
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
# parse_known_args, unlike parse_args, does not stop at an argument it does not
# know: it returns such arguments in `unknown`, which nothing reads. The examples
# keep it for runs launched through mpi4py.futures. A misspelled option is
# therefore ignored without a message, and its default applies.
args, unknown = parser.parse_known_args()
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# The model, in Cobaya's YAML format: the likelihoods; the parameters (a sampled
# parameter has a prior; "ref" is the distribution that random starting points
# are drawn from; "proposal" is read by Cobaya's own MCMC sampler, not by this
# script); and the theory block, where the emulators replace CAMB. Paths such as
# ./external_modules are relative to the Cocoa folder, the folder to run from.
yaml_string=r"""
likelihood:
  planck_2018_highl_plik.TTTEEE_lite: 
    path: ./external_modules/
    clik_file: plc_3.0/hi_l/plik_lite/plik_lite_v22_TTTEEE.clik
  planck_2018_lowl.TT: 
    path: ./external_modules
  # choose only one low ell EE likelihood
  #planck_2018_lowl.EE:
  #  path: ./external_modules
  planck_2018_lowl.EE_sroll2: null 
  #planck_2020_lollipop.lowlE:
  #  data_folder: planck/lollipop
  sn.desy5: 
    path: ./external_modules/data/sn_data
  bao.desi_dr2.desi_bao_all:
    path: ./external_modules/data/ 
  act_dr6_lenslike.ACTDR6LensLike:
    lens_only: False
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
    latex: '\log(10^{10} A_\mathrm{s}'
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
  thetastar:
    prior:
      min: 1
      max: 1.2
    ref:
      dist: norm
      loc: 1.04109
      scale: 0.0004
    proposal: 0.0002
    latex: '100\theta_\mathrm{*}'
    renames: theta
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
    proposal: 0.01
    latex: \tau_\mathrm{reio}
  As:
    derived: 'lambda logA: 1e-10*np.exp(logA)'
    latex: 'A_\mathrm{s}'
  A:
    derived: 'lambda As: 1e9*As'
    latex: '10^9 A_\mathrm{s}'
  mnu:
    value: 0.06
  w0pwa:
    value: -1.0
    latex: 'w_{0,\mathrm{DE}}+w_{a,\mathrm{DE}}'
    drop: true
  w:
    value: -1.0
    latex: 'w_{0,\mathrm{DE}}'
  wa:
    value: 'lambda w0pwa, w: w0pwa - w'
    derived: false
    latex: 'w_{a,\mathrm{DE}}'
  H0:
    latex: H_0
  omegamh2:
    derived: true
    value: 'lambda omegach2, omegabh2, mnu: omegach2+omegabh2+(mnu*(3.046/3)**0.75)/94.0708'
    latex: '\Omega_\mathrm{m} h^2'
  omegam:
    latex: '\Omega_\mathrm{m}'
  rdrag:
    latex: 'r_\mathrm{drag}'
theory:
  emultheta:
    path: ./cobaya/cobaya/theories/
    provides: ['H0', 'omegam']
    extra_args:
      file: ['external_modules/data/emultrf/CMB_TRF/emul_lcdm_thetaH0_GP.joblib']
      extra: ['external_modules/data/emultrf/CMB_TRF/extra_lcdm_thetaH0.npy']
      ord: [['omegabh2','omegach2','thetastar']]
      extrapar: [{'MLA' : "GP"}]
  emulrdrag:
    path: ./cobaya/cobaya/theories/
    provides: ['rdrag']
    extra_args:
      file: ['external_modules/data/emultrf/BAO_SN_RES/emul_lcdm_rdrag_GP.joblib'] 
      extra: ['external_modules/data/emultrf/BAO_SN_RES/extra_lcdm_rdrag.npy'] 
      ord: [['omegabh2','omegach2']]
  emulcmb:
    path: ./cobaya/cobaya/theories/
    extra_args:
      # This version of the emul was not trained with CosmoRec
      eval: [True, True, True, True] #TT,TE,EE,PHIPHI
      device: "cuda"
      ord: [['omegabh2','omegach2','H0','tau','ns','logA','mnu','w','wa'],
            ['omegabh2','omegach2','H0','tau','ns','logA','mnu','w','wa'],
            ['omegabh2','omegach2','H0','tau','ns','logA','mnu','w','wa'],
            ['omegabh2','omegach2','H0','tau','ns','logA','mnu','w','wa']]
      file: ['external_modules/data/emultrf/CMB_TRF/emul_lcdm_CMBTT_CNN.pt',
             'external_modules/data/emultrf/CMB_TRF/emul_lcdm_CMBTE_CNN.pt',
             'external_modules/data/emultrf/CMB_TRF/emul_lcdm_CMBEE_CNN.pt', 
             'external_modules/data/emultrf/CMB_TRF/emul_lcdm_phi_ResMLP.pt']
      extra: ['external_modules/data/emultrf/CMB_TRF/extra_lcdm_CMBTT_CNN.npy',
              'external_modules/data/emultrf/CMB_TRF/extra_lcdm_CMBTE_CNN.npy',
              'external_modules/data/emultrf/CMB_TRF/extra_lcdm_CMBEE_CNN.npy', 
              'external_modules/data/emultrf/CMB_TRF/extra_lcdm_phi_ResMLP.npy']
      extrapar: [{'ellmax' : 5000, 'MLA': 'CNN', 'INTDIM': 4, 'INTCNN': 5120},
                 {'ellmax' : 5000, 'MLA': 'CNN', 'INTDIM': 4, 'INTCNN': 5120},
                 {'ellmax' : 5000, 'MLA': 'CNN', 'INTDIM': 4, 'INTCNN': 5120}, 
                 {'MLA': 'ResMLP', 'INTDIM': 4, 'NLAYER': 4, 
                  'TMAT': 'external_modules/data/emultrf/CMB_TRF/PCA_lcdm_phi.npy'}]
  emulbaosn:
    path: ./cobaya/cobaya/theories/
    stop_at_error: True
    extra_args:
      device: "cuda"
      file:  [None, 'external_modules/data/emultrf/BAO_SN_RES/emul_lcdm_H.pt']
      extra: [None, 'external_modules/data/emultrf/BAO_SN_RES/extra_lcdm_H.npy']    
      ord: [None, ['omegam','H0']]
      extrapar: [{'MLA': 'INT', 'ZMIN' : 0.0001, 'ZMAX' : 3, 'NZ' : 600},
                 {'MLA': 'ResMLP', 'offset' : 0.0, 'INTDIM' : 1, 'NLAYER' : 1,
                  'TMAT': 'external_modules/data/emultrf/BAO_SN_RES/PCA_lcdm_H.npy',
                  'ZLIN': 'external_modules/data/emultrf/BAO_SN_RES/z_lin_lcdm.npy'}]
"""
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Module-level code runs on every MPI process (mpi4py.futures workers load this
# file too), so each process builds its own Cobaya model and runs its
# minimizations with it.
model = get_model(yaml_load(yaml_string))
def chi2(p):
    """Return chi2 = -2 log posterior at one point, or 1e20 for a rejected point.

    Every script of this project calls -2 (log prior + log likelihood) "chi2".
    Cobaya's log prior keeps the normalization of each prior: -ln(max - min) for a
    flat prior, and the Gaussian density for A_planck. chi2 therefore differs from
    the likelihood chi2 = -2 log likelihood by a constant plus the A_planck prior
    penalty; the constant cancels in every chi2 difference.

    1e20 is the rejection value: a point outside the prior, or with a non-finite
    likelihood, gets it, and logprob in min_chi2 maps any chi2 above 1e19 to -inf.

    Arguments:
        p = sampled parameter values in the model's sampled order: a list or an
            array [n_sampled], or a dict whose values are in that order.
    Returns:
        chi2 as a float, or 1e20.
    Raises:
        ValueError when a parameter value is infinite or NaN.
    """
    # A dict is replaced by the list of its values (which must be in sampled
    # order); a list or an array is used as given.
    p = [float(v) for v in p.values()] if isinstance(p, dict) else p
    if np.any(np.isinf(p)) or  np.any(np.isnan(p)):
      raise ValueError(f"At least one parameter value was infinite (CoCoa) param = {p}")
    # zip pairs each sampled-parameter name with the value in the same position
    # and stops at the shorter input, so p must hold one value per parameter.
    point = dict(zip(model.parameterization.sampled_params(), p))
    # make_finite=False keeps -inf (outside the prior) instead of replacing it by
    # the most negative float.
    res1 = model.logprior(point,make_finite=False)
    if np.isinf(res1) or  np.any(np.isnan(res1)):
      return 1e20
    # cached=False recomputes the theory even if Cobaya holds a result for these
    # parameters; return_derived=False skips the derived parameters.
    res2 = model.loglike(point,
                         make_finite=False,
                         cached=False,
                         return_derived=False)
    if np.isinf(res2) or  np.any(np.isnan(res2)):
      return 1e20
    return -2.0*(res1+res2)
def chi2v2(p):
    """Return the chi2 of each likelihood followed by the prior term, at one point.

    Arguments:
        p = sampled parameter values, as for chi2.
    Returns:
        array [n_likelihoods + 1]: -2 log likelihood of each likelihood, in the
        order of model.info()["likelihood"] (the order the output header uses),
        then -2 log prior. For a point inside the prior the entries sum to chi2.
    """
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
             maxfeval=3000, 
             nwalkers=5):
    """Minimize chi2 by annealed emcee sampling in this process; return the best.

    At temperature T the walkers sample the posterior raised to the power 1/T,
    whose log is -chi2/(2T). For a nearly Gaussian posterior they sit about
    n_dim*T above the minimum in chi2. Each stage of the ladder 1.0, 0.25, 0.1,
    0.005, 0.001 draws nwalkers starting points from a Gaussian centered on x0
    with covariance cov*T/4 (one quarter of the variance of the tempered
    posterior when cov is close to the posterior covariance), runs nstw emcee
    steps, and passes its best sample to the next stage. The function returns
    the best of the starting point and the stage optima. No pool is given to
    emcee: the walkers are evaluated one after the other in this process.

    The body reads nstw (emcee steps per walker at each temperature), which is
    not among the parameters, and the signature names maxfeval, which the body
    never reads. prf passes nstw=, so as written the call raises TypeError
    (unexpected keyword argument 'nstw') before any sampling.

    Arguments:
        x0 = starting point [n_sampled], in the model's sampled order.
        cov = covariance [n_sampled, n_sampled] in the same order; the row and the
            column of the fixed parameter are removed before walkers are drawn.
        fixed = index of the parameter held at x0[fixed], or -1 for none.
        maxfeval = not read by the body.
        nwalkers = number of walkers.
    Returns:
        [best point [n_dim], its chi2], where n_dim = n_sampled, or n_sampled - 1
        when a parameter is fixed (the fixed value is not included).
    """
    def mychi2(params, *args):
        """Return chi2/T at one walker position, with the fixed coordinate put back.

        Arguments:
            params = walker position [n_dim], without the fixed coordinate.
            args = (value of the fixed parameter, its index or a negative number
                for none, temperature T).
        Returns:
            chi2/T as a float (1e20/T for a rejected point).
        """
        z, fixed, T = args
        params = np.array(params, dtype='float64')
        if fixed > -1:
            params = np.insert(params, fixed, z)
        return chi2(p=params)/T
    # A fixed parameter leaves x0 and cov; mychi2 puts its value back. Without
    # one, args carries -2 as the index, which mychi2 reads as "none".
    if fixed > -1:
        z      = x0[fixed]
        x0     = np.delete(x0, (fixed))
        args = (z, fixed, 1.0)
        cov = np.delete(cov, (fixed), axis=0)
        cov = np.delete(cov, (fixed), axis=1)
    else:
        args = (0.0, -2.0, 1.0)

    def logprob(params, *args):
        """Return -chi2/(2T), the log of the tempered posterior emcee samples.

        A chi2 above 1e19 (the rejection value of chi2) becomes -inf, and emcee
        never accepts a move to such a position.

        Arguments:
            params, args = as in mychi2.
        Returns:
            the tempered log posterior as a float, or -inf.
        """
        res = mychi2(params, *args)
        if (res > 1.e19 or np.isinf(res) or  np.isnan(res)):
          return -np.inf
        else:
          return -0.5*res
    
    class GaussianStep:
       """Draw points from a Gaussian centered on x with covariance stepsize*cov.

       cov is the covariance of the enclosing min_chi2 call (fixed row and column
       already removed). Draws use NumPy's global random generator, which this
       script never seeds.
       """
       def __init__(self, stepsize=0.2):
           """Store the covariance stepsize*cov.

           Arguments:
               stepsize = factor that multiplies cov; min_chi2 passes T/4.
           """
           self.cov = stepsize*cov
       def __call__(self, x):
           """Return one draw from N(x, stepsize*cov).

           Arguments:
               x = center of the Gaussian [n_dim].
           Returns:
               array [1, n_dim].
           """
           return np.random.multivariate_normal(x, self.cov, size=1)
    
    ndim        = int(x0.shape[0])
    nwalkers    = int(nwalkers)
    nstw        = int(nstw)
    # Temperature ladder; each stage starts its walkers from N(x0, cov*T/4). The
    # main block assumes five temperatures (ntemp = 5).
    temperature = np.array([1.0, 0.25, 0.1, 0.005, 0.001], dtype='float64')
    ntemp       = len(temperature)
    stepsz      = temperature/4.0

    # partial_samples and partial hold the starting point and then the best point
    # of each stage, with their chi2 at T = 1 (the args tuple of this call keeps
    # T = 1).
    partial_samples = [x0]
    partial = [mychi2(x0, *args)]

    for i in range(len(temperature)):
        # Starting walkers of this stage: nwalkers draws from N(x0, cov*T/4).
        # GaussianStep(...)(x0) has shape [1, n_dim]; [0,:] takes its only row.
        x = [] # Initial point
        for j in range(nwalkers):
            x.append(GaussianStep(stepsize=stepsz[i])(x0)[0,:])
        # emcee calls logprob(position, *args), with args = (fixed value, fixed
        # index, temperature of this stage). Moves: DE (differential evolution,
        # 80% of the steps) steps along the difference of two other walkers; DE
        # snooker (20%) uses three other walkers and follows long, curved
        # degeneracies better.
        sampler = emcee.EnsembleSampler(nwalkers, 
                                        ndim, 
                                        logprob, 
                                        args=(args[0], args[1], temperature[i]),
                                        moves=[(emcee.moves.DEMove(), 0.8),
                                               (emcee.moves.DESnookerMove(), 0.2)]) 
        # skip_initial_state_check=True turns off emcee's test that the starting
        # walkers are linearly independent.
        sampler.run_mcmc(np.array(x, dtype='float64'), 
                         nstw, 
                         skip_initial_state_check=True)
        # get_chain(flat=True) stacks every walker at every step into an array
        # [nstw*nwalkers, n_dim]; get_log_prob(flat=True) is the matching tempered
        # log posterior [nstw*nwalkers]. The smallest -log posterior marks the
        # best sample of the stage.
        samples = sampler.get_chain(flat=True, discard=0)
        j = np.argmin(-1.0*np.array(sampler.get_log_prob(flat=True)))
        partial_samples.append(samples[j])
        partial.append(mychi2(samples[j], *args))
        # The next stage starts from this stage's best sample, even when an earlier
        # stage found a smaller chi2; the returned result takes the best of all.
        x0 = copy.deepcopy(samples[j])
        sampler.reset()  
    # The result is the best of the starting point and the stage optima.
    j = np.argmin(np.array(partial))
    result = [partial_samples[j], partial[j]]
    return result
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
def prf(x0, fixed, nstw, nwalkers, cov):
    """Minimize at one grid point: run min_chi2 from x0 converted to float64.

    executor.map calls this function once per grid point, on a worker process.
    As written, the call to min_chi2 raises TypeError (see min_chi2).

    Arguments:
        x0 = starting point [n_sampled]: one row of xf, with the scanned
            parameter at its grid value.
        fixed = index of the scanned parameter.
        nstw = emcee steps per walker at each temperature.
        nwalkers = number of walkers.
        cov = prior covariance [n_sampled, n_sampled].
    Returns:
        [best point without the scanned coordinate, its chi2].
    """
    res =  min_chi2(x0=np.array(x0, dtype='float64'), 
                    cov=cov, 
                    fixed=fixed,
                    nstw=nstw, 
                    nwalkers=nwalkers)
    return res
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# mpi4py.futures: under "python -m mpi4py.futures", rank 0 runs this file as the
# main script and the other ranks wait for tasks. The workers load this file
# under another module name, so the module-level code above (options, model)
# runs on every process, while the block under __main__ runs on rank 0 only.
from mpi4py import MPI
from mpi4py.futures import MPIPoolExecutor

if __name__ == '__main__':
    # 1st: grid over the prior range -------------------------------------------
    # One grid point per MPI process (Get_size() counts rank 0 too), evenly spaced
    # over the central 0.999999 interval of the scanned parameter's prior: its
    # full range for a flat prior.
    executor = MPIPoolExecutor()
    comm = MPI.COMM_WORLD
    numpts = comm.Get_size()

    dim    = model.prior.d()
    index  = args.profile
    bounds = model.prior.bounds(confidence=0.999999)              # Cobaya call
    start  = np.zeros(dim, dtype='float64')
    stop   = np.zeros(dim, dtype='float64')
    for i in range(dim):
      start[i] = bounds[i][0]
      stop[i]  = bounds[i][1]
    param = np.linspace(start = start[index], 
                        stop  = stop[index], 
                        num   = numpts)
    
    # 2nd: print the settings --------------------------------------------------
    # ntemp = 5 must equal the number of temperatures in min_chi2; maxfeval =
    # nstw*ntemp*nwalkers counts the chi2 evaluations of one grid point. Three
    # walkers per sampled parameter (emcee's moves refuse fewer than two).
    names = list(model.parameterization.sampled_params().keys()) # Cobaya Call
    nstw = args.nstw
    ntemp    = 5
    nwalkers = 3*dim
    maxfeval = nstw*ntemp*nwalkers
    print(f"maxfeval={maxfeval}, " 
          f"nstw (evals/Temp/walkers)={nstw}, "
          f"param={names[index]}")
    print(f"profile param values = {param}")

    # 3rd: starting points -----------------------------------------------------
    # xf [numpts, n_sampled]: one random draw from the "ref" distributions, copied
    # to every row, with the scanned column set to the grid value. A row whose
    # chi2 is not finite is redrawn, at most 101 times. The error below fires
    # only when jstop ends at exactly 100 (the 100th redraw succeeded), not when
    # every redraw failed.
    (x0, results) = model.get_valid_point(max_tries=1000, 
                                          ignore_fixed_ref=False,
                                          logposterior_as_dict=True)
    xf = np.tile(x0, (numpts, 1))
    xf[:,index] = param
    for i in range(numpts):
      tmp = chi2(xf[i,:])
      jstop = 0
      while (tmp > 1.e19 or np.isinf(tmp) or  np.isnan(tmp)) and jstop < 101:
        (xf[i,:], results) = model.get_valid_point(max_tries=1000, 
                                                   ignore_fixed_ref=False,
                                                   logposterior_as_dict=True)
        xf[i,index] = param[i]
        tmp = chi2(xf[i,:])
        jstop = jstop + 1
      if jstop == 100:
        raise RuntimeError(f"Can't find an initial point with finite"
                           f" likelihood at fixed param = {param[i]}")
    
    # 4th: one minimization per grid point -------------------------------------
    # functools.partial(prf, fixed=..., nstw=..., nwalkers=..., cov=...) is prf
    # with every argument except x0 filled in. executor.map calls it once per row
    # of xf on the worker processes and yields the results in row order; res is
    # the object array [numpts, 2] of [best point, chi2] pairs.
    cov = model.prior.covmat(ignore_external=False) # cov from prior
    res = np.array(list(executor.map(functools.partial(prf, 
                                                       fixed=index,
                                                       nstw=nstw, 
                                                       nwalkers=nwalkers,
                                                       cov=cov), xf)),dtype="object")
    # The list comprehension puts each grid value p back at the scanned index of
    # its best point (res[:,0]), giving xf [numpts, n_sampled] again.
    xf = np.array([np.insert(row,index,p) for row, p in zip(res[:,0], param)], dtype='float64')
    chi2res = np.array(res[:,1], dtype='float64')
    
    # 5th: append the chi2 of each likelihood and -2 log prior -----------------
    xf = np.column_stack((xf, 
                          np.array([chi2v2(d) for d in xf], dtype='float64')))
    
    # 6th: write <root>chains/<outroot>.<name>.txt -----------------------------
    # Columns: grid value, chi2, then the rows of xf; np.c_[a, b] stacks the
    # one-dimensional arrays a and b as the columns of a two-column table.
    os.makedirs(os.path.dirname(f"{args.root}chains/"), exist_ok=True)
    hd = [names[index], "chi2"] + names
    hd = hd + list(model.info()['likelihood'].keys()) + ["prior"]
    os.makedirs(os.path.dirname(f"{args.root}chains/"),exist_ok=True)
    np.savetxt(f"{args.root}chains/{args.outroot}.{names[index]}.txt",
               np.concatenate([np.c_[param, chi2res], xf],axis=1),
               fmt="%.6e",
               header=f"maxfeval={maxfeval}, param={names[index]}\n"+' '.join(hd),
               comments="# ")
    # 7th: release the worker processes ----------------------------------------
    executor.shutdown()
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------