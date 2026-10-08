r"""Compute the profile of one parameter of the example CAMB model.

A profile fixes one parameter at each value of a grid and, at each value,
minimizes chi2 = -2 log posterior = -2 (log prior + log likelihood) over all the
other parameters. The curve of minimized chi2 against the fixed value is the
profile. For one parameter, the values where it rises by Delta chi2 = 1 and 4
above its minimum bound the 68% and 95% confidence intervals of a Gaussian case;
unlike a marginalized posterior, a profile does not depend on the volume of the
parameters it minimizes over. The priors stay in chi2: flat priors only add a
constant, but the Gaussian prior of A_planck enters the minimized quantity, so
this is a profile of the posterior rather than of the likelihood alone.

The model is the YAML text in yaml_string:
- likelihoods: Planck 2018 high-l TT, TE, EE (plik-lite) and low-l TT; Planck
  2020 low-l EE (LoLLiPoP); DES-Y5 supernovae; DESI DR2 BAO; ACT DR6 CMB lensing;
- cosmology: LCDM with the neutrino mass fixed at 0.06 eV, computed by CAMB, the
  Boltzmann code that solves for the CMB and matter power spectra (here with
  CosmoRec recombination and HMcode-2020 nonlinear corrections);
- sampled parameters (7), in sampled order: logA = ln(10^10 A_s), ns,
  theta_MC_100 = 100 theta_MC (CosmoMC's approximation to the acoustic angle),
  omegabh2, omegach2, tau, and the Planck calibration A_planck, which the Planck
  likelihoods add with a Gaussian prior (mean 1, standard deviation 0.0025).

Each minimization is annealed emcee sampling, as in EXAMPLE_MINIMIZE1.py: a set
of points called walkers samples the posterior raised to the power 1/T while the
temperature T falls along a ladder, so the walkers crowd around the minimum
(min_chi2 gives the details). The module docstring of
external_modules/code/cosmolike_core/cocoa_hybrid_sampling.py explains the same
annealing scheme (ladders, cov*T/3 starting clouds, DE and snooker moves) and the
reasons for its choices.

Steps, numbered as in the main block: (1) covariance; (2) global minimum, read
from --minfile or computed; (3) grid of the profiled parameter, centered on the
minimum; (4) print the grid; (5) result tables; (6, 7) minimizations at the grid
points, walking outward from the minimum so that each starts from the optimum of
its neighbor; (8) derived parameters and chi2 values; (9) output file.

Inputs (command line; paths are relative to the Cocoa folder):
  --root     folder, ending in "/", that receives chains/ and holds --cov
             (default ./projects/example/).
  --outroot  basename of the output file (default test.dat).
  --profile  zero-based index, in sampled order, of the profiled parameter
             (default 1, ns).
  --numpts   number of grid points besides the minimum (default 20); an odd
             value is lowered by one, so the grid is symmetric.
  --factor   half-width of the grid, in standard deviations sqrt(diag(cov)) of
             the profiled parameter (default 3; at most 1 without --cov).
  --nstw     emcee steps per walker at each temperature (default 200; the global
             minimization uses 5*nstw/4).
  --cov      covariance file relative to --root, such as the .covmat of a chain
             of this model; its first n_sampled rows and columns must be the
             sampled parameters in sampled order (the header is not read).
             Without it, the prior covariance is used.
  --minfile  one-row file, relative to the Cocoa folder (not to --root), whose
             first n_sampled columns are the global minimum and whose last column
             is its chi2, computed with this same YAML; this is the layout that
             EXAMPLE_MINIMIZE1.py writes, but that script uses another low-l EE
             likelihood, so its output fails the consistency check here. Without
             --minfile the minimum is computed first, which takes long.

Output: <root>chains/<outroot>.<name>.txt, with <name> the profiled parameter,
one row per grid point: the grid value, chi2, the sampled parameters of the
minimization, H0 [km/s/Mpc], Omega_m, r_drag [Mpc], the chi2 of each likelihood,
and -2 log prior. The header records nstw and names the columns.

Cost: each grid point runs 4 * nstw * nwalkers CAMB evaluations, with nwalkers =
max(3 * n_sampled, number of MPI processes). No random seed is set.

Run from the Cocoa folder after "source start_cocoa.sh". Rank 0 coordinates and
the other MPI processes evaluate the walkers, each with OMP_NUM_THREADS OpenMP
threads for CAMB:

    export OMP_NUM_THREADS=4
    mpirun -n 22 --oversubscribe python ./projects/example/EXAMPLE_PROFILE1.py \
        --root ./projects/example/ --outroot EXAMPLE_PROFILE1 --profile 1 \
        --numpts 10 --nstw 200
"""

import warnings, os, psutil
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
from schwimmbad import MPIPool
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# True once the process has pinned itself to its cores: logprob reads it so that
# each worker process calls enforce_affinity at its first evaluation only.
_affinity_set = False

def enforce_affinity():
    """Pin this MPI process to its own block of CPU cores.

    Every MPI rank runs OMP_NUM_THREADS OpenMP threads (CAMB's). Rank r is pinned
    to cores r*OMP_NUM_THREADS through (r+1)*OMP_NUM_THREADS - 1, so the ranks do
    not compete for cores. schwimmbad's MPIPool runs one Python process per MPI
    rank, so each process sets its own affinity with psutil instead of relying on
    the binding options of mpirun. The rank is read from OMPI_COMM_WORLD_RANK,
    which Open MPI exports (0 when absent). The core numbers assume that all ranks
    share one node whose cores are numbered consecutively. Where psutil cannot
    set an affinity (macOS has no such call), the error is printed and the
    process runs unpinned.

    Side effects:
        changes the CPU affinity of the calling process.
    """
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
# Command-line options; the module docstring gives their meaning. Every option
# uses nargs='?' with const=1: written without a value, an option takes the value
# 1 (for --root that is the integer 1, not a path), so always give a value.
parser = argparse.ArgumentParser(prog='EXAMPLE_EMUL_PROFILE1')
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
parser.add_argument("--factor",
                    dest="factor",
                    help="Factor that set the bounds (multiple of cov matrix)",
                    type=float,
                    nargs='?',
                    const=1.0,
                    default=3.0)
parser.add_argument("--numpts",
                    dest="numpts",
                    help="Number of Points to Compute Minimum",
                    type=int,
                    nargs='?',
                    const=1,
                    default=20)
parser.add_argument("--minfile",
                    dest="minfile",
                    help="Minimization Result",
                    nargs='?',
                    const=1)
parser.add_argument("--cov",
                    dest="cov",
                    help="Chain Covariance Matrix",
                    nargs='?',
                    const=1,
                    default=None)
# parse_known_args, unlike parse_args, does not stop at an argument it does not
# know: it returns such arguments in `unknown`, which nothing reads. A misspelled
# option is therefore ignored without a message, and its default applies.
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
# script; "drop: true" keeps a parameter away from CAMB, which receives As and
# cosmomc_theta computed from it); and the CAMB settings. Paths such as
# ./external_modules are relative to the Cocoa folder, the folder to run from.
yaml_string=r"""
likelihood:
  planck_2018_highl_plik.TTTEEE_lite: 
    path: ./external_modules/
    clik_file: plc_3.0/hi_l/plik_lite/plik_lite_v22_TTTEEE.clik
  planck_2018_lowl.TT: 
    path: ./external_modules
  planck_2020_lollipop.lowlE:
    data_folder: planck/lollipop
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
      min: 2.7
      max: 3.4
    ref:
      dist: norm
      loc: 3.04
      scale: 0.025
    proposal: 0.025
    latex: \log(10^{10} A_\mathrm{s}
    drop: true
  ns:
    prior:
      min: 0.93
      max: 1.01
    ref:
      dist: norm
      loc: 0.96
      scale: 0.0075
    proposal: 0.0075
    latex: n_\mathrm{s}
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
      max: 0.03
    ref:
      dist: norm
      loc: 0.022383
      scale: 0.005
    proposal: 0.005
    latex: \Omega_\mathrm{b} h^2
  omegach2:
    prior:
      min: 0.08
      max: 0.16
    ref:
      dist: norm
      loc: 0.12011
      scale: 0.01
    proposal: 0.01
    latex: \Omega_\mathrm{c} h^2
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
      halofit_version: mead2020
      lmax: 4000
      lens_margin: 1250
      AccuracyBoost: 1.05
      lens_potential_accuracy: 4
      lens_k_eta_reference: 18000
      nonlinear: NonLinear_both
      recombination_model: CosmoRec
      Accuracy.AccurateBB: True
"""
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Module-level code runs on every MPI process, so each process builds its own
# Cobaya model (CAMB, likelihoods and data) and evaluates the walkers it receives
# with it.
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
    res1 = model.logprior(point,
                          make_finite=False)
    if np.isinf(res1) or np.any(np.isnan(res1)):
      return 1.e20
    # cached=False recomputes the theory even if Cobaya holds a result for these
    # parameters; return_derived=False skips the derived parameters.
    res2 = model.loglike(point,
                         make_finite=False,
                         cached=False,
                         return_derived=False)
    if np.isinf(res2) or np.isnan(res2):
      return 1.e20
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
             nstw=200,
             nwalkers=5,
             pool=None):
    """Minimize chi2 by annealed emcee sampling; return the best point found.

    At temperature T the walkers sample the posterior raised to the power 1/T,
    whose log is -chi2/(2T). T = 1 samples the posterior itself. For a nearly
    Gaussian posterior the walkers then sit about n_dim*T above the minimum in
    chi2: about 0.007 at T = 0.001 with seven parameters. Each stage draws
    nwalkers starting points from a Gaussian centered on x0 with covariance
    cov*T/3, runs nstw emcee steps, and passes its best sample to the next stage
    as the new x0. The function returns the best of the starting point and the
    stage optima, so the result is never worse than x0.

    Ladders: a free minimization (fixed = -1) runs T = 1.0, 0.25, 0.1, 0.005,
    0.001. A profile point (fixed >= 0) runs 0.3, 0.1, 0.005, 0.001: its x0 is
    the optimum of a neighboring grid point, so the basin is already found and
    the stages only re-polish after one coordinate moved.

    Why cov*T/3: posterior^(1/T) of a Gaussian posterior with covariance C has
    covariance T*C. When cov is close to the posterior covariance, the starting
    cloud thus has one third of the variance the stage samples (0.58 times its
    standard deviation) and starts inside the region the stage refines.

    Arguments:
        x0 = starting point [n_sampled], in the model's sampled order.
        cov = covariance [n_sampled, n_sampled] in the same order; the row and the
            column of a fixed parameter are removed before walkers are drawn.
        fixed = index of the parameter held at x0[fixed], or -1 for none.
        nstw = emcee steps per walker at each temperature.
        nwalkers = number of walkers.
        pool = schwimmbad MPIPool whose worker processes evaluate the walkers, or
            None to evaluate them in this process.
    Returns:
        best point [n_dim], where n_dim = n_sampled, or n_sampled - 1 when a
        parameter is fixed (the fixed value is not included).
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

        The first call in each process (with a pool, the first evaluation on each
        worker) also pins the process to its cores with enforce_affinity and
        prints the time of that evaluation with the process's MPI rank. A chi2
        above 1e19 (the rejection value of chi2) becomes -inf, and emcee never
        accepts a move to such a position.

        Arguments:
            params, args = as in mychi2.
        Returns:
            the tempered log posterior as a float, or -inf.
        Side effects:
            sets the module-level flag _affinity_set on the first call.
        """
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
       """Draw points from a Gaussian centered on x with covariance stepsize*cov.

       cov is the covariance of the enclosing min_chi2 call (fixed row and column
       already removed). Draws use NumPy's global random generator, which this
       script never seeds.
       """
       def __init__(self, stepsize=0.2):
           """Store the covariance stepsize*cov.

           Arguments:
               stepsize = factor that multiplies cov; min_chi2 passes T/3.
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
    # Temperature ladders (free minimization, profile point); each stage starts
    # its walkers from N(x0, cov*T/3), see the docstring for both choices.
    if fixed == -1:
      temperature = np.array([1.0, 0.25, 0.1, 0.005, 0.001], dtype='float64')
    else:
      temperature = np.array([0.3, 0.1, 0.005, 0.001], dtype='float64')
    stepsz      = temperature/3.0

    # partial_samples and partial hold the starting point and then the best point
    # of each stage, with their chi2 at T = 1 (the args tuple of this call keeps
    # T = 1).
    partial_samples = [x0]
    partial = [mychi2(x0, *args)]

    for i in range(len(temperature)):
        # Starting walkers of this stage: nwalkers draws from N(x0, cov*T/3).
        # GaussianStep(...)(x0) has shape [1, n_dim]; [0,:] takes its only row.
        x = [] # Initial point
        for j in range(nwalkers):
            x.append(GaussianStep(stepsize=stepsz[i])(x0)[0,:]) 
        # emcee calls logprob(position, *args), with args = (fixed value, fixed
        # index, temperature of this stage); pool spreads the evaluations over the
        # MPI workers. Moves: DE (differential evolution, 80% of the steps) steps
        # along the difference of two other walkers; DE snooker (20%) uses three
        # other walkers and follows long, curved degeneracies better.
        sampler = emcee.EnsembleSampler(nwalkers=nwalkers, 
                                        ndim=ndim, 
                                        log_prob_fn=logprob, 
                                        args=(args[0], args[1], temperature[i]),
                                        moves=[(emcee.moves.DEMove(), 0.8),
                                               (emcee.moves.DESnookerMove(), 0.2)],
                                        pool=pool)
        # skip_initial_state_check=True turns off emcee's test that the starting
        # walkers are linearly independent.
        sampler.run_mcmc(np.array(x,dtype='float64'), 
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
        # stage found a smaller chi2; the returned point is the best of all.
        x0 = copy.deepcopy(samples[j])
        sampler.reset()
    # The result is the best of the starting point and the stage optima.
    j = np.argmin(np.array(partial))
    return partial_samples[j]
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
def prf(x0, nstw, cov, fixed=-1, nwalkers=5, pool=None):
    """Run min_chi2 from x0, converted to a float64 array; return its result.

    Arguments:
        x0, nstw, cov, fixed, nwalkers, pool = as in min_chi2.
    Returns:
        the best point, as min_chi2 returns it (without the fixed coordinate).
    """
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
# Every MPI process enters this block. With schwimmbad 0.4.2, the version the
# Cocoa environment pins, MPIPool() keeps every process except rank 0 inside its
# constructor: those workers evaluate the walker positions that rank 0 sends, and
# they exit when rank 0 closes the pool at the end of the with block. Only rank 0
# runs the lines below; the is_master() test repeats that guard.
if __name__ == '__main__':
    with MPIPool() as pool:
        if not pool.is_master():
            pool.wait()
            sys.exit(0)
        # Rank 0 pins itself here; each worker pins itself at its first
        # evaluation (logprob).
        enforce_affinity() # enforce affinity (so Hybrid MPI-OpenMP works)!

        # dim = model.prior.d() = number of sampled parameters. Walkers: three per
        # sampled parameter (emcee's moves refuse fewer than two per parameter),
        # and at least as many as MPI processes; Get_size() counts every process,
        # rank 0 included.
        dim      = model.prior.d()     
        nwalkers = max(3*dim, pool.comm.Get_size())
        nstw = args.nstw

        # 1st: covariance ------------------------------------------------------
        # Without --cov: the prior covariance (diagonal, each entry the variance of
        # one prior, (max - min)^2/12 for a flat prior), and the grid half-width is
        # capped at one prior standard deviation (factor <= 1). With --cov: the
        # first n_sampled rows and columns of the file, which must be the sampled
        # parameters in sampled order. sigma = standard deviation of each one.
        if args.cov is None:
          cov = model.prior.covmat(ignore_external=False) # cov from prior
          factor = min(1.0, args.factor)
        else:
          cov = np.loadtxt(args.root+args.cov)[0:model.prior.d(),0:model.prior.d()]
          factor = args.factor
        sigma = np.sqrt(np.diag(cov))

        # 2nd: global minimum --------------------------------------------------
        # With --minfile: the first n_sampled columns of its row are x0 and the
        # last column is chi20, its chi2. Otherwise: a free minimization from a
        # random draw of the "ref" distributions (max_tries bounds the number of
        # draws), with 5/4 of nstw steps per walker at each temperature.
        if args.minfile is not None: # load minimum from running MCMC
          x0 = np.loadtxt(args.minfile)
          chi20 = x0[-1]
          x0 = x0[0:model.prior.d()]
        else: # Compute the minimum (slow)
          (x0, results) = model.get_valid_point(max_tries=1000, 
                                     ignore_fixed_ref=False,
                                     logposterior_as_dict=True)
          res = np.array(list(prf(x0=x0, 
                                  nstw=int(5.*nstw/4.), 
                                  nwalkers=nwalkers,
                                  pool=pool,
                                  cov=cov,
                                  fixed=-1)), dtype="object")
          x0 = np.array(res, dtype='float64')[0:model.prior.d()]
          chi20 = chi2(x0)
          print(f"Global Min: params = {x0}, and chi2 = {chi20}")

        # The minimum must reproduce its own chi2 within 0.02; a larger difference
        # means --minfile came from another model (another YAML, another
        # parameter order). A minimum computed above passes trivially.
        if (abs(chi2(x0)-chi20)>0.02):
          raise ValueError("Inconsistency Min and Profile setups")

        # 3rd: grid of the profiled parameter ----------------------------------
        # Range from x0 - factor*sigma to x0 + factor*sigma for every parameter;
        # the loop below clips it to the prior.
        start = np.zeros(model.prior.d(), dtype='float64')
        stop  = np.zeros(model.prior.d(), dtype='float64')
        start = x0 - factor*sigma
        stop  = x0 + factor*sigma
        
        # Clip to the central 0.999999 interval of each prior: the full range of a
        # flat prior, about 4.9 standard deviations around the mean of a Gaussian
        # prior (A_planck).
        bounds0 = model.prior.bounds(confidence=0.999999)
        for i in range(model.prior.d()):
            if (start[i] < bounds0[i][0]):
              start[i] = bounds0[i][0]
            if (stop[i] > bounds0[i][1]):
              stop[i] = bounds0[i][1]

        # Only the clipped width of the profiled parameter is kept: the grid is
        # centered on the minimum with half that width. After a clip on one side
        # only, the grid therefore reaches past the prior on that side, and the
        # grid points there get chi2 = 1e20.
        half_range = (stop[args.profile] - start[args.profile]) / 2.0
       
        # Grid: an even number of points from x0 - half_range to x0 + half_range,
        # plus the minimum's own value inserted in the middle, so the minimum is
        # the center row. Example: --numpts 10 (or 11, lowered to 10), x0 = 0.966,
        # half_range = 0.012 gives 0.954, 0.95667, ..., 0.96467, then 0.966 at
        # index 5, then 0.96733, ..., 0.978: 11 sorted points.
        numpts = args.numpts-1 if args.numpts%2 == 1 else args.numpts 
      
        param  = np.linspace(start = x0[args.profile] - half_range,
                             stop  = x0[args.profile] + half_range,
                             num = numpts)
        numpts=numpts+1
        param = np.insert(param, numpts//2, x0[args.profile])
        
        # 4th: print the grid --------------------------------------------------
        names = list(model.parameterization.sampled_params().keys()) # Cobaya Call
        print(f"nstw (evals/Temp/walkers)={args.nstw}, "
              f" param={names[args.profile]}\n"
              f"profile param values = {param}")
        
        # 5th: result tables ---------------------------------------------------
        # xf [numpts, n_sampled]: one row per grid point, x0 with the profiled
        # column set to the grid value; chi2res [numpts]: the minimized chi2. The
        # center row is the minimum itself.
        xf = np.tile(x0, (numpts, 1))
        xf[:,args.profile] = param

        chi2res = np.zeros(numpts)  
        chi2res[numpts//2] = chi20
        
        # 6th: minimize from the center to the last grid point -----------------
        # Each grid point starts from the optimum of the previous one (tmp), which
        # is close to its own optimum. prf returns the point without the profiled
        # coordinate; np.insert puts the grid value back.
        tmp = np.array(xf[numpts//2,:], dtype='float64')
        for i in range(numpts//2+1,numpts): 
            tmp[args.profile] = param[i]
            res = prf(tmp, 
                      fixed=args.profile,
                      nstw=int(nstw), 
                      nwalkers=nwalkers,
                      pool=pool,
                      cov=cov)
            xf[i,:] = np.insert(res, args.profile, param[i])
            tmp = np.array(xf[i,:],dtype='float64')
            chi2res[i] = chi2(xf[i,:])
            print(f"Partial ({i+1}/{numpts}): params={tmp}, and chi2={chi2res[i]}")
        
        # 7th: minimize from the center to the first grid point ----------------
        tmp = np.array(xf[numpts//2,:], dtype='float64')
        for i in range(numpts//2-1, -1, -1):
            tmp[args.profile] = param[i]
            res = prf(tmp, 
                      fixed=args.profile,
                      nstw=int(nstw), 
                      nwalkers=nwalkers,
                      pool=pool,
                      cov=cov)
            xf[i,:] = np.insert(res, args.profile, param[i])
            tmp = np.array(xf[i,:],dtype='float64')
            chi2res[i] = chi2(xf[i,:])
            print(f"Partial ({i+1}/{numpts}): params={tmp}, and chi2={chi2res[i]}")
        
        # 8th: append derived parameters and chi2 values -----------------------
        # H0 [km/s/Mpc], Omega_m and r_drag [Mpc] from CAMB, through the derived
        # parameters of Cobaya's logposterior; then the chi2 of each likelihood
        # and -2 log prior, for every grid point.
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
                              np.array([chi2v2(d) for d in xf], dtype='float64')))

        # 9th: write <root>chains/<outroot>.<name>.txt -------------------------
        # Columns: grid value, chi2, then the rows of xf; np.c_[a, b] stacks the
        # one-dimensional arrays a and b as the columns of a two-column table.
        os.makedirs(os.path.dirname(f"{args.root}chains/"), exist_ok=True)
        hd = [names[args.profile], "chi2"] + names + ['H0', 'omegam', 'rdrag']
        hd = hd + list(model.info()['likelihood'].keys()) + ["prior"]
        np.savetxt(f"{args.root}chains/{args.outroot}.{names[args.profile]}.txt",
                   np.concatenate([np.c_[param, chi2res], xf], axis=1),
                   fmt="%.6e",
                   header=f"nstw (evals/Temp/walker)={args.nstw}, param={names[args.profile]}\n"+' '.join(hd),
                   comments="# ")
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------