r"""Sample the example emulator model with emcee; write the chain in GetDist format.

Markov chain Monte Carlo (MCMC) draws samples whose density follows the
posterior, so histograms of the samples estimate the parameter constraints.
emcee is an ensemble MCMC sampler: a set of points called walkers moves through
parameter space together, and each proposed step of a walker is built from the
positions of other walkers. Here 80% of the steps are differential-evolution (DE)
moves, along the difference of two other walkers, and 20% are DE snooker moves,
built from three other walkers, which follow long, curved degeneracies better.
Unlike the minimizers of this project, this script samples the posterior itself
(temperature 1); scripts/EXAMPLE_PLOT_COMPARE_CHAINS.py overlays the result on
the Metropolis-Hastings, Nautilus and PolyChord chains of the same model.

The model is the YAML in yaml_string: the likelihoods and parameters of
EXAMPLE_EMUL_MINIMIZE1.py (Planck 2018 high-l TT, TE, EE plik-lite, low-l TT and
low-l EE SRoll2; DES-Y5 supernovae; DESI DR2 BAO; ACT DR6 lensing; LCDM with the
neutrino mass fixed at 0.06 eV), except that the ACT DR6 lensing likelihood here
has lens_only: True (False there). Emulators (neural networks and Gaussian
processes trained on Boltzmann-code outputs) replace CAMB. Seven parameters are
sampled: logA, ns, thetastar = 100 theta_*, omegabh2, omegach2, tau, and the
Planck calibration A_planck (Gaussian prior, mean 1, standard deviation 0.0025).

Convergence: emcee estimates the integrated autocorrelation time tau of each
parameter, the number of steps after which a walker no longer remembers where
it was. The chain is trustworthy when each walker runs about 50 tau steps; with
tau of order 200 that means --maxfeval above 10000 * nwalkers. The first
5*max(tau) steps are dropped as burn-in, and every int(min(tau)/2)-th step is
kept. The header of the .1.txt file prints tau.

Inputs (command line; paths are relative to the Cocoa folder):
  --root      folder, ending in "/", whose chains/ subfolder receives the output;
              chains/ must exist before the run, because emcee creates its
              checkpoint there before the script creates the folder (default
              ./projects/example/).
  --outroot   basename of the output files (default test.dat).
  --maxfeval  total number of likelihood evaluations (default 5000); each walker
              runs maxfeval/nwalkers steps, with nwalkers = max(3 * n_sampled,
              number of MPI processes).
  --progress  shows emcee's progress bar when true. type=bool turns any
              non-empty text into True, even "False"; the option without a value
              gives None, which shows no bar.

Output files in <root>chains/:
  <outroot>.h5        emcee's HDF5 checkpoint with every step of every walker.
                      If it already exists with the same number of walkers and
                      parameters, the run appends to it, starting again from new
                      random points; another walker count stops with ValueError.
  <outroot>.1.txt     the chain in GetDist's text format: weight (1), the log
                      posterior, the sampled parameters, and chi2 = -2 log
                      posterior. GetDist reads the second column as -log
                      posterior, so its likelihood statistics (such as the
                      best-fit sample) come out with the wrong sign; the
                      marginalized constraints use only the weights.
  <outroot>.ranges    the 0.999999 prior interval of each sampled parameter.
  <outroot>.paramnames  name and LaTeX label of each column after the second.
  <outroot>.covmat    the covariance GetDist computes from the chain, chi2
                      included, with the names in the header.

Run from the Cocoa folder after "source start_cocoa.sh". Rank 0 coordinates and
the other MPI processes evaluate the walkers:

    mpirun -n 12 --oversubscribe python ./projects/example/EXAMPLE_EMUL_EMCEE1.py \
        --root ./projects/example/ --outroot EXAMPLE_EMUL_EMCEE1 --maxfeval 80000
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
from emcee.autocorr import AutocorrError
from cobaya.yaml import yaml_load
from cobaya.model import get_model
from getdist import IniFile
from getdist import loadMCSamples
from schwimmbad import MPIPool
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Command-line options; the module docstring gives their meaning. Options with
# nargs='?' and const=1 take the value 1 when written without a value (for
# --root that is the integer 1, not a path), so always give a value.
parser = argparse.ArgumentParser(prog='EXAMPLE_EMUL_EMCEE')

parser.add_argument("--maxfeval",
                    dest="maxfeval",
                    help="Minimizer: maximum number of likelihood evaluations",
                    type=int,
                    nargs='?',
                    const=1,
                    default=5000)
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
parser.add_argument("--progress",
                    dest="progress",
                    help="Show Emcee Progress",
                    nargs='?',
                    type=bool,
                    default=False)
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
# parameter has a prior; "ref" is the distribution that the starting walkers are
# drawn from; "proposal" is read by Cobaya's own MCMC sampler, not by this
# script); and the theory block, where the emulators replace CAMB. Paths such as
# ./external_modules are relative to the Cocoa folder, the folder to run from.
yaml_string = r"""
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
# Module-level code runs on every MPI process, so each process builds its own
# Cobaya model (theory, likelihoods and data) and evaluates the walkers it receives
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
    likelihood, gets it, and logprob in chain maps any chi2 above 1e19 to -inf.

    Arguments:
        p = sampled parameter values in the model's sampled order: a list or an
            array [n_sampled], or a dict whose values are in that order (emcee
            passes a dict here, see chain).
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
      return 1.e20
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

    Not called by this script; it matches the function of the other scripts.

    Arguments:
        p = sampled parameter values, as for chi2.
    Returns:
        array [n_likelihoods + 1]: -2 log likelihood of each likelihood, in the
        order of model.info()["likelihood"], then -2 log prior.
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
def chain(x0,
          ndim,
          nwalkers,
          cov,
          names,
          maxfeval=3000, 
          pool=None,
          checkpoint=None):    
    """Run emcee from x0 and return the chain after burn-in and thinning.

    The sampler stores every step in the HDF5 file checkpoint (emcee's
    HDFBackend). After maxfeval steps per walker, emcee estimates the integrated
    autocorrelation time tau of each parameter: the number of steps after which
    a walker's position no longer depends on where it was. The first 5*max(tau)
    steps are dropped as burn-in, the stretch where the walkers still remember
    their random starting points, and every int(min(tau)/2)-th step is kept
    (thinning, which lowers the correlation between kept rows).

    Arguments:
        x0 = starting walkers [nwalkers, ndim].
        ndim = number of sampled parameters.
        nwalkers = number of walkers.
        cov = prior covariance; not used.
        names = sampled-parameter names, in sampled order.
        maxfeval = steps per walker; the main block passes --maxfeval/nwalkers.
        pool = schwimmbad MPIPool whose worker processes evaluate the walkers, or
            None to evaluate them in this process.
        checkpoint = path of the HDF5 file.
    Returns:
        [table [n_kept, ndim + 3], tau [ndim]]: the columns of the table are the
        weight (1), the log posterior, the sampled parameters, and chi2 = -2 log
        posterior; tau is in steps.
    Side effects:
        creates or extends the HDF5 checkpoint; prints tau.
    """
    def logprob(params, *args):
        """Return the log posterior -chi2/2 that emcee samples, or -inf.

        emcee passes params as a dict {name: value}, because the sampler is built
        with parameter_names; chi2 turns it back into the values in sampled
        order. A chi2 above 1e19 (the rejection value of chi2) becomes -inf, and
        emcee never accepts a move to such a position.

        Arguments:
            params = dict {sampled-parameter name: value} of one walker.
            args = not used.
        Returns:
            the log posterior as a float, or -inf.
        """
        res = chi2(params)
        if (res > 1.e19 or np.isinf(res) or  np.isnan(res)):
          return -np.inf
        else:
          return -0.5*res

    # The backend writes every step of every walker to the file as the run
    # proceeds. A file that already holds a chain of the same shape is extended,
    # not replaced.
    backend = emcee.backends.HDFBackend(checkpoint)

    # Moves: DE (differential evolution, 80% of the steps) steps along the
    # difference of two other walkers; DE snooker (20%) uses three other walkers.
    # pool spreads the evaluations over the MPI workers. parameter_names makes
    # emcee pass logprob a dict {name: value} instead of an array.
    sampler = emcee.EnsembleSampler(nwalkers=nwalkers, 
                                    ndim=ndim, 
                                    log_prob_fn=logprob, 
                                    parameter_names=names,
                                    moves=[(emcee.moves.DEMove(), 0.8),
                                           (emcee.moves.DESnookerMove(), 0.2)],
                                    pool=pool,
                                    backend=backend)
    # maxfeval steps per walker. skip_initial_state_check=True turns off emcee's
    # test that the starting walkers are linearly independent.
    sampler.run_mcmc(x0, 
                     maxfeval, 
                     skip_initial_state_check=True, 
                     progress=args.progress)
    
    # tau [ndim], in steps. quiet=True makes emcee warn, instead of raising
    # AutocorrError, when the chain is shorter than 50 tau, the length emcee
    # needs to trust its own estimate; has_walkers=True averages the
    # autocorrelation over the walkers.
    tau = sampler.get_autocorr_time(quiet=True, has_walkers=True)
    print(f"Partial Result: tau = {tau}, nwalkers={nwalkers}")

    # Burn-in and thinning, in steps. thin is 0 when min(tau) < 2, and get_chain
    # then stops with ValueError (a slice step cannot be zero).
    burn_in = int(5*np.max(tau))
    thin    = int(0.5 * np.min(tau))
    xf      = sampler.get_chain(flat=True, discard=burn_in, thin=thin)
    lnpf    = sampler.get_log_prob(flat=True, discard=burn_in, thin=thin)
    # Every kept sample has weight 1: emcee samples are unweighted.
    weights = np.ones((len(xf),1), dtype='float64')
    local_chi2    = -2*lnpf
    
    return [np.concatenate([weights,
                           lnpf[:,None], 
                           xf, 
                           local_chi2[:,None]], axis=1), 
            tau]

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
        
        dim      = model.prior.d()                                      # Cobaya call
        bounds   = model.prior.bounds(confidence=0.999999)              # Cobaya call
        names    = list(model.parameterization.sampled_params().keys()) # Cobaya Call
        # Walkers: three per sampled parameter (emcee's moves refuse fewer than two
        # per parameter), and at least as many as MPI processes; Get_size() counts
        # every process, rank 0 included. maxevals = steps per walker.
        nwalkers = max(3*dim,pool.comm.Get_size())
        maxevals = int(args.maxfeval/(nwalkers))
        print(f"\n\n\n"
              f"maxfeval={args.maxfeval}, "
              f"nwalkers={nwalkers}, "
              f"maxfeval per walker = {maxevals}"
              f"\n\n\n")
        # starting walkers -----------------------------------------------------
        # nwalkers independent draws from the "ref" distributions of the YAML,
        # each with a finite posterior (max_tries bounds the number of draws).
        x0 = [] # Initial point x0
        for j in range(nwalkers):
          (tmp_x0, tmp) = model.get_valid_point(max_tries=10000, 
                                                ignore_fixed_ref=False,
                                                logposterior_as_dict=True)
          x0.append(tmp_x0[0:dim])
        x0 = np.array(x0, dtype='float64')
        
        # prior covariance (passed to chain, which does not use it) ------------
        cov = model.prior.covmat(ignore_external=False) # cov from prior

        # HDF5 checkpoint: <root>chains/<outroot>.h5 ---------------------------
        # emcee creates it inside <root>chains/, which must already exist; the
        # folder is created only after the run.
        checkpoint = f"{args.root}chains/{args.outroot}.h5"

        # run emcee ------------------------------------------------------------
        res = chain(x0=np.array(x0, dtype='float64'),
                    ndim=dim,
                    nwalkers=nwalkers,
                    cov=cov, 
                    names=names,
                    maxfeval=maxevals,
                    pool=pool,
                    checkpoint=checkpoint)

        # write <outroot>.1.txt, the chain in GetDist's text format ------------
        # One row per kept sample: weight, log posterior, sampled parameters, chi2.
        # GetDist reads the second column as -log posterior; this file holds +log
        # posterior there, so GetDist's likelihood statistics have the wrong sign
        # (its marginalized constraints use only the weights). The header line
        # records nwalkers, maxfeval and the array tau.
        os.makedirs(os.path.dirname(f"{args.root}chains/"),exist_ok=True)
        hd=f"nwalkers={nwalkers}, maxfeval={args.maxfeval}, max tau={res[1]}\n"
        np.savetxt(f"{args.root}chains/{args.outroot}.1.txt",
                   res[0],
                   fmt="%.7e",
                   header=hd + ' '.join(names),
                   comments="# ")
        # write <outroot>.ranges -----------------------------------------------
        # One line per sampled parameter: name, lower and upper end of its 0.999999
        # prior interval. GetDist reads this file to know where each prior stops,
        # so its smoothed densities do not spread past a prior edge. The first
        # line, which starts with #, is a comment listing the chain columns.
        # rows is a list of (name, low, high) tuples; writelines receives a
        # generator expression that formats one text line per tuple.
        hd = ["weights","lnp"] + names + ["chi2*"]
        rows = [(str(n),float(l),float(h)) for n,l,h in zip(names, bounds[:,0], bounds[:,1])]
        with open(f"{args.root}chains/{args.outroot}.ranges", "w") as f: 
          f.write(f"# {' '.join(hd)}\n")
          f.writelines(f"{n} {l:.5e} {h:.5e}\n" for n, l, h in rows)

        # write <outroot>.paramnames -------------------------------------------
        # One line per column after the second: the name and its LaTeX label from
        # the YAML. In GetDist's convention the trailing * of chi2* marks a derived
        # parameter. names gains chi2* here, so the covariance header includes it.
        param_info = model.info()['params']
        latex  = [param_info[x]['latex'] for x in names]
        names.append("chi2*")
        latex.append("\\chi^2")
        np.savetxt(f"{args.root}chains/{args.outroot}.paramnames", 
                   np.column_stack((names,latex)),
                   fmt="%s")
    
        # write <outroot>.covmat -----------------------------------------------
        # GetDist reloads the chain just written (ignore_rows 0: burn-in is already
        # removed) and computes the covariance of the parameter columns, chi2
        # included.
        samples = loadMCSamples(f"{args.root}chains/{args.outroot}",
                                settings={'ignore_rows': u'0.0'})
        np.savetxt(f"{args.root}chains/{args.outroot}.covmat",
                   np.array(samples.cov(), dtype='float64'),
                   fmt="%.5e",
                   header=' '.join(names),
                   comments="# ")
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------