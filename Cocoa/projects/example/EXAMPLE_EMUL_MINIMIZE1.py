r"""Find the best fit of the example emulator model by annealed emcee sampling.

A minimizer searches parameter space for the point where a function is smallest.
Here the function is chi2 = -2 log posterior = -2 (log prior + log likelihood), so
the point found is the best fit (maximum posterior) of the model in yaml_string:

- likelihoods: Planck 2018 high-l TT, TE, EE (plik-lite), low-l TT and low-l EE
  (SRoll2); DES-Y5 supernovae; DESI DR2 BAO; ACT DR6 CMB lensing;
- cosmology: LCDM, with the neutrino mass fixed at 0.06 eV (w = -1, wa = 0);
- sampled parameters (7): logA = ln(10^10 A_s), ns, thetastar = 100 theta_*,
  omegabh2, omegach2, tau, and the Planck calibration A_planck, which the Planck
  likelihoods add with a Gaussian prior (mean 1, standard deviation 0.0025).

Emulators replace the Boltzmann code (CAMB). An emulator is a neural network or a
Gaussian process trained on Boltzmann-code outputs and evaluated much faster than
the code: emulcmb predicts the CMB TT, TE, EE and lensing-potential spectra,
emultheta maps (omegabh2, omegach2, thetastar) to H0 and Omega_m, emulrdrag
predicts the drag-epoch sound horizon r_drag, and emulbaosn predicts H(z), from
which the BAO and supernova distances follow. device "cuda" falls back to the
Apple GPU (mps) or to the CPU on a machine without a CUDA GPU.

Method. emcee is an ensemble sampler: a set of points called walkers moves through
parameter space, and each proposed step of a walker is built from the positions
of other walkers. Annealing samples the posterior raised to the power 1/T: at
temperature T = 1 the walkers explore the posterior, and at small T they crowd
around the minimum. The temperature ladder 1.0, 0.25, 0.1, 0.005, 0.001 therefore
goes from exploring to polishing (min_chi2 gives the details). The module
docstring of external_modules/code/cosmolike_core/cocoa_hybrid_sampling.py
explains the same annealing scheme (ladder, cov*T/3 starting clouds, DE and
snooker moves) and the reasons for its choices; that module is the shared version
run by the hybrid (EMUL2) Cosmolike projects, and this script is a standalone one.

Inputs (command line; paths are relative to the Cocoa folder):
  --root     folder, ending in "/", that receives chains/ and holds the --cov and
             --minfile files. The default, ./projects/lsst_y1/, belongs to
             another project: pass ./projects/example/.
  --outroot  basename of the output file (default example_min1).
  --nstw     emcee steps per walker at each temperature (default 200).
  --cov, --minfile
             restart files, both relative to --root: a covariance matrix and an
             earlier output of this script. With both, the run starts at the
             saved minimum and runs only the cold stages 0.1, 0.005, 0.001.
             With either one missing, both are ignored: the run starts from a
             random draw of the "ref" distributions, with the prior covariance.
  --emoves   emcee move set: 0 (default) = differential evolution (DE) 80% and
             DE snooker 20%; 1 = DE only; 2 = walk; 3 = KDE (min_chi2 explains).

Output: <root>chains/<outroot>.txt with one row: the best point (the sampled
parameters, in sampled order), the chi2 of each likelihood, the prior term
-2 log prior, and the total chi2. The first header line records nstw, the second
names the columns.

Cost: about 5 * nstw * nwalkers likelihood evaluations (3 * nstw * nwalkers for a
restart). nwalkers = max(3 * n_sampled, number of MPI processes): 21 for the seven
parameters here with at most 21 processes. No random seed is set, so two runs
start from different walkers.

Run from the Cocoa folder after "source start_cocoa.sh". Rank 0 coordinates and
the other MPI processes evaluate the walkers:

    mpirun -n 21 --oversubscribe python ./projects/example/EXAMPLE_EMUL_MINIMIZE1.py \
        --root ./projects/example/ --outroot EXAMPLE_EMUL_MIN1 --nstw 200
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
from schwimmbad import MPIPool
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Command-line options; the module docstring gives their meaning. Every option
# uses nargs='?' with const=1: written without a value, an option takes the value
# 1 (for --root that is the integer 1, not a path), so always give a value.
parser = argparse.ArgumentParser(prog='EXAMPLE_MINIMIZE1')
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
parser.add_argument("--cov",
                    dest="cov",
                    help="Chain Covariance Matrix",
                    nargs='?',
                    const=1,
                    default=None)
parser.add_argument("--minfile",
                    dest="minfile",
                    help="Previous Minimization Result",
                    nargs='?',
                    const=1,
                    default=None)
parser.add_argument("--emoves",
                    dest="emoves",
                    help="Type of Emcee moves",
                    type=int,
                    nargs='?',
                    const=1,
                    default=0)
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
             nstw=200,
             nwalkers=5,
             pool=None,
             emoves=0,
             lowT=True):
    """Minimize chi2 by annealed emcee sampling; return the best point and its chi2.

    At temperature T the walkers sample the posterior raised to the power 1/T,
    whose log is -chi2/(2T). T = 1 samples the posterior itself. For a nearly
    Gaussian posterior the walkers then sit about n_dim*T above the minimum in
    chi2: about 0.007 at T = 0.001 with seven parameters. Each stage of the ladder
    draws nwalkers starting points from a Gaussian centered on x0 with covariance
    cov*T/3, runs nstw emcee steps, and passes its best sample to the next stage
    as the new x0. The function returns the best of the stage optima.

    Why cov*T/3: posterior^(1/T) of a Gaussian posterior with covariance C has
    covariance T*C. When cov is close to the posterior covariance, the starting
    cloud thus has one third of the variance the stage samples (0.58 times its
    standard deviation) and starts inside the region the stage refines. With the
    prior covariance the clouds shrink from 0.58 prior standard deviations at
    T = 1 to 0.018 at T = 0.001.

    Arguments:
        x0 = starting point [n_sampled], in the model's sampled order.
        cov = covariance [n_sampled, n_sampled] in the same order; the row and the
            column of a fixed parameter are removed before walkers are drawn.
        fixed = index of a parameter held at x0[fixed], or -1 for none (this
            script passes -1).
        nstw = emcee steps per walker at each temperature.
        nwalkers = number of walkers.
        pool = schwimmbad MPIPool whose worker processes evaluate the walkers, or
            None to evaluate them in this process.
        emoves = emcee move set: 0 = DE 80% and DE snooker 20%, 1 = DE, 2 = walk,
            3 = KDE.
        lowT = True: ladder 1.0, 0.25, 0.1, 0.005, 0.001, for a start from a
            random point. False: ladder 0.1, 0.005, 0.001, for a restart from a
            saved minimum, whose basin is already found.
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
    # Temperature ladder; each stage starts its walkers from N(x0, cov*T/3), see
    # the docstring for both choices.
    if lowT == True:
      nstw        = int(nstw)
      temperature = np.array([1.0, 0.25, 0.1, 0.005, 0.001], dtype='float64')
      stepsz      = temperature/3.0
    else:
      nstw        = int(nstw)
      temperature = np.array([0.1, 0.005, 0.001], dtype='float64')
      stepsz      = temperature/3.0

    # emcee moves, the rule that builds a walker's proposal from other walkers.
    # DE (differential evolution) steps along the difference of two other walkers;
    # DE snooker uses three other walkers and follows long, curved degeneracies
    # better; walk draws a Gaussian step with the covariance of the other half of
    # the walkers; KDE draws from a kernel density estimate of the other half. A
    # list of (move, weight) pairs picks one move at random for each step.
    if emoves == 0:
      moves = [(emcee.moves.DEMove(), 0.8), (emcee.moves.DESnookerMove(), 0.2)]
    elif emoves == 1:
      moves = emcee.moves.DEMove()
    elif emoves == 2:
      moves = emcee.moves.WalkMove()
    else:
      moves = emcee.moves.KDEMove()

    # partial_samples[i] and partial[i] = best point of stage i and its chi2 at
    # T = 1 (the args tuple of this call keeps T = 1).
    partial_samples = []
    partial = []
    for i in range(len(temperature)):
        # Starting walkers of this stage: nwalkers draws from N(x0, cov*T/3).
        # GaussianStep(...)(x0) has shape [1, n_dim]; [0,:] takes its only row.
        x = [] # Initial point
        for j in range(nwalkers):
            x.append(GaussianStep(stepsize=stepsz[i])(x0)[0,:])  
        # emcee calls logprob(position, *args), with args = (fixed value, fixed
        # index, temperature of this stage); pool spreads the evaluations over the
        # MPI workers. Each step updates the walkers in two halves (DE, walk, KDE)
        # or four quarters (DE snooker), one group after the other, so at most
        # nwalkers/2 evaluations run at the same time.
        sampler = emcee.EnsembleSampler(nwalkers, 
                                        ndim, 
                                        logprob, 
                                        args=(args[0], args[1], temperature[i]),
                                        moves=moves,
                                        pool=pool)    
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
        # The progress line shows the best point over the stages run so far.
        j = np.argmin(np.array(partial))
        print(f"Partial ({i+1}/{len(temperature)}): "
              f"params = {partial_samples[j]}, and chi2 = {partial[j]}")
    # The result is the best of the stage optima, not the last stage's point.
    j = np.argmin(np.array(partial))
    result = [partial_samples[j], partial[j]]
    return result
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
def prf(x0, nstw, cov, fixed=-1, nwalkers=5, pool=None, emoves=0, lowT=True):
    """Run min_chi2 from x0, converted to a float64 array; return its result.

    Arguments:
        x0, nstw, cov, fixed, nwalkers, pool, emoves, lowT = as in min_chi2.
    Returns:
        [best point, its chi2], as min_chi2 returns them.
    """
    res =  min_chi2(x0=np.array(x0, dtype='float64'), 
                    fixed=fixed,
                    cov=cov, 
                    nstw=nstw, 
                    nwalkers=nwalkers,
                    pool=pool,
                    emoves=emoves,
                    lowT=lowT)
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
        dim      = model.prior.d()     
        nstw = args.nstw
        emoves = args.emoves #emcee moves
        # dim = model.prior.d() = number of sampled parameters. Walkers: three per
        # sampled parameter (emcee's moves refuse fewer than two per parameter),
        # and at least as many as MPI processes; Get_size() counts every process,
        # rank 0 included. The KDE move (emoves = 3) estimates a density from half
        # of the walkers, so it gets six per parameter.
        if emoves == 3:
          nwalkers = max(6*dim, pool.comm.Get_size()) # KDE MOVES REQUIRE MORE WALKERS
        else:
          nwalkers = max(3*dim, pool.comm.Get_size())
         
        # 1st: starting point and covariance -----------------------------------
        # With both --cov and --minfile (paths relative to --root) the run restarts:
        # - cov = the first n_sampled rows and columns of the --cov file, which must
        #   be the sampled parameters in sampled order (its header is not read);
        # - x0 = the first n_sampled columns of the --minfile row, whose last
        #   column is its chi2. Recomputing that chi2 must agree within 0.02; a
        #   larger difference means the model differs from the run that wrote the
        #   file (another YAML, other emulator files, another parameter order).
        # Otherwise x0 is a random draw from the "ref" distributions (max_tries
        # bounds the number of draws until the posterior is finite), and cov is
        # the prior covariance: diagonal, each entry the variance of one prior
        # ((max - min)^2/12 for a flat prior).
        if args.cov is not None and args.minfile is not None:
          cov = np.loadtxt(args.root+args.cov)[0:dim,0:dim]
          tmp = np.loadtxt(args.root+args.minfile)
          x0 = tmp[0:model.prior.d()]
          if (abs(chi2(x0)-tmp[-1])>0.02):
            raise ValueError("Inconsistency min file with current setup")
          lowT = False
        else:
          cov = model.prior.covmat(ignore_external=False) # cov from prior 
          (x0, results) = model.get_valid_point(max_tries=1000, 
                                                ignore_fixed_ref=False,
                                                logposterior_as_dict=True)
          lowT = True

        # 2nd: annealed minimization -------------------------------------------
        # min_chi2 reimplements the Procoli method (arXiv:2401.14225).
        res = np.array(list(prf(np.array(x0, dtype='float64'), 
                                fixed=-1, 
                                nstw=nstw,
                                nwalkers=nwalkers,
                                pool=pool,
                                cov=cov,
                                emoves=emoves,
                                lowT=lowT)), dtype="object")
        # res = [best point, chi2] as an object array (its two entries have
        # different shapes); xf is the best point as a one-row table [1, n_sampled].
        xf = np.array([res[0]],dtype='float64')
        
        # 3rd: append each likelihood's chi2, -2 log prior, total chi2 ---------
        xf = np.column_stack((xf, 
                              np.array([chi2v2(d) for d in xf], dtype='float64'),
                              res[1]))
        
        # 4th: write <root>chains/<outroot>.txt --------------------------------
        names = list(model.parameterization.sampled_params().keys()) # Cobaya Call
        names = names+list(model.info()['likelihood'].keys())+["prior"]+["chi2"]
        # dirname of a path that ends in "/" is the folder itself: <root>chains.
        os.makedirs(os.path.dirname(f"{args.root}chains/"),exist_ok=True)
        np.savetxt(f"{args.root}chains/{args.outroot}.txt", 
                   xf,
                   fmt="%.7e",
                   header=f"nswt (evals/Temp/walker)={nstw}\n"+' '.join(names),
                   comments="# ")
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------