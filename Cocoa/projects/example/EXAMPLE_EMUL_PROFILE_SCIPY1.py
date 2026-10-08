r"""Profile one parameter of the example emulator model with scipy's Nelder-Mead.

A profile fixes one parameter at each value of a grid and, at each value,
minimizes chi2 = -2 log posterior = -2 (log prior + log likelihood) over all the
other parameters. The curve of minimized chi2 against the fixed value is the
profile. For one parameter, the values where it rises by Delta chi2 = 1 and 4
above its minimum bound the 68% and 95% confidence intervals of a Gaussian case.
The priors stay in chi2: flat priors only add a constant, but the Gaussian prior
of A_planck enters the minimized quantity.

This is the second profile method of the example project. Where
EXAMPLE_EMUL_PROFILE1.py runs annealed emcee sampling with MPI, this script runs
scipy's Nelder-Mead minimizer at each grid point, in one process. Nelder-Mead
needs no derivatives: it keeps a simplex of n + 1 points in n dimensions and
replaces its worst point by reflecting, expanding or contracting it through the
others, shrinking the simplex toward its best point when no such move improves.
It finds the nearest local minimum, so it suits a smooth posterior with few
parameters and closely spaced grid points, since each grid point starts from the
optimum of its neighbor. scripts/EXAMPLE_PLOT_PROFILE1_COMP.py compares the two
methods.

The model is the YAML in yaml_string, the same as in EXAMPLE_EMUL_PROFILE1.py:
Planck 2018 high-l TT, TE, EE (plik-lite), low-l TT and low-l EE (SRoll2), DES-Y5
supernovae, DESI DR2 BAO and ACT DR6 lensing, in LCDM with the neutrino mass
fixed at 0.06 eV. Emulators (neural networks and Gaussian processes trained on
Boltzmann-code outputs) replace CAMB. Seven sampled parameters, in sampled order:
logA, ns, thetastar = 100 theta_*, omegabh2, omegach2, tau, and the Planck
calibration A_planck (Gaussian prior, mean 1, standard deviation 0.0025).

Steps, numbered as in the main block: (1) covariance; (2) global minimum; (3)
grid of the profiled parameter, centered on the minimum; (4) print the grid; (5)
result tables; (6, 7) Nelder-Mead minimizations, walking outward from the
minimum; (8) derived parameters and chi2 values; (9) output file.

Inputs (command line; paths are relative to the Cocoa folder):
  --root      folder, ending in "/", that receives chains/ and holds --cov
              (default ./projects/example/).
  --outroot   basename of the output file (default test.dat).
  --profile   zero-based index, in sampled order, of the profiled parameter
              (default 1, ns).
  --numpts    number of grid points besides the minimum (default 20); an odd
              value is lowered by one, so the grid is symmetric.
  --factor    half-width of the grid, in standard deviations sqrt(diag(cov)) of
              the profiled parameter; an integer here (default 3).
  --maxfeval  passed to scipy as maxiter, the cap on Nelder-Mead iterations per
              grid point (default 5000). An iteration costs one or two chi2
              evaluations, more when the simplex shrinks, and no separate cap on
              evaluations is set.
  --cov       required: covariance file relative to --root, such as the .covmat
              of the EXAMPLE_EMUL_MCMC1.yaml chain; its first n_sampled rows and
              columns must be the sampled parameters in sampled order.
  --minfile   required: one-row file, relative to the Cocoa folder (not to
              --root), whose first n_sampled columns are the global minimum and
              whose last column is its chi2, such as the output of
              EXAMPLE_EMUL_MINIMIZE1.py. That chi2 is used for the center grid
              point without being recomputed.
Without --cov or --minfile the script stops with a TypeError (None is not
subscriptable).

Output: <root>chains/<outroot>.<name>.txt, with <name> the profiled parameter,
one row per grid point: the grid value, chi2, the sampled parameters of the
minimization, H0 [km/s/Mpc], Omega_m, the chi2 of each likelihood, and -2 log
prior. The header records maxfeval and names the columns.

Run from the Cocoa folder after "source start_cocoa.sh":

    mpirun -n 1 python ./projects/example/EXAMPLE_EMUL_PROFILE_SCIPY1.py \
        --root ./projects/example/ --cov chains/EXAMPLE_EMUL_MCMC1.covmat \
        --outroot EXAMPLE_EMUL_PROFILE1M2 --factor 3 --maxfeval 5000 --numpts 10 \
        --profile 1 --minfile ./projects/example/chains/EXAMPLE_EMUL_MIN1.txt
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
from scipy import optimize
from cobaya.yaml import yaml_load
from cobaya.model import get_model
from getdist import IniFile
import sys, platform, os
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Command-line options; the module docstring gives their meaning. Options with
# nargs='?' and const=1 take the value 1 when written without a value, so always
# give one. --minfile and --cov use nargs=1: each is a list holding one path.
parser = argparse.ArgumentParser(prog='EXAMPLE_EMUL_PROFILE1')

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
                    type=int,
                    nargs='?',
                    const=1,
                    default=3)
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
                    nargs=1)
parser.add_argument("--cov",
                    dest="cov",
                    help="Chain Covariance Matrix",
                    nargs=1)
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
# parameter has a prior; "ref" and "proposal" are read by Cobaya's samplers, not
# by this script); and the theory block, where the emulators replace CAMB. Paths
# such as ./external_modules are relative to the Cocoa folder, the folder to run
# from.
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
model = get_model(yaml_load(yaml_string))
def chi2(p):
    """Return chi2 = -2 log posterior at one point, or 1e20 for a rejected point.

    Every script of this project calls -2 (log prior + log likelihood) "chi2".
    Cobaya's log prior keeps the normalization of each prior: -ln(max - min) for a
    flat prior, and the Gaussian density for A_planck. chi2 therefore differs from
    the likelihood chi2 = -2 log likelihood by a constant plus the A_planck prior
    penalty; the constant cancels in every chi2 difference.

    1e20 is the rejection value, returned for a point outside the prior, a
    non-finite likelihood, and (unlike the emcee scripts, which raise) a
    non-finite parameter value, so the minimizer always receives a number.

    Arguments:
        p = sampled parameter values in the model's sampled order: a list or an
            array [n_sampled], or a dict whose values are in that order.
    Returns:
        chi2 as a float, or 1e20.
    """
    # A dict is replaced by the list of its values (which must be in sampled
    # order); a list or an array is used as given.
    p = [float(v) for v in p.values()] if isinstance(p, dict) else p
    if np.any(np.isinf(p)) or  np.any(np.isnan(p)):
      return 1e20
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
# Standalone copy of the emultheta emulator of the theory block, used after the
# profile for the derived parameters of the grid points: a Gaussian process that
# maps (omegabh2, omegach2, thetastar) to H0 and, with omegamh2, to Omega_m.
from cobaya.theories.emultheta.emultheta2 import emultheta
etheta = emultheta(extra_args={ 
    'device': "cuda",
    'file': ['external_modules/data/emultrf/CMB_TRF/emul_lcdm_thetaH0_GP.joblib'],
    'extra':['external_modules/data/emultrf/CMB_TRF/extra_lcdm_thetaH0.npy'],
    'ord':  [['omegabh2','omegach2','thetastar']],
    'extrapar': [{'MLA' : "GP"}]})
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# The block below runs when the file is executed as a script. There is no MPI
# pool: every chi2 is evaluated in this one process.
if __name__ == '__main__':
    maxevals = args.maxfeval

    # 1st: covariance (required) -----------------------------------------------
    # The first n_sampled rows and columns of <root> + --cov, for example the
    # .covmat of the EXAMPLE_EMUL_MCMC1.yaml chain; they must be the sampled
    # parameters in sampled order. sigma = standard deviation of each one.
    cov = np.loadtxt(args.root+args.cov[0])[0:model.prior.d(),0:model.prior.d()]
    factor = args.factor
    sigma = np.sqrt(np.diag(cov))

    # 2nd: global minimum (required) -------------------------------------------
    # x0 = the first n_sampled columns of the --minfile row; chi20 = its last
    # column, used as the chi2 of the center grid point without recomputing it.
    print(args.minfile[0])
    x0 = np.loadtxt(args.minfile[0])
    chi20 = x0[-1]
    x0 = x0[0:model.prior.d()]
    
    # 3rd: grid of the profiled parameter --------------------------------------
    # Range from x0 - factor*sigma to x0 + factor*sigma for every parameter; the
    # loop below clips it to the prior.
    start = np.zeros(model.prior.d(), dtype='float64')
    stop  = np.zeros(model.prior.d(), dtype='float64')
    start = x0 - factor*sigma
    stop  = x0 + factor*sigma
    
    # Clip to the central 0.999999 interval of each prior: the full range of a
    # flat prior, about 4.9 standard deviations around the mean of a Gaussian
    # prior (A_planck). bounds0 also bounds the Nelder-Mead search below.
    bounds0 = model.prior.bounds(confidence=0.999999)
    for i in range(model.prior.d()):
        if (start[i] < bounds0[i][0]):
          start[i] = bounds0[i][0]
        if (stop[i] > bounds0[i][1]):
          stop[i] = bounds0[i][1]

    # Only the clipped width of the profiled parameter is kept: the grid is
    # centered on the minimum with half that width. After a clip on one side
    # only, the grid therefore reaches past the prior on that side, and the grid
    # points there get chi2 = 1e20.
    half_range = (stop[args.profile] - start[args.profile]) / 2.0
   
    # Grid: an even number of points from x0 - half_range to x0 + half_range, plus
    # the minimum's own value inserted in the middle, so the minimum is the
    # center row. Example: --numpts 10 (or 11, lowered to 10), x0 = 0.966,
    # half_range = 0.012 gives 0.954, 0.95667, ..., 0.96467, then 0.966 at index
    # 5, then 0.96733, ..., 0.978: 11 sorted points.
    numpts = args.numpts-1 if args.numpts%2 == 1 else args.numpts 
  
    param  = np.linspace(start = x0[args.profile] - half_range,
                         stop  = x0[args.profile] + half_range,
                         num = numpts)
    numpts=numpts+1
    param = np.insert(param, numpts//2, x0[args.profile])

    # 4th: print the grid ------------------------------------------------------
    names = list(model.parameterization.sampled_params().keys()) # Cobaya Call
    print(f"maxfeval={args.maxfeval}, param={names[args.profile]}")
    print(f"profile param values = {param}")
    
    # 5th: result tables -------------------------------------------------------
    # xf [numpts, n_sampled]: one row per grid point, x0 with the profiled column
    # set to the grid value; chi2res [numpts]: the minimized chi2. The center row
    # is the minimum itself.
    xf = np.tile(x0, (numpts, 1))
    xf[:,args.profile] = param
    
    chi2res = np.zeros(numpts)
    chi2res[numpts//2] = chi20
    
    # 6th: minimize from the center to the last grid point ---------------------
    # Each grid point starts from the optimum of the previous one (tmp), which is
    # close to its own optimum. scipy builds the starting simplex from x0 = tmp
    # by moving one coordinate per extra point by 5% of its value. Options:
    # maxiter caps the iterations (--maxfeval); bounds keeps the simplex inside
    # the prior intervals bounds0; tol = 0.01 sets both stopping tolerances: the
    # simplex points must agree within 0.01 in every parameter and within 0.01
    # in chi2. 0.01 exceeds the posterior width of most parameters here (about
    # 1e-4 for omegabh2), so the chi2 tolerance sets the precision. The
    # profiled coordinate of res.x can drift, since chi2_local ignores it; the
    # loop resets it to the grid value.
    tmp = np.array(xf[numpts//2,:], dtype='float64')
    for i in range(numpts//2+1,numpts):
        tmp[args.profile] = param[i]
        def chi2_local(x):
          """Return chi2 at x, with the profiled coordinate reset to the grid value.

          Nelder-Mead moves all n_sampled coordinates. chi2_local overwrites the
          profiled one, on the copy of x that scipy passes, with the grid value
          tmp[args.profile] before evaluating, so moves along that coordinate leave
          chi2 unchanged and the search runs over the other parameters. tmp is read
          from the enclosing loop at call time.

          Arguments:
              x = point [n_sampled] proposed by Nelder-Mead.
          Returns:
              chi2 as a float, or 1e20.
          """
          x[args.profile] = tmp[args.profile]
          return chi2(x)
        res = optimize.minimize(chi2_local, 
                                x0=tmp, 
                                options={'maxiter': maxevals},
                                method='Nelder-Mead',
                                bounds=bounds0,
                                tol=0.01)
        tmp = res.x
        tmp[args.profile] = param[i]
        xf[i,:] = tmp[:]
        chi2res[i] = chi2(tmp)
        print(f"Partial ({i+1}/{numpts}): params = {tmp}, and chi2 = {chi2res[i]}")
    
    # 7th: minimize from the center to the first grid point --------------------
    tmp = np.array(xf[numpts//2,:], dtype='float64')
    for i in range(numpts//2-1, -1, -1):
        tmp[args.profile] = param[i]
        def chi2_local(x):
          """Return chi2 at x, with the profiled coordinate reset to the grid value.

          Nelder-Mead moves all n_sampled coordinates. chi2_local overwrites the
          profiled one, on the copy of x that scipy passes, with the grid value
          tmp[args.profile] before evaluating, so moves along that coordinate leave
          chi2 unchanged and the search runs over the other parameters. tmp is read
          from the enclosing loop at call time.

          Arguments:
              x = point [n_sampled] proposed by Nelder-Mead.
          Returns:
              chi2 as a float, or 1e20.
          """
          x[args.profile] = tmp[args.profile]
          return chi2(x)
        res = optimize.minimize(chi2_local, 
                                x0=tmp, 
                                options={'maxiter': maxevals},
                                method='Nelder-Mead',
                                bounds=bounds0,
                                tol=0.01)
        tmp = res.x
        tmp[args.profile] = param[i]
        xf[i,:] = tmp[:]
        chi2res[i] = chi2(tmp)
        print(f"Partial ({i+1}/{numpts}): params = {tmp}, and chi2 = {chi2res[i]}")       

    # 8th: append derived parameters and chi2 values ---------------------------
    # etheta gives H0 [km/s/Mpc] and Omega_m; the list comprehension below builds
    # one result dict per row of xf. Columns 2, 3 and 4 of xf are thetastar,
    # omegabh2 and omegach2 because they are the third to fifth parameters of the
    # YAML params block; reordering that block breaks these indices. omegamh2
    # adds the neutrino density mnu*(3.046/3)^0.75/94.0708 for mnu = 0.06 eV, the
    # value and the CAMB conversion of omegamh2 in the YAML. The last block of
    # columns is the chi2 of each likelihood and -2 log prior.
    tmp = [
        etheta.calculate({
            'thetastar': row[2],
            'omegabh2':  row[3],
            'omegach2':  row[4],
            'omegamh2':  row[3] + row[4] + (0.06*(3.046/3)**0.75)/94.0708
        })
        for row in xf
      ]
    xf = np.column_stack((xf, 
                          np.array([d['H0'] for d in tmp], dtype='float64'), 
                          np.array([d['omegam'] for d in tmp], dtype='float64'),
                          np.array([chi2v2(d) for d in xf], dtype='float64')))
    
    # 9th: write <root>chains/<outroot>.<name>.txt -----------------------------
    # Columns: grid value, chi2, then the rows of xf; np.c_[a, b] stacks the
    # one-dimensional arrays a and b as the columns of a two-column table.
    os.makedirs(os.path.dirname(f"{args.root}chains/"),exist_ok=True)
    hd = [names[args.profile], "chi2"] + names + ["H0", 'omegam']
    hd = hd + list(model.info()['likelihood'].keys()) + ["prior"]    
    np.savetxt(f"{args.root}chains/{args.outroot}.{names[args.profile]}.txt",
               np.concatenate([np.c_[param, chi2res],xf],axis=1),
               fmt="%.6e",
               header=f"maxfeval={args.maxfeval}, param={names[args.profile]}\n"+' '.join(hd),
               comments="# ")
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------