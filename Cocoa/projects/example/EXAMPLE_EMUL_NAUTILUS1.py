r"""Sample the example emulator model with the Nautilus nested sampler.

Nested sampling estimates the Bayesian evidence Z, the integral of likelihood
times prior over parameter space (the normalization of the posterior, used to
compare models), and returns weighted posterior samples. It keeps a set of live
points drawn from the prior and keeps replacing the live point of lowest
likelihood by a new point of higher likelihood, so the live set climbs through
nested shells of rising likelihood. Nautilus trains neural networks on the points
evaluated so far to propose the new points.

The model is the YAML in yaml_string, with the likelihoods and parameters of
EXAMPLE_EMUL_MINIMIZE1.py: Planck 2018 high-l TT, TE, EE (plik-lite), low-l TT
and low-l EE (SRoll2), DES-Y5 supernovae, DESI DR2 BAO and ACT DR6 lensing, in
LCDM with the neutrino mass fixed at 0.06 eV. Emulators (neural networks and
Gaussian processes trained on Boltzmann-code outputs) replace CAMB. Seven
parameters are sampled: logA, ns, thetastar = 100 theta_*, omegabh2, omegach2,
tau, and the Planck calibration A_planck (Gaussian prior, mean 1, standard
deviation 0.0025).

How the script maps Cobaya onto Nautilus. The Nautilus prior is uniform on the
central 0.999999 interval of each Cobaya prior: the full range of a flat prior,
and about 4.9 standard deviations around the mean of the Gaussian A_planck
prior. The function Nautilus receives as the likelihood is the full Cobaya log
posterior, log prior + log likelihood. The posterior samples are correct, since
the uniform Nautilus prior only adds a constant. The evidence is not: the prior
density enters twice, and the log Z written in the header equals the
log-evidence minus ln V, where V is the volume of that box (the product of the
seven interval widths, V = 3.1e-7 here, so log Z comes out 15.0 too high).

Inputs (command line; paths are relative to the Cocoa folder):
  --root       folder, ending in "/", whose chains/ subfolder receives the
               output; chains/ must exist before the run, because Nautilus
               writes its checkpoint there during the run. The default,
               ./projects/lsst_y1/, belongs to another project: pass
               ./projects/example/.
  --outroot    basename of the output files (default example_nautilus1).
  --nlive      number of live points (default 1000).
  --nnetworks  number of neural networks Nautilus trains (default 4).
  --flive      the run stops adding shells once the live points hold less than
               this fraction of the evidence (default 0.01).
  --neff       then it keeps sampling until the effective sample size of the
               weighted posterior reaches this value (default 10000).
  --maxfeval   cap on the likelihood evaluations, counted across resumed runs
               (default 100000).

Output files in <root>chains/:
  <outroot>_checkpoint.hdf5  Nautilus's checkpoint; a new run with the same
                      --root and --outroot resumes from it.
  <outroot>.1.txt     the chain in GetDist's text format: weight, Nautilus's log
                      likelihood (here the log posterior), the sampled
                      parameters, and chi2 = -2 log posterior. GetDist reads the
                      second column as -log posterior, so its likelihood
                      statistics (such as the best-fit sample) come out with the
                      wrong sign; the marginalized constraints use only the
                      weights. The header records nlive, maxfeval and log Z.
  <outroot>.ranges    the 0.999999 prior interval of each sampled parameter.
  <outroot>.paramnames  name and LaTeX label of each column after the second.
  <outroot>.covmat    the covariance GetDist computes from the chain, chi2
                      included, with the names in the header.

Run from the Cocoa folder after "source start_cocoa.sh"; python -m
mpi4py.futures starts the worker processes that evaluate the likelihood:

    mpirun -n 12 --oversubscribe python -m mpi4py.futures \
        ./projects/example/EXAMPLE_EMUL_NAUTILUS1.py --root ./projects/example/ \
        --outroot EXAMPLE_EMUL_NAUTILUS1 --maxfeval 450000 --nlive 2048 \
        --neff 15000 --flive 0.01 --nnetworks 5
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
import argparse, random
import numpy as np
from cobaya.yaml import yaml_load
from cobaya.model import get_model
from nautilus import Prior, Sampler
from getdist import loadMCSamples
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Command-line options; the module docstring gives their meaning. Every option
# uses nargs='?' with const=1: written without a value, an option takes the value
# 1 (for --root that is the integer 1, not a path), so always give a value.
parser = argparse.ArgumentParser(prog='LSST_Y1_PROJECT_NAUTILUS1')
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
                    default="example_nautilus1")
parser.add_argument("--nlive",
                    dest="nlive",
                    help="Number of live points ",
                    type=int,
                    nargs='?',
                    const=1,
                    default=1000)
parser.add_argument("--maxfeval",
                    dest="maxfeval",
                    help="Minimizer: maximum number of likelihood evaluations",
                    type=int,
                    nargs='?',
                    const=1,
                    default=100000)
parser.add_argument("--neff",
                    dest="neff",
                    help="Minimum effective sample size. ",
                    type=int,
                    nargs='?',
                    const=1,
                    default=10000)
parser.add_argument("--flive",
                    dest="flive",
                    help="Maximum fraction of the evidence contained in the live set before building the initial shells terminates",
                    type=float,
                    nargs='?',
                    const=1,
                    default=0.01)
parser.add_argument("--nnetworks",
                    dest="nnetworks",
                    help="Number of Neural Networks",
                    type=int,
                    nargs='?',
                    const=1,
                    default=4)
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
# Module-level code runs on every MPI process (mpi4py.futures workers load this
# file too), so each process builds its own Cobaya model and evaluates the
# likelihood calls it receives with it.
model = get_model(yaml_load(yaml_string))
def chi2(p):
    """Return chi2 = -2 log posterior at one point, or 1e20 for a rejected point.

    Every script of this project calls -2 (log prior + log likelihood) "chi2".
    Cobaya's log prior keeps the normalization of each prior: -ln(max - min) for a
    flat prior, and the Gaussian density for A_planck. chi2 therefore differs from
    the likelihood chi2 = -2 log likelihood by a constant plus the A_planck prior
    penalty; the constant cancels in every chi2 difference.

    1e20 is the rejection value: a point outside the prior, or with a non-finite
    likelihood, gets it, and likelihood() maps any chi2 above 1e19 to -inf.

    Arguments:
        p = sampled parameter values in the model's sampled order: a list or an
            array [n_sampled], or a dict whose values are in that order (Nautilus
            passes a dict here).
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

def likelihood(params):
  """Return log prior + log likelihood = -chi2/2, the value Nautilus samples.

  Nautilus uses this value as its log likelihood. It contains Cobaya's log prior,
  which the uniform Nautilus prior does not replace, so the prior density counts
  twice in the evidence (see the module docstring); the posterior shape is
  unaffected. A chi2 above 1e19 (the rejection value of chi2) becomes -inf.

  Arguments:
      params = dict {sampled-parameter name: value} proposed by Nautilus (it
          passes a dict because its prior is a nautilus.Prior), in sampled order.
  Returns:
      the log posterior as a float, or -inf.
  """
  res = chi2(params)
  if (res > 1.e19 or np.isinf(res) or  np.isnan(res)):
    return -np.inf
  else:
    return -0.5*res
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
from mpi4py.futures import MPIPoolExecutor

if __name__ == '__main__':
    print(f"nlive={args.nlive}, output={args.root}chains/{args.outroot}")
    # Nautilus prior: for each sampled parameter, in sampled order, a uniform
    # distribution on the central 0.999999 interval of its Cobaya prior.
    NautilusPrior = Prior()                                       # Nautilus Call 
    dim    = model.prior.d()                                      # Cobaya call
    bounds = model.prior.bounds(confidence=0.999999)              # Cobaya call
    names  = list(model.parameterization.sampled_params().keys()) # Cobaya Call
    for b, name in zip(bounds, names):
      NautilusPrior.add_parameter(name, dist=(b[0], b[1]))

    # pool=MPIPoolExecutor() sends the likelihood calls to the worker processes.
    # resume=True continues from the checkpoint file when it exists.
    sampler = Sampler(NautilusPrior, 
                      likelihood,  
                      filepath=f"{args.root}chains/{args.outroot}_checkpoint.hdf5", 
                      n_dim=dim,
                      pool=MPIPoolExecutor(),
                      n_live=args.nlive,
                      n_networks=args.nnetworks,
                      resume=True)
    # discard_exploration=True drops the points of the first, exploration phase
    # from the posterior and the evidence, as Nautilus requires for unbiased
    # estimates.
    sampler.run(f_live=args.flive,
                n_eff=args.neff,
                n_like_max=args.maxfeval,
                verbose=True,
                discard_exploration=True)
    # points [n, n_sampled]; log_w [n], the log weights, normalized so that the
    # weights sum to one; log_l [n], the value likelihood() returned at each
    # point, that is, the log posterior.
    points, log_w, log_l = sampler.posterior()
    
    # write <outroot>.1.txt, the chain in GetDist's text format ----------------
    # One row per point: weight, log posterior, sampled parameters, chi2. GetDist
    # reads the second column as -log posterior; this file holds +log posterior
    # there, so GetDist's likelihood statistics have the wrong sign (its
    # marginalized constraints use only the weights). The header records log Z,
    # which is offset by -ln V (module docstring).
    os.makedirs(os.path.dirname(f"{args.root}chains/"),exist_ok=True)
    np.savetxt(f"{args.root}chains/{args.outroot}.1.txt",
               np.column_stack((np.exp(log_w), log_l, points, -2*log_l)),
               fmt="%.5e",
               header=f"nlive={args.nlive}, maxfeval={args.maxfeval}, log-Z ={sampler.log_z}\n"+' '.join(names),
               comments="# ")
    
    # write <outroot>.ranges ---------------------------------------------------
    # One line per sampled parameter: name, lower and upper end of its 0.999999
    # prior interval. GetDist reads this file to know where each prior stops, so
    # its smoothed densities do not spread past a prior edge. rows is a list of
    # (name, low, high) tuples; writelines receives a generator expression that
    # formats one text line per tuple.
    rows = [(str(n),float(l),float(h)) for n,l,h in zip(names,bounds[:,0],bounds[:,1])]
    with open(f"{args.root}chains/{args.outroot}.ranges", "w") as f: 
      f.writelines(f"{n} {l:.5e} {h:.5e}\n" for n, l, h in rows)

    # write <outroot>.paramnames -----------------------------------------------
    # One line per column after the second: the name and its LaTeX label from the
    # YAML. In GetDist's convention the trailing * of chi2* marks a derived
    # parameter. names gains chi2* here, so the covariance header includes it.
    param_info = model.info()['params']
    latex  = [param_info[x]['latex'] for x in names]
    names.append("chi2*")
    latex.append("\\chi^2")
    np.savetxt(f"{args.root}chains/{args.outroot}.paramnames", 
               np.column_stack((names,latex)),
               fmt="%s")

    # write <outroot>.covmat ---------------------------------------------------
    # GetDist reloads the weighted chain just written (ignore_rows 0: nested
    # sampling has no burn-in) and computes the covariance of the parameter
    # columns, chi2 included.
    samples = loadMCSamples(f"{args.root}chains/{args.outroot}",
                            settings={'ignore_rows': u'0.0'})
    np.savetxt(f"{args.root}chains/{args.outroot}.covmat",
               np.array(samples.cov(), dtype='float64'),
               fmt="%.7e",
               header=' '.join(names),
               comments="# ")
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------