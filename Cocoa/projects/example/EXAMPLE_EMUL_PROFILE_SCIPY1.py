import warnings
import os
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
args, unknown = parser.parse_known_args()
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
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
    p = [float(v) for v in p.values()] if isinstance(p, dict) else p
    if np.any(np.isinf(p)) or  np.any(np.isnan(p)):
      return 1e20
    point = dict(zip(model.parameterization.sampled_params(), p))
    res1 = model.logprior(point,make_finite=False)
    if np.isinf(res1) or  np.any(np.isnan(res1)):
      return 1e20
    res2 = model.loglike(point,
                         make_finite=False,
                         cached=False,
                         return_derived=False)
    if np.isinf(res2) or  np.any(np.isnan(res2)):
      return 1e20
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
if __name__ == '__main__':
    maxevals = args.maxfeval

    # 1st: load the cov. matrix (from running EXAMPLE_EMUL_MCMC1.yaml) ---------
    cov = np.loadtxt(args.root+args.cov[0])[0:model.prior.d(),0:model.prior.d()]
    factor = args.factor
    sigma = np.sqrt(np.diag(cov))

    # 2nd: read minimum --------------------------------------------------------
    x0 = np.loadtxt(args.minfile[0])
    chi20 = x0[-1]
    x0 = x0[0:model.prior.d()]
    
    # 3rd: set the parameter profile range -------------------------------------
    start = np.zeros(model.prior.d(), dtype='float64')
    stop  = np.zeros(model.prior.d(), dtype='float64')
    start = x0 - factor*sigma
    stop  = x0 + factor*sigma
    
    # We need to respect the YAML priors
    bounds0 = model.prior.bounds(confidence=0.999999)
    for i in range(model.prior.d()):
        if (start[i] < bounds0[i][0]):
          start[i] = bounds0[i][0]
        if (stop[i] > bounds0[i][1]):
          stop[i] = bounds0[i][1]

    half_range = (stop[args.profile] - start[args.profile]) / 2.0
   
    numpts = args.numpts-1 if args.numpts%2 == 1 else args.numpts 
  
    param  = np.linspace(start = x0[args.profile] - half_range,
                         stop  = x0[args.profile] + half_range,
                         num = numpts)
    numpts=numpts+1
    param = np.insert(param, numpts//2, x0[args.profile])

    # 4th Print to the terminal ------------------------------------------------
    names = list(model.parameterization.sampled_params().keys()) # Cobaya Call
    print(f"maxfeval={args.maxfeval}, param={names[args.profile]}")
    print(f"profile param values = {param}")
    
    # 5th: Set the vectors that will hold the final result ---------------------
    xf = np.tile(x0, (numpts, 1))
    xf[:,args.profile] = param
    
    chi2res = np.zeros(numpts)
    chi2res[numpts//2] = chi20
    
    # 6th: run from midpoint to right ------------------------------------------
    tmp = np.array(xf[numpts//2,:], dtype='float64')
    for i in range(numpts//2+1,numpts):
        tmp[args.profile] = param[i]
        def chi2_local(x):
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
    
    # 7th: run from midpoint to left -------------------------------------------
    tmp = np.array(xf[numpts//2,:], dtype='float64')
    for i in range(numpts//2-1, -1, -1):
        tmp[args.profile] = param[i]
        def chi2_local(x):
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

    # 8th Append derived parameters --------------------------------------------
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
    
    # 9th Save output file -----------------------------------------------------
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