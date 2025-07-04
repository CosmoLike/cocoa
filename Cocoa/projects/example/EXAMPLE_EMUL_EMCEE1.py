import warnings
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
    message=r".*invalid value encountered in subtract.*",
    module=r"emcee\.moves\.mh"
)
warnings.filterwarnings(
    "ignore",
    category=RuntimeWarning,
    message=r".*overflow encountered in exp.*"
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
info_txt = r"""
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
    latex: \log(10^{10} A_\mathrm{s})
  ns:
    prior:
      min: 0.92
      max: 1.05
    ref:
      dist: norm
      loc: 0.96605
      scale: 0.005
    proposal: 0.005
    latex: n_\mathrm{s}
  thetastar:
    prior:
      min: 1
      max: 1.2
    ref:
      dist: norm
      loc: 1.04109
      scale: 0.0004
    proposal: 0.0002
    latex: 100\theta_\mathrm{*}
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
    latex: \Omega_\mathrm{b} h^2
  omegach2:
    prior:
      min: 0.06
      max: 0.2
    ref:
      dist: norm
      loc: 0.12011
      scale: 0.03
    proposal: 0.03
    latex: \Omega_\mathrm{c} h^2
  tau:
    prior:
      dist: norm
      loc: 0.0544
      scale: 0.0073
    ref:
      dist: norm
      loc: 0.055
      scale: 0.006
    proposal: 0.003
    latex: \tau_\mathrm{reio}
  As:
    derived: 'lambda logA: 1e-10*np.exp(logA)'
    latex: A_\mathrm{s}
  A:
    derived: 'lambda As: 1e9*As'
    latex: 10^9 A_\mathrm{s}
  mnu:
    value: 0.06
  w0pwa:
    value: -1.0
    latex: w_{0,\mathrm{DE}}+w_{a,\mathrm{DE}}
    drop: true
  w:
    value: -1.0
    latex: w_{0,\mathrm{DE}}
  wa:
    value: 'lambda w0pwa, w: w0pwa - w'
    derived: false
    latex: w_{a,\mathrm{DE}}
  H0:
    derived: true
    latex: H_0
  omegamh2:
    derived: true
    value: 'lambda omegach2, omegabh2, mnu: omegach2+omegabh2+(mnu*(3.046/3)**0.75)/94.0708'
    latex: \Omega_\mathrm{m} h^2
  omegam:
    derived: true
    latex: \Omega_\mathrm{m}
  rdrag:
    derived: true
    latex: r_\mathrm{drag}

theory:
  emultheta:
    path: ./cobaya/cobaya/theories/
    stop_at_error: True
    provides: ['H0', 'omegam']
    extra_args:
      file: ['external_modules/data/emultrf/CMB_TRF/thetaH0GP.joblib']
      extra: ['external_modules/data/emultrf/CMB_TRF/extrainfotheta.npy']
      ord: [['omegabh2','omegach2','thetastar']]
      extrapar: [{'MLA' : "GP"}]
  emulrdrag:
    path: ./cobaya/cobaya/theories/
    stop_at_error: True
    provides: ['rdrag']
    extra_args:
      file: ['external_modules/data/emultrf/BAO_SN_RES/rdragGP.joblib'] 
      extra: ['external_modules/data/emultrf/BAO_SN_RES/extrainfordrag.npy'] 
      ord: [['omegabh2','omegach2']]
  emulcmb:
    path: ./cobaya/cobaya/theories/
    stop_at_error: True
    extra_args:
      # This version of the emul was not trained with CosmoRec
      eval: [True, True, True, False] #TT,TE,EE,PHIPHI
      device: "cuda"
      file: ['external_modules/data/emultrf/CMB_TRF/chiTTAstautrf1dot2milnewlhcevansqrtrescalec16.pt',
             'external_modules/data/emultrf/CMB_TRF/chiTEAstautrf1dot2milnewlhcevansqrtrescalec16.pt',
             'external_modules/data/emultrf/CMB_TRF/chiEEAstautrf1dot2milnewlhcevansqrtrescalec16.pt',
             None] 
      extra: ['external_modules/data/emultrf/CMB_TRF/extrainfo_lhs_tt_96.npy',
              'external_modules/data/emultrf/CMB_TRF/extrainfo_lhs_te_96.npy',
              'external_modules/data/emultrf/CMB_TRF/extrainfo_lhs_ee_96.npy',
              None]
      ord: [['omegabh2','omegach2','H0','tau','ns','logA','mnu','w','wa'],
            ['omegabh2','omegach2','H0','tau','ns','logA','mnu','w','wa'],
            ['omegabh2','omegach2','H0','tau','ns','logA','mnu','w','wa'],
            None]
      extrapar: [{'ellmax' : 5000, 'MLA': 'TRF', 'NCTRF': 16, 'INTDIM': 4, 'INTTRF': 5120},
                 {'ellmax' : 5000, 'MLA': 'TRF', 'NCTRF': 16, 'INTDIM': 4, 'INTTRF': 5120},
                 {'ellmax' : 5000, 'MLA': 'TRF', 'NCTRF': 16, 'INTDIM': 4, 'INTTRF': 5120},
                 None]
  emulbaosn:
    path: ./cobaya/cobaya/theories/
    stop_at_error: True
    extra_args:
      eval: [True, True] # dL(z), H(z)
      device: "cuda"
      method: ["INT", "NN"] # with INT, we will integrate H(z)
      file:  [None,'external_modules/data/emultrf/BAO_SN_RES/chisquareH.pt']
      extra: [None,'external_modules/data/emultrf/BAO_SN_RES/extrainfoH.npy']
      tmat:  [None,'external_modules/data/emultrf/BAO_SN_RES/PCAH.npy']
      zlin:  [None,'external_modules/data/emultrf/BAO_SN_RES/zlinlcdm.npy']       
      ord: [None, ['omegam','H0']]
      extrapar: [{'ZMIN' : 0.0001, 'ZMAX' : 3, 'NZ' : 1200},
                 {'offset' : 0.0000, 'INTDIM' : 1, 'NLAYER' : 1}]
"""
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

info  = yaml_load(info_txt)
model = get_model(info)

def chi2(p):
    point = dict(zip(model.parameterization.sampled_params(),
                 model.prior.sample(ignore_external=True)[0]))
    point.update({'logA': p[0], 
                  'ns':  p[1],
                  'thetastar': p[2], 
                  'omegabh2': p[3], 
                  'omegach2': p[4], 
                  'tau': p[5],
                  'A_planck': p[6]})
    logposterior = model.logposterior(point, as_dict=True)
    if 'loglikes' in logposterior.keys(): 
      res1 = logposterior["loglikes"].get("planck_2018_highl_plik.TTTEEE_lite",-1e200)
      res2 = logposterior["loglikes"].get("planck_2018_lowl.TT",-1e200)
      res3 = logposterior["loglikes"].get("planck_2018_lowl.EE",-1e200)
      res4 = logposterior["loglikes"].get("sn.desy5",-1e200)
      res5 = logposterior["loglikes"].get("bao.desi_dr2.desi_bao_all",-1e200)
      res  = -2.0*(res1 + res2 + res3 + res4 + res5)
    else:
      res = -2.0*(1e20)
    return res

x0 = np.array([
    3.045845885,         # logA
    9.652308970e-01,     # ns
    1.0410562,           # thetastar
    2.246801442e-02,     # omegabh2
    1.198257361e-01 ,    # omegach2
    5.433339482e-02,     # tau
    1.00138
], dtype='float64')

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Code below does not require changes ------------------------------------------
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

parser.add_argument("--cov",
                    dest="cov",
                    help="Chain Covariance Matrix",
                    nargs='?',
                    const=1) # zero or one

# need to use parse_known_args because of mpifuture 
args, unknown = parser.parse_known_args()

maxfeval = args.maxfeval
oroot    = args.root + "chains/" + args.outroot
cov_file = args.root + args.cov
cov      = np.loadtxt(cov_file)[0:len(x0),0:len(x0)]

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

def chain(maxfeval=3000, 
          pool=None):

    def mychi2(params, *args):
        return chi2(p=params)

    def log_prior(params):
        return 1.0
    
    def logprob(params, *args):
        lp = log_prior(params)
        if not np.isfinite(lp):
            return -np.inf
        else:
            return -0.5*mychi2(params, *args) + lp

    class GaussianStep:
       def __init__(self, stepsize=0.2):
           self.cov = stepsize*cov
       def __call__(self, x):
           return np.random.multivariate_normal(x, self.cov, size=1) 

    ndim        = int(x0.shape[0])
    nwalkers    = 2*ndim

    sampler = emcee.EnsembleSampler(int(nwalkers), 
                                    ndim, 
                                    logprob, 
                                    moves=[(emcee.moves.GaussianMove(cov=cov), 1.)],
                                    pool=pool)
    x = [] # Initial point
    for j in range(nwalkers):
        x.append(GaussianStep(stepsize=0.25)(x0)[0,:])
    x = np.array(x,dtype='float64')
    sampler.run_mcmc(x, maxfeval)

    xf   = sampler.get_chain(flat=True, discard=0)
    lnpf = sampler.get_log_prob(flat=True, discard=0)
    return np.concatenate([xf, lnpf[:, None]],axis=1)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

if __name__ == '__main__':
    with MPIPool() as pool:
        if not pool.is_master():
            pool.wait()
            sys.exit(0)

        print(f"maxfeval={maxfeval}")

        res = chain(maxfeval=maxfeval, pool=pool)

        rnd = random.randint(0,9999)
        print("Output file = ", oroot + "_" + str(rnd) + ".txt")
        np.savetxt(oroot + "_" + str(rnd) + ".txt",
                   res,
                   header=f"maxfeval={maxfeval}",
                   comments="# ")

#HOW TO CALL THIS SCRIPT
#mpirun -n 5 --oversubscribe --mca pml ^ucx  \
#  --mca btl vader,tcp,self --bind-to core:overload-allowed \
#  --rank-by slot --map-by core:pe=${OMP_NUM_THREADS}  \
#  python ./projects/example/EXAMPLE_EMUL_EMCEE1.py \
#  --maxfeval 15000 --outroot "example_emul_emcee1" \
#  --cov="EXAMPLE_EMCEE.covmat"