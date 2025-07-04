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
x = np.array([
    3.045845885,         # logA
    9.652308970e-01,     # ns
    1.0410562,           # thetastar
    2.246801442e-02,     # omegabh2
    1.198257361e-01 ,    # omegach2
    5.433339482e-02,     # tau
    1.00138
], dtype='float64')

bounds0 = np.array([
  [+1.610, +3.91],
  [+0.920, +1.05],
  [+1.02,  +1.06],
  [+0.010, +0.04],
  [+0.060, +0.20],
  [+0.040, +0.09],
  [+0.500, +2.00]
], dtype='float64') 

name  = [ 
    "logAs",     # As
    "ns",        # ns
    "thetastar", # thetastar
    "omegabh2",  # omegabh2
    "omegach2",  # omegach2
    "tau",       # tau
    "A_planck",  # A_planck
]

start  = np.zeros(len(x), dtype='float64')
stop   = np.zeros(len(x), dtype='float64')

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

def chi2v2(p):
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
      res1 = logposterior["loglikes"].get("planck_2018_highl_plik.TTTEEE_lite",-1e20)
      res2 = logposterior["loglikes"].get("planck_2018_lowl.TT",-1e20)
      res3 = logposterior["loglikes"].get("planck_2018_lowl.EE",-1e20)
      res4 = logposterior["loglikes"].get("sn.desy5",-1e20)
      res5 = logposterior["loglikes"].get("bao.desi_dr2.desi_bao_all",-1e20)
      return [-2.0*res1, 
              -2.0*res2, 
              -2.0*res3,
              -2.0*res4,
              -2.0*res5]
    else:
      res [-2.0*(1e20), -2.0*(1e20), -2.0*(1e20), -2.0*(1e20), -2.0*(1e20)]

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

parser.add_argument("--profile",
                    dest="profile",
                    help="Which Parameter to Profile",
                    type=int,
                    nargs='?',
                    const=1,
                    default=1)

parser.add_argument("--numpts",
                    dest="numpts",
                    help="Number of Points to Compute Minimum",
                    type=int,
                    nargs='?',
                    const=1,
                    default=20)

parser.add_argument("--factor",
                    dest="factor",
                    help="Factor that set the bounds (multiple of cov matrix)",
                    type=int,
                    nargs='?',
                    const=1,
                    default=3)

parser.add_argument("--cov",
                    dest="cov",
                    help="Chain Covariance Matrix",
                    nargs='?',
                    const=1) # zero or one

parser.add_argument("--nwalkers",
                    dest="nwalkers",
                    help="Number of emcee walkers",
                    nargs='?',
                    const=1)

# need to use parse_known_args because of mpifuture 
args, unknown = parser.parse_known_args()

maxfeval = args.maxfeval
oroot    = args.root + "chains/" + args.outroot
index    = args.profile
numpts   = args.numpts
nwalkers = args.nwalkers

cov_file = args.root + args.cov
cov      = np.loadtxt(cov_file)[0:len(x),0:len(x)]
sigma    = np.sqrt(np.diag(cov))
start    = x - args.factor*sigma
stop     = x + args.factor*sigma
for i in range(len(x)):
    if (start[i] < bounds0[i][0]):
      start[i] = bounds0[i][0]
    if (stop[i] > bounds0[i][1]):
      stop[i] = bounds0[i][1]

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

def min_chi2(x0, 
             bounds, 
             fixed=-1, 
             maxfeval=3000, 
             cov=cov,
             nwalkers=5):

    def mychi2(params, *args):
        z, fixed, T = args
        params = np.array(params, dtype='float64')
        if fixed > -1:
            params = np.insert(params, fixed, z)
        return chi2(p=params)/T

    if fixed > -1:
        z      = x0[fixed]
        x0     = np.delete(x0, (fixed))
        bounds = np.delete(bounds, (fixed), axis=0)
        args = (z, fixed, 1.0)
    else:
        args = (0.0, -2.0, 1.0)

    def log_prior(params):
        return 1.0
    
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
    
    ndim        = int(x0.shape[0])
    nwalkers    = int(2*x0.shape[0])
    nsteps      = maxfeval
    temperature = np.array([1.0, 0.25, 0.1, 0.005, 0.001], dtype='float64')
    stepsz      = temperature/4.0

    mychi2(x0, *args) # first call takes a lot longer (when running on cuda)
    start_time = time.time()
    mychi2(GaussianStep(stepsize=0.1)(x0)[0,:], *args)
    elapsed_time = time.time() - start_time

    print(f"Emcee: nwalkers = {nwalkers}, "
          f"nTemp = {len(temperature)}, "
          f"feval (per walker) = {maxfeval}, "
          f"feval (per Temp) = {nwalkers*maxfeval}, "
          f"feval = {nwalkers*maxfeval*len(temperature)}")
    print(f"Emcee: Like Eval Time: {elapsed_time:.4f} secs, "
          f"Eval Time = {elapsed_time*nwalkers*maxfeval*len(temperature)/3600.:.4f} hours.")

    partial_samples = []
    partial = []

    for i in range(len(temperature)):
        x = [] # Initial point
        for j in range(nwalkers):
            x.append(GaussianStep(stepsize=stepsz[i])(x0)[0,:])
        x = np.array(x,dtype='float64')

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
        tchi2 = mychi2(samples[j], *args)
        partial.append(tchi2)
        x0 = copy.deepcopy(samples[j])
        sampler.reset()
    # min chi2 from the entire emcee runs
    j = np.argmin(np.array(partial))
    result = [partial_samples[j], partial[j]]
    return result

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

def prf(x0, index, maxfeval, bounds, nwalkers=5):
    t0 = np.array(x0, dtype='float64')
    t1 = np.array(bounds, dtype="float64") # np.array do a deep copy. Deep copy necessary 
                                           # line to avoid weird bug that changes on bounds
                                           # propagate from different iterations (same MPI core)
    t1[:,0] += t0
    t1[:,1] += t0

    for i in range(len(x)):
        if (t1[i][0] < bounds0[i][0]):
          t1[i][0] = bounds0[i][0]
        if (t1[i][1] > bounds0[i][1]):
          t1[i][1] = bounds0[i][1]

    res =  min_chi2(x0=t0, 
                    bounds=t1, 
                    fixed=index, 
                    maxfeval=maxfeval, 
                    nwalkers=nwalkers)
    return res

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
from cobaya.theories.emultheta.emultheta2 import emultheta

etheta = emultheta(extra_args={ 
    'device': "cuda",
    'file': ['external_modules/data/emultrf/CMB_TRF/thetaH0GP.joblib'],
    'extra':['external_modules/data/emultrf/CMB_TRF/extrainfotheta.npy'],
    'ord':  [['omegabh2','omegach2','thetastar']],
    'extrapar': [{'MLA' : "GP"}]})

# want to add rd information
from cobaya.theories.emulrdrag.emulrdrag2 import emulrdrag

# LCDM emulator. On this mass range, DE does not influence rd
erd = emulrdrag(extra_args={ 
    'file': ['external_modules/data/emultrf/BAO_SN_RES/rdragGP.joblib'] ,
    'extra': ['external_modules/data/emultrf/BAO_SN_RES/extrainfordrag.npy'],
    'ord':  [['omegabh2','omegach2']]})

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

from mpi4py.futures import MPIPoolExecutor

if __name__ == '__main__':
    
    print(f"nwalkers={nwalkers}, maxfeval={maxfeval}, param={index}")

    executor = MPIPoolExecutor()
    
    param = np.linspace(start=start[index], stop=stop[index], num=numpts)
    print(f"profile param values = {param}")
    x0 = np.tile(x, (param.size, 1))
    x0[:,index] = param
    bounds = np.c_[- 25.0*args.factor*sigma,+ 25.0*args.factor*sigma]

    res = np.array(list(executor.map(functools.partial(prf, 
                                                       index=index,
                                                       maxfeval=maxfeval, 
                                                       bounds=bounds, 
                                                       nwalkers=nwalkers),x0)),dtype="object")

    x0 = np.array([np.insert(row,index,p) for row, p in zip(res[:,0],param)],dtype='float64')

    # Append H0 and omegam (begins) --------------------
    tmp = [
        etheta.calculate({
            'thetastar': row[2],
            'omegabh2':  row[3],
            'omegach2':  row[4],
            'omegamh2':  row[3] + row[4] + (0.06*(3.046/3)**0.75)/94.0708
        })
        for row in x0
      ]
    h  = np.array([d['H0']     for d in tmp], dtype='float64')
    om = np.array([d['omegam'] for d in tmp], dtype='float64')
    x0 = np.column_stack((x0, h, om))
    # Append H0 and omegam (ends) --------------------

    # Append rd (begins) --------------------
    tmp = [
        erd.calculate({
            'thetastar': row[2],
            'omegabh2':  row[3],
            'omegach2':  row[4],
            'omegamh2':  row[3] + row[4] + (0.06*(3.046/3)**0.75)/94.0708,
            'H0':        row[-2]
        })
        for row in x0
      ]
    rdrag = np.array([d['rdrag'] for d in tmp], dtype='float64')
    x0 = np.column_stack((x0, rdrag))
    # Append rd (ends) --------------------

    # Append individual chi2 (begins) --------------------
    tmp = np.array([chi2v2(d) for d in x0], dtype='float64')
    x0 = np.column_stack((x0,tmp[:,0],tmp[:,1],tmp[:,2],tmp[:,3],tmp[:,4]))
    # Append individual chi2 (ends) --------------------

    rnd = random.randint(0,9999)
    out = oroot + "_" + str(rnd) + "_" + name[index] 
    print("Output file = ", out + ".txt")
    np.savetxt(out+"_ref"+".txt",
               np.concatenate([np.c_[param,res[:,1]],x0],axis=1),
               header=f"nwalkers={nwalkers}, maxfeval={maxfeval}, param={name[index]}",
               comments="# ")
    executor.shutdown()

#HOW TO CALL THIS SCRIPT
#mpirun -n 5 --oversubscribe --mca pml ^ucx  \
#  --mca btl vader,tcp,self --bind-to core:overload-allowed \
#  --rank-by slot --map-by core:pe=${OMP_NUM_THREADS}  \
#  python -m mpi4py.futures ./projects/example/EXAMPLE_EMUL_PROFILE3.py \
#  --nwalkers 5 --profile 1 --maxfeval 15 --numpts 4  \
#  --outroot "example_emul_profile3" --factor 5 \
#  --cov 'EXAMPLE_EMUL_MCMC2.covmat'