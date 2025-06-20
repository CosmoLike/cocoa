import warnings
from sklearn.exceptions import InconsistentVersionWarning
warnings.filterwarnings("ignore", category=InconsistentVersionWarning)
warnings.filterwarnings(
    "ignore",
    message=".*column is deprecated.*",
    module=r"sacc\.sacc"
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
      file:  ['external_modules/data/emultrf/BAO_SN_RES/chisquaredllog.pt',
              'external_modules/data/emultrf/BAO_SN_RES/chisquareH.pt']
      extra: ['external_modules/data/emultrf/BAO_SN_RES/extrainfodl.npy',
              'external_modules/data/emultrf/BAO_SN_RES/extrainfoH.npy']
      tmat:  ['external_modules/data/emultrf/BAO_SN_RES/PCAdl.npy',
              'external_modules/data/emultrf/BAO_SN_RES/PCAH.npy']
      zlin:  ['external_modules/data/emultrf/BAO_SN_RES/zlinlcdm.npy',
              'external_modules/data/emultrf/BAO_SN_RES/zlinlcdm.npy']       
      ord: [['omegam','H0'],
            ['omegam','H0']]
      extrapar: [{'offset' : 4400.0, 'INTDIM' : 4, 'NLAYER' : 4},
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

name  = [ 
    "logAs",     # As
    "ns",        # ns
    "thetastar", # thetastar
    "omegabh2",  # omegabh2
    "omegach2",  # omegach2
    "tau",       # tau
    "A_planck",  # A_planck
]

bounds = np.zeros((len(x),2), dtype='float64')
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
    res1 = logposterior["loglikes"].get("planck_2018_highl_plik.TTTEEE_lite",-1e200)
    res2 = logposterior["loglikes"].get("planck_2018_lowl.TT",-1e200)
    res3 = logposterior["loglikes"].get("planck_2018_lowl.EE",-1e200)
    res4 = logposterior["loglikes"].get("sn.desy5",-1e200)
    res5 = logposterior["loglikes"].get("bao.desi_dr2.desi_bao_all",-1e200)
    res  = -2.0*(res1 + res2 + res3 + res4 + res5)
    return res


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Code below does not require changes ------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

parser = argparse.ArgumentParser(prog='EXAMPLE_EMUL_PROFILE1')

parser.add_argument("--tol",
                    dest="tolerance",
                    help="Minimizer Tolerance",
                    type=float,
                    nargs='?',
                    const=1,
                    default=0.05)

parser.add_argument("--maxfeval",
                    dest="maxfeval",
                    help="Minimizer: maximum number of likelihood evaluations",
                    type=int,
                    nargs='?',
                    const=1,
                    default=5000)

parser.add_argument("--maxiter",
                    dest="maxiter",
                    help="Minimizer: maximum number of minimizer iterations",
                    type=int,
                    nargs='?',
                    const=1,
                    default=3)

parser.add_argument("--minmethod",
                    dest="minmethod",
                    help="Minimizer: minimizer method",
                    type=int,
                    nargs='?',
                    const=1,
                    default=5)

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

parser.add_argument("--ref",
                    dest="ref",
                    help="Refinement Factor (min method 2)",
                    type=int,
                    nargs='?',
                    const=1,
                    default=1) # zero or one

parser.add_argument("--cov",
                    dest="cov",
                    help="Chain Covariance Matrix",
                    nargs='?',
                    const=1) # zero or one

# need to use parse_known_args because of mpifuture 
args, unknown = parser.parse_known_args()

tol         = args.tolerance
maxfeval    = args.maxfeval
maxiter     = args.maxiter
min_method  = args.minmethod
oroot       = args.root + "chains/" + args.outroot
index       = args.profile
numpts      = args.numpts
ref         = args.ref

cov_file = args.root + args.cov
cov      = np.loadtxt(cov_file)[0:len(x),0:len(x)]
sigma    = np.sqrt(np.diag(cov))
start    = x - args.factor*sigma
stop     = x + args.factor*sigma
bounds   = np.c_[x - 100.0*args.factor*sigma, x + 100.0*args.factor*sigma]

#for i in range(len(x)):
#    start[i]     = x[i] -args.factor*np.sqrt(cov[i,i])
#    stop[i]      = x[i] +args.factor*np.sqrt(cov[i,i])
#    bounds[i][0] = x[i] -100.0*args.factor*np.sqrt(cov[i,i])  
#    bounds[i][1] = x[i] +100.0*args.factor*np.sqrt(cov[i,i])

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

def min_chi2(x0, 
             bounds, 
             min_method, 
             fixed=-1, 
             tol=0.01, 
             maxfeval=3000, 
             maxiter=10, 
             cov=cov):

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

    if min_method == 1:    
    
        x = copy.deepcopy(x0)

        mychi2(x, *args) # first call takes a lot longer (when running on cuda)
        start_time = time.time()
        mychi2(GaussianStep(stepsize=0.1)(x)[0,:], *args)
        elapsed_time = time.time() - start_time
        print(f"MN: Like Eval Time: {elapsed_time:.4f} secs, "
              f"Eval Time = {elapsed_time*maxiter*maxfeval/3600.:.4f} hours.")

        partial_samples = [x]
        partial = [mychi2(x, *args)]
        for i in range(maxiter):
            tmp = iminuit.minimize(fun=mychi2, 
                                   x0=x, 
                                   args=args, 
                                   bounds=bounds, 
                                   method="migrad", 
                                   tol=tol,
                                   options = {'stra'  : 2, 
                                              'maxfun': maxfeval})
            partial_samples.append(tmp.x)
            partial.append(tmp.fun)
            j = np.argmin(np.array(partial))
            x = GaussianStep(stepsize=0.005)(partial_samples[j])[0,:]
            #print(f"MN: i = {i}, chi2 = {partial[j]}, param = {args[0]}")
        j = np.argmin(np.array(partial))
        result = [partial_samples[j], partial[j]] 
    
    elif min_method == 2: # adapted from PROCOLI       
    
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
            
            sampler.run_mcmc(x, nsteps, skip_initial_state_check=False)
            samples = sampler.get_chain(flat=True, thin=1, discard=0)

            j = np.argmin(-1.0*np.array(sampler.get_log_prob(flat=True)))
            partial_samples.append(samples[j])
            tchi2 = mychi2(samples[j], *args)
            partial.append(tchi2)
            x0 = copy.deepcopy(samples[j])
            sampler.reset()
            #print(f"emcee: i = {i}, chi2 = {tchi2}, param = {args[0]}")
        
        # min chi2 from the entire emcee runs
        j = np.argmin(np.array(partial))
        result = [partial_samples[j], partial[j]]
    
    else:
        raise RuntimeError("Unknown Mimimizer Type") 
    return result

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

def prf(x0, index, min_method, maxiter, maxfeval):
    res =  min_chi2(x0=np.array(x0, dtype='float64'), 
                    bounds=np.array(bounds, dtype='float64'), 
                    min_method=min_method, 
                    fixed=index, 
                    tol=tol, 
                    maxfeval=maxfeval, 
                    maxiter=maxiter)
    return res

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

from mpi4py.futures import MPIPoolExecutor

if __name__ == '__main__':
    
    print(f"min_method={min_method}, maxiter={maxiter}, maxfeval={maxfeval}, tol={tol}, param={index}")

    executor = MPIPoolExecutor()
    
    param = np.linspace(start=start[index], stop=stop[index], num=numpts)
    print(f"profile param values = {param}")

    x0 = []
    for i in range(numpts):
        tmp = copy.deepcopy(x)
        tmp[index] = param[i]
        x0.append(tmp)
    x0 = np.array(x0, dtype='float64')

    res = np.asarray(list(executor.map(functools.partial(prf, index=index, 
        min_method=min_method, maxiter=maxiter, maxfeval=maxfeval), x0)), dtype="object")

    rnd = random.randint(0,1000)
    out = oroot + "_" + str(rnd) + "_method_" + str(min_method) + "_" + name[index] 
    print("Output file = ", out + ".txt")
    np.savetxt(out + ".txt", np.c_[param, res[:,1]])

    if ref > 0 and min_method != 1:
        x0 = []
        for i in range(numpts): 
            x0.append(np.insert(res[i,0], index, param[i]))
        x0 = np.array(x0, dtype=object)
        
        res = np.array(list(executor.map(functools.partial(prf, index=index, 
            min_method=1, maxiter=maxiter, maxfeval=max(10*maxfeval,250)), x0)), dtype="object")
        
        print("Output file = ", out + "_ref" + ".txt")
        np.savetxt(out + "_ref" + ".txt", np.c_[param, res[:,1]])
    executor.shutdown()

#HOW TO CALL THIS SCRIPT
#  mpirun -n 5 --oversubscribe --mca pml ^ucx  \
#    --mca btl vader,tcp,self --bind-to core:overload-allowed \
#    --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS}  \
#    python -m mpi4py.futures ./projects/example/EXAMPLE_EMUL_PROFILE3.py \
#    --tol 0.05 --profile 1 --maxiter 2 --maxfeval 50 --numpts 4 \
#    --outroot "example_emul_profile3" --minmethod 2 --factor 5 --ref 1 \
#    --cov 'EXAMPLE_EMUL_MCMC2.covmat'