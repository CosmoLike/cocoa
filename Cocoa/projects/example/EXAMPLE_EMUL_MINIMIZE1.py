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
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
parser = argparse.ArgumentParser(prog='EXAMPLE_MINIMIZE1')
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
                    const=1) # zero or one
parser.add_argument("--nwalkers",
                    dest="nwalkers",
                    help="Number of emcee walkers",
                    nargs='?',
                    const=1)
# need to use parse_known_args because of mpifuture 
args, unknown = parser.parse_known_args()
maxfeval    = args.maxfeval
oroot       = args.root + "chains/" + args.outroot
nwalkers    = args.nwalkers
cov_file = args.root + args.cov
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
model = get_model(yaml_load("./projects/example/EXAMPLE_EMUL_EVALUATE1.yaml"))
name  = list(model.parameterization.sampled_params().keys())
def chi2(p):
    point = dict(zip(model.parameterization.sampled_params(),
                 model.prior.sample(ignore_external=True)[0]))
    names = list(model.parameterization.sampled_params().keys())
    point.update({name: val for name, val in zip(names, p)})
    res1 = model.logprior(point,make_finite=True)
    res2 = model.loglike(point,make_finite=True,cached=False,return_derived=False)
    return -2.0*(res1+res2)
def chi2v2(p):
    point = dict(zip(model.parameterization.sampled_params(),
                 model.prior.sample(ignore_external=True)[0]))
    names=list(model.parameterization.sampled_params().keys())
    point.update({name: val for name, val in zip(names, p)})
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
             bounds,
             cov, 
             fixed=-1, 
             maxfeval=3000,
             nwalkers=5,
             pool=None):
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
        
        cov = np.delete(cov, (fixed), axis=0)
        cov = np.delete(cov, (fixed), axis=1)
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
    
    class GaussianStep:
       def __init__(self, stepsize=0.2):
           self.cov = stepsize*cov
       def __call__(self, x):
           return np.random.multivariate_normal(x, self.cov, size=1)   
    
    ndim        = int(x0.shape[0])
    nwalkers    = int(nwalkers)
    nsteps      = maxfeval
    temperature = np.array([1.0, 0.25, 0.1, 0.005, 0.001], dtype='float64')
    stepsz      = temperature/4.0

    mychi2(x0, *args) # first call takes a lot longer (when running on cuda)

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
                                        moves=[(emcee.moves.GaussianMove(cov=GScov),1.)],
                                        pool=pool)
        
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

def prf(x0, index, maxfeval, bounds, cov, nwalkers=5, pool=None):
    t0 = np.array(x0, dtype='float64')
    t1 = np.array(bounds, dtype="float64") # np.array do a deep copy. Deep copy necessary 
                                           # line to avoid weird bug that changes on bounds
                                           # propagate from different iterations (same MPI core)
    res =  min_chi2(x0=t0, 
                    bounds=t1, 
                    fixed=index, 
                    maxfeval=maxfeval, 
                    nwalkers=nwalkers, 
                    pool=pool,
                    cov=cov)
    return res
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if __name__ == '__main__':
    with MPIPool() as pool:
        if not pool.is_master():
            pool.wait()
            sys.exit(0)
        print(f"nwalkers={nwalkers}, maxfeval={maxfeval}")
        (x, results) = model.get_valid_point(max_tries=10000, 
                                             ignore_fixed_ref=False,
                                             logposterior_as_dict=True)
        bounds0 = model.prior.bounds(confidence=0.999999)
        cov = np.loadtxt(cov_file)[0:model.prior.d(),0:model.prior.d()]
        res = np.array(list(prf(np.array(x, dtype='float64'), 
                               index=-1, 
                               maxfeval=maxfeval, 
                               bounds=bounds0, 
                               nwalkers=nwalkers,
                               pool=pool,
                               cov=cov)), dtype="object")
        x0 = np.array([res[0]],dtype='float64')
        # Append derived (begins) ----------------------------------------------
        x0 = np.column_stack((x0, 
                              np.array([chi2v2(d) for d in x0],dtype='float64'),
                              res[1]))
        # Append derived (ends) ------------------------------------------------
        # --- saving file begins -----------------------------------------------
        names = list(model.parameterization.sampled_params().keys()) # Cobaya Call
        names = names+list(model.info()['likelihood'].keys())+["prior"]+["chi2"]
        rnd = random.randint(0,1000)
        print("Output file = ", oroot + "_" + str(rnd) + ".txt")
        np.savetxt(oroot + "_" + str(rnd) +".txt", 
                   x0,
                   fmt="%.6e",
                   header=f"nwalkers={nwalkers}, maxfeval={maxfeval}\n"+' '.join(names),
                   comments="# ")
        # --- saving file ends -------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
#HOW TO CALL THIS SCRIPT
#export nmpi=5
#mpirun -n ${nmpi} --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self \
#  --bind-to core:overload-allowed --mca mpi_yield_when_idle 1 \
#  --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS} \
#  python ./projects/example/EXAMPLE_EMUL_MINIMIZE1.py --root ./projects/example/ \
#  --cov 'EXAMPLE_EMUL_MCMC1.covmat' --outroot "example_min1" --nwalkers 5 --maxfeval 7000