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
    message=r".*overflow encountered in exp.*"
)
import argparse, random
import numpy as np
from cobaya.yaml import yaml_load
from cobaya.model import get_model
from nautilus import Prior, Sampler
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
parser = argparse.ArgumentParser(prog='EXAMPLE_NAUTILUS1')
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
args, unknown = parser.parse_known_args()
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
model = get_model(yaml_load("./projects/example/EXAMPLE_EMUL_EVALUATE1.yaml"))
def likelihood(p):
    point = dict(zip(model.parameterization.sampled_params(),
                 model.prior.sample(ignore_external=True)[0]))
    names=list(model.parameterization.sampled_params().keys())
    point.update({ name: p[name].item() for name in names })
    res1 = model.logprior(point,make_finite=True)
    res2 = model.loglike(point,make_finite=True,cached=False,return_derived=False)
    return res1+res2
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
from mpi4py.futures import MPIPoolExecutor
if __name__ == '__main__':
    #rnd = random.randint(0,1000)
    rnd = 100
    cfile= args.root + "chains/" + args.outroot +  "_" + str(rnd) + "_checkpoint" ".hdf5"
    print(f"nlive={args.nlive}, output={cfile}")
    # Here we need to pass Cobaya Prior to Nautilus
    NautilusPrior = Prior()                              # Nautilus Call 
    dim    = model.prior.d()                             # Cobaya call
    bounds = model.prior.bounds(confidence=0.999999)     # Cobaya call
    names  = list(model.parameterization.sampled_params().keys()) # Cobaya Call
    for b, name in zip(bounds, names):
      NautilusPrior.add_parameter(name, dist=(b[0], b[1]))
    # Start Nautilus
    sampler = Sampler(NautilusPrior, 
                      likelihood,  
                      filepath=cfile, 
                      n_dim=dim,
                      pool=MPIPoolExecutor(),
                      n_live=args.nlive,
                      resume=True)
    sampler.run(f_live=args.flive,
                n_eff=args.neff,
                n_like_max=args.maxfeval,
                verbose=True,
                discard_exploration=True)
    points, log_w, log_l = sampler.posterior()
    # --- saving file begins --------------------
    cfile   = args.root + "chains/" + args.outroot +  "_" + str(rnd) + ".1.txt"
    print("Output file = ", cfile)
    np.savetxt(cfile,
               np.column_stack((np.exp(log_w), log_l, points, -2*log_l)),
               fmt="%.5e",
               header=f"nlive={args.nlive}, maxfeval={args.maxfeval}\n"+' '.join(names),
               comments="# ")
    # Now we need to save a paramname files --------------------
    cfile   = args.root + "chains/" + args.outroot +  "_" + str(rnd) + ".paramnames"
    param_info = model.info()['params']
    names2 = names.copy()
    latex  = [param_info[x]['latex'] for x in names2]
    names2.append("chi2*")
    latex.append("\\chi^2")
    np.savetxt(cfile, np.column_stack((names2,latex)),fmt="%s")
    # Now we need to save a range files --------------------
    cfile   = args.root + "chains/" + args.outroot +  "_" + str(rnd) + ".ranges"
    rows = [(str(n),float(l),float(h)) for n,l,h in zip(names,bounds[:,0],bounds[:,1])]
    with open(cfile, "w") as f: 
      f.writelines(f"{n} {l:.5e} {h:.5e}\n" for n, l, h in rows)
    # --- saving file ends --------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
#HOW TO CALL THIS SCRIPT
#mpirun -n 5 --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self \
#   --bind-to core:overload-allowed --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS} \
#   python -m mpi4py.futures ./projects/example/EXAMPLE_EMUL_NAUTILUS1.py \
#   --root ./projects/example/ --outroot "EXAMPLE_NAUTILUS1"  \
#   --maxfeval 1000000 --nlive 500 --neff 10000 --flive 0.01