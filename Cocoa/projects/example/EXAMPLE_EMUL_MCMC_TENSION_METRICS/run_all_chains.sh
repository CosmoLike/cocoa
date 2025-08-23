#!/bin/bash

declare -a CHAINS=("EXAMPLE_EMUL_MCMC1" 
                   "PANTHEON/EXAMPLE_EMUL_MCMC1" 
                   "PANTHEON/EXAMPLE_EMUL_MCMC2" 
                   "PANTHEON/EXAMPLE_EMUL_MCMC3" 
                   "PANTHEON/EXAMPLE_EMUL_MCMC4" 
                   "PANTHEONPLUS/EXAMPLE_EMUL_MCMC1" 
                   "PANTHEONPLUS/EXAMPLE_EMUL_MCMC2" 
                   "PANTHEONPLUS/EXAMPLE_EMUL_MCMC3" 
                   "PANTHEONPLUS/EXAMPLE_EMUL_MCMC4" 
                   "DESY5/EXAMPLE_EMUL_MCMC1" 
                   "DESY5/EXAMPLE_EMUL_MCMC2" 
                   "DESY5/EXAMPLE_EMUL_MCMC3" 
                   "DESY5/EXAMPLE_EMUL_MCMC4"
                   "UNION3/EXAMPLE_EMUL_MCMC1" 
                   "UNION3/EXAMPLE_EMUL_MCMC2" 
                   "UNION3/EXAMPLE_EMUL_MCMC3" 
                   "UNION3/EXAMPLE_EMUL_MCMC4")

for chain in "${CHAINS[@]}"; do
  root="./projects/example/EXAMPLE_EMUL_MCMC_TENSION_METRICS"
  mpirun -n 4 --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self \
     --bind-to core --map-by numa --report-bindings \
     cobaya-run ${root}/${chain}.yaml -r
done