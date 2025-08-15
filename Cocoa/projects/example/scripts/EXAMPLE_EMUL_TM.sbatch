#!/bin/bash

#SBATCH --job-name=MCMC
#SBATCH --output=MCMC-%A-%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --ntasks-per-socket=4
#SBATCH --cpus-per-task=1  
#SBATCH --time=12:00:00
#SBATCH --partition=high_priority
#SBATCH --qos=user_qos_timeifler
#SBATCH --account=timeifler


echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo Slurm job NAME is $SLURM_JOB_NAME
echo Slurm job ID is $SLURM_JOBID
echo Number of task is $SLURM_NTASKS
echo Number of cpus per task is $SLURM_CPUS_PER_TASK

cd $SLURM_SUBMIT_DIR

# conda and (.local) envs should be already loaded prior to submitting a slurm script
if [ -z "${ROOTDIR}" ]; then
  sleep $(( 10 + SLURM_ARRAY_TASK_ID*20 )) # help avoid different scripts running start_cocoa.sh simultaneously
                                           # clash can still happen when one script runs start_cocoa.sh, the other
                                           # script is reading a file (start_cocoa.sh invokes stop_cocoa.sh)
                                           # which delete links and ROOTDIR so the running script can find file/likelihood
  source ~/.bashrc
  conda activate cocoa
  source start_cocoa.sh
fi

export OMP_PROC_BIND=close
if [ -n "$SLURM_CPUS_PER_TASK" ]; then
  export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
else
  export OMP_NUM_THREADS=1
fi

mkdir -p "${ROOTDIR:?}/projects/example/chains/TM/PANTHEON"
mkdir -p "${ROOTDIR:?}/projects/example/chains/TM/PANTHEONPLUS"
mkdir -p "${ROOTDIR:?}/projects/example/chains/TM/DESY5"
mkdir -p "${ROOTDIR:?}/projects/example/chains/TM/UNION3"

# MCMC REQUIRED FOR TENSION METRICS ON `EXAMPLE_TENSION_METRICS.ipynb`
declare -a YAML=( "EXAMPLE_EMUL_MCMC1.yaml" 
                  "EXAMPLE_EMUL_MCMC2.yaml" 
                  "PANTHEON/EXAMPLE_EMUL_MCMC1.yaml" 
                  "PANTHEON/EXAMPLE_EMUL_MCMC2.yaml" 
                  "PANTHEON/EXAMPLE_EMUL_MCMC3.yaml" 
                  "PANTHEONPLUS/EXAMPLE_EMUL_MCMC1.yaml" 
                  "PANTHEONPLUS/EXAMPLE_EMUL_MCMC2.yaml" 
                  "PANTHEONPLUS/EXAMPLE_EMUL_MCMC3.yaml" 
                  "DESY5/EXAMPLE_EMUL_MCMC1.yaml" 
                  "DESY5/EXAMPLE_EMUL_MCMC2.yaml" 
                  "DESY5/EXAMPLE_EMUL_MCMC3.yaml" 
                  "UNION3/EXAMPLE_EMUL_MCMC1.yaml" 
                  "UNION3/EXAMPLE_EMUL_MCMC2.yaml" 
                  "UNION3/EXAMPLE_EMUL_MCMC3.yaml" )

mpirun -n ${SLURM_NTASKS} --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self \
  --bind-to core:overload-allowed --rank-by slot \
  --map-by slot:pe=${OMP_NUM_THREADS} --mca mpi_yield_when_idle 1 \
  cobaya-run ./projects/example/EXAMPLE_EMUL_MCMC_TENSION_METRICS/${YAML[${SLURM_ARRAY_TASK_ID}]} -f