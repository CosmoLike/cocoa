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

sleep $(( RANDOM%61 + SLURM_ARRAY_TASK_ID*10 )) # avoid different scripts running start_cocoa.sh 
                                                # at the same time
cd $SLURM_SUBMIT_DIR
source ~/.bashrc
conda activate cocoa
source start_cocoa.sh

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