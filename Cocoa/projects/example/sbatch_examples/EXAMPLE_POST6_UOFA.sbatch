#!/bin/bash

#SBATCH --job-name=SO_POST1V2
#SBATCH --output=SO_POST1V2-%j.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --ntasks-per-socket=4
#SBATCH --cpus-per-task=10
#SBATCH --partition=standard
#SBATCH --account=cosmolike
#SBATCH --time=120:00:00
#SBATCH --exclusive 

# Clear the environment from any previously loaded modules
module purge > /dev/null 2>&1
module load anaconda
conda init bash
source ~/.bashrc 

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo Slurm job NAME is $SLURM_JOB_NAME
echo Slurm job ID is $SLURM_JOBID

cd $SLURM_SUBMIT_DIR
conda activate cocoa
source start_cocoa

export OMP_PROC_BIND=close
if [ -n "$SLURM_CPUS_PER_TASK" ]; then
  export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
else
  export OMP_NUM_THREADS=1
fi

mpirun -n ${SLURM_NTASKS} --oversubscribe --report-bindings --mca btl tcp,self --bind-to core --map-by numa:pe=${OMP_NUM_THREADS} cobaya-run ./projects/so/EXAMPLE_POST1V2.yaml -f