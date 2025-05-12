#!/bin/bash
#SBATCH --job-name=EP1
#SBATCH --output=EP1-%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=90
#SBATCH --ntasks-per-node=90
#SBATCH --cpus-per-task=1
#SBATCH --time=60:00:00
#SBATCH --partition=standard
#SBATCH --account=cosmolike
#SBATCH --exclusive
#SBATCH --array=0-5

# Clear the environment from any previously loaded modules
module purge > /dev/null 2>&1
source ~/.bashrc 

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo Slurm job NAME is $SLURM_JOB_NAME
echo Slurm job ID is $SLURM_JOBID
echo Number of task is $SLURM_NTASKS
echo Number of cpus per task is $SLURM_CPUS_PER_TASK
echo SLURM_SUBMIT_DIR is $SLURM_SUBMIT_DIR

cd $SLURM_SUBMIT_DIR

conda activate cocoa

source start_cocoa.sh

export OMP_NUM_THREADS=1
export tmp=$((${SLURM_NTASKS}-1))

mpirun -n ${SLURM_NTASKS} --oversubscribe --mca pml ^ucx \
  --mca btl vader,tcp,self --bind-to core:overload-allowed \
  --rank-by slot --map-by core:pe=${OMP_NUM_THREADS} \
  python -m mpi4py.futures ./projects/example/EXAMPLE_PROFILE1.py \
  --tol 0.05 --profile ${SLURM_ARRAY_TASK_ID} --maxiter 5 \
  --maxfeval 10000 --mpi tmp \
  --outroot "profile1" --minmethod 1 --factor 4
