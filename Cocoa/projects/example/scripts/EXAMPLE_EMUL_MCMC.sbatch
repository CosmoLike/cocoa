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

sleep $(( 10 + SLURM_ARRAY_TASK_ID*20 )) # avoid different scripts running start_cocoa.sh at the same time
cd $SLURM_SUBMIT_DIR
source ~/.bashrc
conda activate cocoa
source start_cocoa.sh

export OMP_NUM_THREADS=1

mpirun -n ${SLURM_NTASKS} --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self \
  --bind-to core:overload-allowed --rank-by slot \
  --map-by slot:pe=${OMP_NUM_THREADS} --mca mpi_yield_when_idle 1 \
  cobaya-run ./projects/example/EXAMPLE_EMUL_MCMC${SLURM_ARRAY_TASK_ID}.yaml -f