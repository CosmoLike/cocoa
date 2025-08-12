#!/bin/bash

#SBATCH --job-name=MCMC
#SBATCH --output=MCMC-%A-%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --ntasks-per-socket=4   # trying to reserve the socket that has 48 cores
#SBATCH --cpus-per-task=12      # trying to reserve the socket that has 48 cores
#SBATCH --time=110:00:00
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

mpirun -n ${SLURM_NTASKS} --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self \
  --bind-to core:overload-allowed --mca mpi_yield_when_idle 1 \
  --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS} \
  cobaya-run ./projects/example/EXAMPLE_MCMC${SLURM_ARRAY_TASK_ID}.yaml -r