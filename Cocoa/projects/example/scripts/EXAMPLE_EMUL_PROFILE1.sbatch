#!/bin/bash

#SBATCH --job-name=PROF1
#SBATCH --output=PROF1-%A-%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=14
#SBATCH --ntasks-per-node=14
#SBATCH --ntasks-per-socket=14
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
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
conda activate cocoa
source start_cocoa.sh

export OMP_PROC_BIND=close
if [ -n "$SLURM_CPUS_PER_TASK" ]; then
  export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
else
  export OMP_NUM_THREADS=1
fi

ulimit -u 2000000 # require line when nmpi is high
mpirun -n ${SLURM_NTASKS} --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self \
       --bind-to core --map-by core --report-bindings --mca mpi_yield_when_idle 1 \
      python ./projects/example/EXAMPLE_EMUL_PROFILE1.py \
      --root ./projects/example/ --cov 'chains/EXAMPLE_EMUL_MCMC1.covmat' \
      --outroot "EXAMPLE_EMUL_PROFILE1" --factor 3 --maxfeval 25000 --numpts 10 \
      --profile ${SLURM_ARRAY_TASK_ID} \
      --minfile="./projects/example/chains/EXAMPLE_EMUL_MIN1.txt"