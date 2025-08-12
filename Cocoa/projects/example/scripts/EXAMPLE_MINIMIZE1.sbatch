#!/bin/bash

#SBATCH --job-name=MIN1C
#SBATCH --output=MIN1C-%A.out
#SBATCH --nodes=1
#SBATCH --ntasks=22
#SBATCH --ntasks-per-node=22
#SBATCH --ntasks-per-socket=11
#SBATCH --cpus-per-task=4
#SBATCH --time=120:00:00
#SBATCH --partition=high_priority
#SBATCH --qos=user_qos_timeifler
#SBATCH --account=timeifler
#SBATCH --exclusive

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
source /home/u8/vivianmiranda/miniconda/etc/profile.d/conda.sh
conda activate cocoa
source start_cocoa.sh

export OMP_PROC_BIND=close
if [ -n "$SLURM_CPUS_PER_TASK" ]; then
  export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
else
  export OMP_NUM_THREADS=1
fi

mpirun -n ${SLURM_NTASKS} --oversubscribe --mca pml ^ucx  \
  --mca btl vader,tcp,self --bind-to core:overload-allowed \
  --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS}  \
  python ./projects/example/EXAMPLE_MINIMIZE1.py --nstw 150 \
    --root ./projects/example/ --outroot "example_lcdm_camb_min1"