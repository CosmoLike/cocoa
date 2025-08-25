#!/bin/bash

#SBATCH --job-name=POLY
#SBATCH --output=POLY-%A-%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=90
#SBATCH --ntasks-per-node=90
#SBATCH --ntasks-per-socket=45 
#SBATCH --cpus-per-task=1  
#SBATCH --time=36:00:00
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

export OMP_NUM_THREADS=1
ulimit -u 2000000 # require line when nmpi is high

mpirun -n ${SLURM_NTASKS} --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self \
  --bind-to core:overload-allowed --rank-by slot \
  --map-by slot:pe=${OMP_NUM_THREADS} --mca mpi_yield_when_idle 1 \
  cobaya-run ./projects/example/EXAMPLE_EMUL_POLY${SLURM_ARRAY_TASK_ID}.yaml -r