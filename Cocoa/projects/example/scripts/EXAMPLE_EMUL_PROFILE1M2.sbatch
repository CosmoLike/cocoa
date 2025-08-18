#!/bin/bash

#SBATCH --job-name=PROF1M2
#SBATCH --output=PROF1-%A.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks-per-socket=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00
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

export OMP_NUM_THREADS=1

for ((i=0; i<7; i++)); do
mpirun -n 1 --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self --rank-by slot \
    --bind-to core:overload-allowed --map-by slot --mca mpi_yield_when_idle 1 \
    python ./projects/example/EXAMPLE_EMUL_PROFILE_SCIPY1.py \
        --root ./projects/example/ --cov 'chains/EXAMPLE_EMUL_MCMC1.covmat' \
        --outroot "EXAMPLE_EMUL_PROFILE1M2" --factor 3 --maxfeval 5000 --numpts 10 \
        --profile $i --minfile="./projects/example/chains/EXAMPLE_EMUL_MIN1.txt"
done