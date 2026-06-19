#!/bin/bash

# Slurm job options (name, compute nodes, job time)
#SBATCH --job-name=cfd
#SBATCH --output=%x-%j.out
#SBATCH --time=0:5:0
#SBATCH --ntasks=1
# You should normally set this to your own budget, e.g. m25oc-s1234567
#SBATCH --cpus-per-task=128
#SBATCH --partition=standard
#SBATCH --qos=short

# Set the number of threads as appropriate. For MPI programs we
# increase the number of tasks; for OpenMP we run as single process
# (ntasks=1) but create multiple threads at runtime. The number of
# threads is controlled by OMP_NUM_THREADS

export OMP_NUM_THREADS=4

# A few other variables just to ensure everything runs optimally

export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores

# Launch the serial job

srun --unbuffered --distribution=block:block --hint=nomultithread ./cfd 8 20000
