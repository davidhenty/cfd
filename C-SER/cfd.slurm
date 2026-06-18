#!/bin/bash

# Slurm job options (name, compute nodes, job time)
#SBATCH --job-name=sharpen
#SBATCH --time=0:5:0
#SBATCH --ntasks=1
# You should normally set this to your own budget, e.g. m25oc-s1234567
#SBATCH --cpus-per-task=1
#SBATCH --partition=standard
#SBATCH --qos=short

# Set the number of threads to 1
#   This prevents any threaded system libraries from automatically
#   using threading.
export OMP_NUM_THREADS=1

export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

# Launch the serial job
#   Using 1 thread
srun --unbuffered --distribution=block:block --hint=nomultithread ./cfd 8 20000
