#!/bin/bash

# Slurm job options (name, compute nodes, job time)
#SBATCH --job-name=cfd
#SBATCH --output=%x-%j.out
#SBATCH --time=0:05:0
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --tasks-per-node=4
#SBATCH --cpus-per-task=1
# You should normally set this to your own budget, e.g. m25oc-s1234567
#SBATCH --partition=standard
#SBATCH --qos=short

# Set the number of threads to 1
#   This prevents any threaded system libraries from automatically
#   using threading.
export OMP_NUM_THREADS=1

export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

# Launch the parallel job
#   srun picks up the distribution from the sbatch options
srun --unbuffered --distribution=block:block --hint=nomultithread ./cfd 8 20000
