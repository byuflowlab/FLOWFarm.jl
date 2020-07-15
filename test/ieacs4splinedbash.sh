#!/bin/bash

#SBATCH --time=01:00:00   # walltime
#SBATCH --ntasks=288   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=1024M   # memory per CPU core
#SBATCH -J "IEACS4 Splined Boundary"   # job name

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
#export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
julia -p $SLURM_CPUS_ON_NODE example_opt_4_ieacs4.jl
