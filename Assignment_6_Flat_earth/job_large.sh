#!/usr/bin/env bash
#SBATCH --job-name=FWC
#SBATCH --partition=modi_short
#SBATCH --nodes=1
#SBATCH --ntasks=8
##SBATCH --exclusive

mpiexec apptainer exec \
   ~/modi_images/ucphhpc/hpc-notebook:latest \
   ./fwc_parallel --iter 1000 --model models/large.hdf5