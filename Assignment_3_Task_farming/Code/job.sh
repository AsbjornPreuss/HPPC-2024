#!/usr/bin/env bash
#SBATCH --job-name=TaskFarm
#SBATCH --partition=modi_devel
#SBATCH --nodes=1
#SBATCH --ntasks=64
##SBATCH --exclusive

mpiexec apptainer exec \
   ~/modi_images/ucphhpc/hpc-notebook:latest \
   ./task_farm_HEP >> HEP_out_5

