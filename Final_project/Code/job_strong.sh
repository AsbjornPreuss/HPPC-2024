#!/usr/bin/env bash
#SBATCH --job-name=FWC
#SBATCH --partition=modi_HPPC
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=slurm_strong_%j.out
#SBATCH --exclusive

mpiexec apptainer exec \
~/modi_images/ucphhpc/hpc-notebook:latest \
./parallel_hsb --temp $Temperature --nspins $nspins --flips $flips --magnet 10 --writeout 0;


## TODO: CHECK --exclusive, MAKE SURE --magnet 10 IS AN APPROPRIATE VALUE.