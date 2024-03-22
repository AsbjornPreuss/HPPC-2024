#!/usr/bin/env bash
#SBATCH --job-name=HotMountainD
#SBATCH --partition=modi_HPPC
#SBATCH --nodes=1
#SBATCH --ntasks=27
#SBATCH -o data_slurm.out
##SBATCH --exclusive
output_file="data";

Temperature=0.451
nspins=1000000
flips=10000000

echo ./parallel_hsb --temp $Temperature --nspins $nspins --flips $flips;
mpiexec apptainer exec \
~/modi_images/ucphhpc/hpc-notebook:latest \
./parallel_hsb --temp $Temperature --nspins $nspins --flips $flips --ofile Data/${output_file}_${Temperature}_${nspins}.txt --magnet 0.05; 
