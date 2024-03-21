#!/usr/bin/env bash
#SBATCH --job-name=FWC
#SBATCH --partition=modi_HPPC
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -o seq_slurm.out
##SBATCH --exclusive
output_file="data";

echo ./seq_hsb --temp $Temperature --ofile ${output_file}_${Temperature}_${nspins}.txt;
mpiexec apptainer exec \
~/modi_images/ucphhpc/hpc-notebook:latest \
./parallel_hsb --temp $Temperature --nspins $nspins --flips $flips --ofile Data/${output_file}_${Temperature}_${nspins}.txt --magnet 10; 
## TODO: MAKE SURE --magnet 10 IS AN APPROPRIATE VALUE.