#!/usr/bin/env bash
#SBATCH --job-name=FWC
#SBATCH --partition=modi_HPPC
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -o seq_slurm.out
##SBATCH --exclusive
output_file="seq_out";
#.txt";
#echo "Position_x Position_y Position_z Spin_x Spin_y Spin_z Spin_energy" > $output_file;
echo $Temperature $nspins;
mpiexec apptainer exec \
~/modi_images/ucphhpc/hpc-notebook:latest \
./seq_hsb --temp $Temperature --nspins $nspins --ofile Data/${output_file}_${Temperature}_${nspins}.txt --magnet 0.05;
