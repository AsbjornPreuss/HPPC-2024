#!/usr/bin/env bash
#SBATCH --job-name=FWC
#SBATCH --partition=modi_HPPC
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --exclusive
output_file="seq_out";

Temperature=0.451
Temperature=0.001
nspins=27000000
flips=17280000
#echo "Position_x Position_y Position_z Spin_x Spin_y Spin_z Spin_energy" > $output_file;
echo $Temperature $nspins;
mpiexec apptainer exec \
~/modi_images/ucphhpc/hpc-notebook:latest \
./parallel_hsb_ranks8 --temp $Temperature --nspins $nspins --flips $flips --writeout 0 --magnet 0.05;
