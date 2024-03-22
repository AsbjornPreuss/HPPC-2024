#!/bin/sh
#SBATCH --output=scaling.out

rm -f slurm*

# Strong scaling
for side in 10 30 100 300 1000 3000; do  # MISSING: PICK APPROPRIATE SIZE FOR STRONG SCALING
	nspins=$(($side * $side * $side))
for temperature in 0.001; do #MISSING: SET TEMP TO APPROPRIATE AN VALUE BEFORE BENCHMARKING
    for flip in 100; do #MISSING: PICK APPROPRIATE NUMBER OF FLIPS
	for nranks in 1; do #MISSING: CHOOSE NRANKS SUCH THAT WE GET A LINEAR SCALING BETWEEN WORKLOAD (nspins) AND NRANKS
		sbatch --ntasks $nranks --export=Temperature=$temperature,nspins=$nspins,flips=$flip ./job_strong.sh;
	done
    done
done
done

# Weak scaling
for side in 10 30 100 300 1000 3000; do
	nspins=$(($side * $side * $side))
for temperature in 0.001; do #MISSING: SET TEMP TO APPROPRIATE AN VALUE BEFORE BENCHMARKING
    for flip in 100; do #MISSING: PICK APPROPRIATE NUMBER OF FLIPS
	for nranks in 1; do #MISSING: CHOOSE NRANKS SUCH THAT WE GET A LINEAR SCALING BETWEEN WORKLOAD (nspins) AND NRANKS
		sbatch --ntasks $nranks --export=Temperature=$temperature,nspins=$nspins,flips=$flip ./job_weak.sh;
    done
	done
done
done