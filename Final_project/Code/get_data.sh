#!/bin/sh

mkdir -p "Data"

for side in 10 30 100; do # 300 1000 3000; do # MISSING: PICK APPROPRIATE SYSTEM SIZE
	nspins=$(($side * $side * $side))
for temperature in 0.001; do # MISSING: PICK APPROPRIATE TEMPERATURES
    for flip in 100; do # MISSING: PICK APPROPRIATE NUMBER OF FLIPS
	for nranks in 8; do # MISSING: PICK APPROPRIATE NUMBER OF RANKS
		sbatch --ntasks $nranks --export=Temperature=$temperature,nspins=$nspins,flips=$flip ./job_data.sh;
    done
	done
done
done
