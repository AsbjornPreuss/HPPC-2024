#!/bin/sh

for side in 10 30 100 300 1000 3000; do
	nspins=$(($side * $side * $side))
	echo $nspins
for temperature in 0.001; do
	for nranks in 8; do
		sbatch --ntasks $nranks --export=Temperature=$temperature,nspins=$nspins ./job.sh;
	done
done
done
