#!/bin/sh

for nspins in 1000; do
for temperature in 0.001; do
	for nranks in 8; do
		sbatch --ntasks $nranks --export=Temperature=$temperature,nspins=$nspins ./job.sh;
	done
done
done
