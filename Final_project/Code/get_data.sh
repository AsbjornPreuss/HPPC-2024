#!/bin/sh

mkdir -p "Data"

for side in 9 30 99; do # 300 1000 3000; do # MISSING: PICK APPROPRIATE SYSTEM SIZE
	nspins=$(($side * $side * $side))
for temperature in $(seq 0.001 0.05 0.5); do
	for flip in $(($side * $side * $side * 1000)); do # MISSING: PICK APPROPRIATE NUMBER OF FLIPS
		sed -i "10s/.*/Temperature=${temperature}/" ./job_data.sh
		sed -i "11s/.*/nspins=${nspins}/" ./job_data.sh
		sed -i "12s/.*/flips=${flip}/" ./job_data.sh
		sbatch ./job_data.sh;
done
done
done
