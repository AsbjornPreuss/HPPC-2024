#!/usr/bin/bash

echo "Final_energy Elapsed_time Temperature B_field System_size No_of_ranks Version" > collected_weak.txt;
echo "Final_energy Elapsed_time Temperature B_field System_size No_of_ranks Version" > collected_strong.txt;

for f in $(ls slurm_weak*); do
	cat $f | tail -n +2 | tee -a collected_weak.txt >> /dev/null;
done

for f in $(ls slurm_strong*); do
	cat $f | tail -n +2 | tee -a collected_strong.txt >> /dev/null;
done