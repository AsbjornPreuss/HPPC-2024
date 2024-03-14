#!/usr/bin/env bash

for Ntasks in 1 2 4 8 32 56 63 126 252; do
    echo $Ntasks;
    Nnodes=$(($Ntasks / 64 + 1));
    for i in $(seq 0 1 8); do
        sbatch --ntasks=$Ntasks --nodes=$Nnodes job_small.sh;
    done
done

for Ntasks in 1 2 4 6 8 10 12 16 20 24 32 36 40 48 60 64 72 80 96 100 120 144 160 180 192 200 240; do 
    echo $Ntasks;
    Nnodes=$(($Ntasks / 64 + 1));
    for i in $(seq 0 1 8); do
        sbatch --ntasks=$Ntasks --nodes=$Nnodes job_medium.sh;
    done
done

for Ntasks in 1 2 4 6 8 10 12 16 20 24 32 36 40 48 64 72 80 96 120 128 144 160 192 200 240 256; do 
    echo $Ntasks;
    Nnodes=$(($Ntasks / 64 + 1));
    for i in $(seq 0 1 8); do
        sbatch --ntasks=$Ntasks --nodes=$Nnodes job_large.sh;
    done
done

