#!/usr/bin/env bash

for Ntasks in 1 2 4 9 16 25 36 64 81 100 121 144 169 196 225 256; do
    echo $Ntasks;
    Nnodes=$(($Ntasks / 64 + 1));
    for i in $(seq 0 1 32); do
        sbatch --ntasks=$Ntasks --nodes=$Nnodes job_small.sh;
    done
    for i in $(seq 0 1 16); do
        sbatch --ntasks=$Ntasks --nodes=$Nnodes job_medium.sh;
    done
    for i in $(seq 0 1 8); do
        sbatch --ntasks=$Ntasks --nodes=$Nnodes job_large.sh;
    done
done
