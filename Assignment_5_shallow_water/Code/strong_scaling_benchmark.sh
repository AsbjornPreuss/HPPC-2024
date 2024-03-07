#!/bin/bash
# First we do strong scaling. Runs for 20 different grid sizes at constant gang and vector size.
# Strong scaling is defined as how the solution time varies with the number of processors for a fixed total problem size
# Therefore we run at N_size = 4096.
output_file="Strong_scaling.txt";
echo "Version Checksum Elapsed_time grid_size N_gangs vec_length" > $output_file;
N_size=4096;
echo $N_size;
echo "Running strong scaling";
replacement_line="constexpr size_t NX =";
replacement_line+=$N_size;
replacement_line+=", NY = ";
replacement_line+=$N_size;
replacement_line+=";";
sed -i "12s/.*/$replacement_line/" ./sw_sequential.cpp;
sed -i "12s/.*/$replacement_line/" ./sw_parallelv1.cpp;
sed -i "12s/.*/$replacement_line/" ./sw_parallelv2.cpp;
sed -i "12s/.*/$replacement_line/" ./sw_parallelv3.cpp;
# Scan the amount of gangs we are using.
for gang in $(seq 1 5 141); do
    replacement_line="constexpr int number_of_gangs_in_run = ";
    replacement_line+=$gang;
    replacement_line+=";";
    sed -i "13s/.*/$replacement_line/" ./sw_parallelv1.cpp;
    sed -i "13s/.*/$replacement_line/" ./sw_parallelv2.cpp;
    sed -i "13s/.*/$replacement_line/" ./sw_parallelv3.cpp;
    make;
    
    #Run each version 5 times, for statistics
    for k in $(seq 0 1 5); do
        if [ $gang -le 1 ]; then
            ./sw_sequential >> $output_file;
        fi
        ./sw_parallelv1 >> $output_file;
        ./sw_parallelv2 >> $output_file;
        ./sw_parallelv3 >> $output_file;
    done
done





