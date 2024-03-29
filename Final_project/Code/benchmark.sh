#!/bin/sh

for side in 108; do
	nspins=$(($side * $side * $side))
	echo $nspins
    for temperature in 0.001; do
    	for flip in 1728000; do
    		sed -i "10s/.*/Temperature=${temperature}/" ./job8.sh
    		sed -i "11s/.*/nspins=${nspins}/" ./job8.sh
    		sed -i "12s/.*/flips=${flip}/" ./job8.sh

            sed -i "10s/.*/Temperature=${temperature}/" ./job27.sh
    		sed -i "11s/.*/nspins=${nspins}/" ./job27.sh
    		sed -i "12s/.*/flips=${flip}/" ./job27.sh
      
            sed -i "10s/.*/Temperature=${temperature}/" ./job64.sh
    		sed -i "11s/.*/nspins=${nspins}/" ./job64.sh
    		sed -i "12s/.*/flips=${flip}/" ./job64.sh

            sed -i "10s/.*/Temperature=${temperature}/" ./job.sh
    		sed -i "11s/.*/nspins=${nspins}/" ./job.sh
    		sed -i "12s/.*/flips=${flip}/" ./job.sh

          	for i in 1 2 3 4 5; do
                sbatch ./job8.sh;
        		sbatch ./job27.sh;
        		sbatch ./job64.sh;
            done
            sbatch ./job.sh;
    	done
    done
done
