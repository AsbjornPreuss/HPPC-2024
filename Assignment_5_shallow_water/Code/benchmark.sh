#!/bin/bash

for i in $(seq 1000 1000 8000); do
	echo $i;
	./sw_parallel --iter $i;
done



