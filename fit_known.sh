#!/bin/bash

#-------------------------------------
# qsub .sh file
#-------------------------------------

for((j=1;j<=100;++j))
do
	for ((i=1;i <=6;++i))
	do 
		Rscript Simulations_known.R "$j" "$i"
	done
done