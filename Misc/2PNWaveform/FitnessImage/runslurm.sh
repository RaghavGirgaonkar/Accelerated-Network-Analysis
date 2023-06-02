#!/bin/bash

N=$1
S=$2
#echo $N
#echo $S
nindex=`echo "$N/sqrt($S)" | bc`
#echo $nindex
for i in $(seq 1 ${nindex} ${N})
do
        for j in $(seq 1 ${nindex} ${N})
        do
		sbatch gwpso$i$j.slurm
        done
done

