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
		sed -i 's#^matlab .*$#matlab -batch  "cd /work/09197/raghav/ls6/Accelerated-Network-Analysis/2PNWaveform/AllPSO; fitValimg  /work/09197/raghav/ls6/allparamfiles'"$i"''"$j"'.json"#' gwpso.slurm
            	cp gwpso.slurm gwpso$i$j.slurm
        done
done

