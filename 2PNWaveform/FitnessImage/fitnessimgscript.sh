#!/bin/bash

#Usage
function usage {
	echo "This is a bash script to create a fitness image for a specified signal in a specified Tau Space range"
	echo " "
	echo "Usage: ./bash fitnessimgscript.sh -n <Gridpoints per dimension (eg 2000)> -s <Total number of sectors (please make sure this a square value)> -m <Mass 1> -M <Mass 2> -a <Starting val of Tau0 range> -b <ending val of tau0 range> -c <Starting val of tau1.5 range> -d <Ending val of tau1.5 range>"
	echo " "
	echo "n: Number of gridpoints per dimension"
	echo "s: Total number of sectors to divide data into"
	echo "m: Value of Mass 1 in Solar Masses"
	echo "M: Value of Mass 2 in Solar Masses"
	echo "a: Starting val of tau0 space"
	echo "b: Ending val of tau0 space range"
	echo "c: Start of range of tau1.5 space"
	echo "d: End range of Tau1.5 space"
	echo "h: Display this help message"
	echo " "
	exit 1
}

#Define list of arguments expected in the input 
optstring=":n:s:m:M:a:b:c:d:h"

while getopts ${optstring} arg; do
  case ${arg} in
    n) 
       N=${OPTARG}
       #echo "${LTA}" 
       echo "Number of gridpoints in each dimension is $N"
       ;;
    s)
       S=${OPTARG}
       echo "Total number of sectors is ${S}" 
       ;;
    m)
       m=${OPTARG}
       echo "Mass 1 is ${m}" 
       ;;
    M)
       M=${OPTARG}
       echo "Mass 2 is ${M}"
       ;;
    a)
       a=${OPTARG}
      # echo "Mass 1 is ${m}"
       ;;
    b)
       b=${OPTARG}
       #echo "Mass 2 is ${M}"
       ;;
    c)
       c=${OPTARG}
       #echo "Mass 1 is ${m}"
       ;;
    d)
       d=${OPTARG}
       #echo "Mass 2 is ${M}"
       ;;
    h)
       usage
       ;; 
    :)
       echo "Invalid option: $OPTARG" 1>&2
       usage
       exit 1
       ;;

    \?)
      echo "Invalid option: -${OPTARG}."
      usage
      ;;
  esac
done

#Update masses and Fitness Image Parameters
sed -i '/masses/c\"masses\": ['"$m"','"$M"'],' signal.json
sed -i '/ngridpoints/c\"ngridpoints\": '"$N"',' signal.json
sed -i '/nsectors/c\"nsectors\": '"$S"',' signal.json
sed -i '/tau0_range/c\"tau0_range\": ['"$a"','"$b"'],' signal.json
sed -i '/tau1p5_range/c\"tau1p5_range\": ['"$c"','"$d"']' signal.json

#Get size of Sector
nindex=`echo "$N/sqrt($S)" | bc`

for i in $(seq 1 ${nindex} ${N})
do
	for j in $(seq 1 ${nindex} ${N})
	do
		sed -i '/startidx_tau0/c\"startidx_tau0\": '"$i"',' signal.json
		sed -i '/startidx_tau1p5/c\"startidx_tau1p5\": '"$j"',' signal.json
		sed -i '/fitvalimgfile/c\"fitvalimgfile\": "'"$SCRATCH"'/files'"$i"''"$j"'.mat"' files.json
		cp signal.json signal$i$j.json
		cp files.json files$i$j.json
		sed -i '/signalparamfile/c\"signalparamfile\": "'"$WORK"'/signal'"$i"''"$j"'.json",' allparamfiles.json
		sed -i '/filenames/c\"filenames\": "'"$WORK"'/files'"$i"''"$j"'.json"' allparamfiles.json
		cp allparamfiles.json allparamfiles$i$j.json
	done
done

chmod +x create_slurmfiles.sh
chmod +x runslurm.sh
#Create Slurm Files
./create_slurmfiles.sh $N $S

#Run Slurm Files
./runslurm.sh $N $S
