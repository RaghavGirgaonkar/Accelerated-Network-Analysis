#!/bin/bash
#----------------------------------------------------
# Sample Slurm job script

#SBATCH -J BNS           # Job name
#SBATCH -o /path/bns.o%j       # Name of stdout output file
#SBATCH -e /path/bns.e%j       # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 1               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 10:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH --mail-user=email@utrgv.edu

module load matlab
matlab -batch  "addpath($SDMBIGDAT19/CODES); rungwpso_bns allparamfiles.json"
