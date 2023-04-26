#!/bin/bash
#----------------------------------------------------
# Sample Slurm job script

#SBATCH -J sample           # Job name
#SBATCH -o /path/sample.o%j       # Name and Path of stdout output file
#SBATCH -e /path/sample.e%j       # Name and Path of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 1               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 10:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH --mail-user=username@utrgv.edu

module load matlab
matlab -batch  "addpath($SDMBIGDAT19/CODES); cd /working_dir; rungwpso  /path/to/jsonfiles/allparamfiles.json"
