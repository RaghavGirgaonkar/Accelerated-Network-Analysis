#!/bin/bash
#----------------------------------------------------
# Sample Slurm job script

#SBATCH -J tauog           # Job name
#SBATCH -o /scratch/09197/raghav/PSD-Estimation/tauog.o%j       # Name of stdout output file
#SBATCH -e /scratch/09197/raghav/PSD-Estimation/tauog.e%j       # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 1               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 18:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH --mail-user=raghav.girgaonkar01@utrgv.edu

module load matlab
matlab -batch  "cd /work/09197/raghav/ls6/MFComparison-Testing_iMac; rungwpso  /work/09197/raghav/ls6/MFComparison-Testing_iMac/allparamfiles_og.json"
