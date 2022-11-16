#!/bin/bash
# --------------------------------------------------------------
### PART 1: Requests resources to run your job.
# --------------------------------------------------------------
### Optional. Set the job name
#SBATCH --job-name=gwpso
### Optional. Set the output filename.
### SLURM reads %x as the job name and %j as the job ID
#SBATCH --output=$SCRATCH/gwpso.out
#SBATCH --error=$SCRATCH/gwpso.err
### REQUIRED. Specify the PI group for this job
#SBATCH --account=manager
### Optional. Request email when job begins and ends
### Specify high priority jobs
#SBATCH --qos=user_qos_manager
# SBATCH --mail-type=ALL
### Optional. Specify email address to use for notification
# SBATCH --mail-user=raghav.girgaonkar01@utrgv.edu
### REQUIRED. Set the partition for your job.
#SBATCH --partition=standard
### REQUIRED. Set the number of cores that will be used for this job.
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
### REQUIRED. Set the memory required for this job.
#SBATCH --mem=3gb
### REQUIRED. Specify the time required for this job, hhh:mm:ss
#SBATCH --time=02:30:00
# --------------------------------------------------------------
### PART 2: Executes bash commands to run your job
# --------------------------------------------------------------
### SLURM Inherits your environment. cd $SLURM_SUBMIT_DIR not needed
pwd; hostname; date
echo "CPUs per task: $SLURM_CPUS_PER_TASK"
### Load required modules/libraries if needed
module load matlab
matlab -nodisplay -nosplash -softwareopengl < $WORK/Accelerated-Network-Analysis/2PNWaveform/rungwpso.m > $SCRATCH/out_gwpso.txt

