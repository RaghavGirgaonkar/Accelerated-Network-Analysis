# Accelerated-Network-Analysis

## Running Chirp-time Space PSO on GW170817 Data

The ``` examples ``` folder contains the L1 datafile for the GW170817. The following settings could be changed,
- Change the number of PSO runs ``` nruns ``` and iterations  ``` maxSteps ``` in the ``` pso.json ``` file.
- Inject your own CBC signal in the ``` rungwpso.m ``` script

The script can be run in MATLAB by 
``` rungwpso allparamfiles.json ```
