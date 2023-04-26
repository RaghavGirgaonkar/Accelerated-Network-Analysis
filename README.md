# Accelerated-Network-Analysis

## Installation Requirements

- Install/clone [SDMBIGDAT19](https://github.com/mohanty-sd/SDMBIGDAT19)
- If SDMBIGDAT19 is installed at ```$SDMBIGDAT19```, add the following path to MATLAB
``` addpath($SDMBIGDAT19/CODES); ``` 

## Running the PSO Script

Chirp-time or Mass PSO can be ran using the ```rungwpso.m``` script. Before running, ensure following steps are complete:
- Change value of the sampling frequency ```sampling_freq``` and data segment length in seconds ```T_sig_len``` in ```signal.json``` accordingly.
- If injecting custom CBC signal, injection parameters can also be changed in ```signal.json```.
- In ```files.json``` specify the full paths of the data file ```datafile``` and the outplot struct ```output_struct_location``` and plot names ```bestfitplot```
  ```psoresultplot``` and ```bestlocplot```.
- If using custom PSD file provide corresponding file path ```psdfile``` in ```files.json```. If creating colored noise from custom noise realizations
 specify full path for the noise file ```noisefile```.
- Change PSO type ```type``` in ```pso.json``` to either ```tau``` or ```mass``` for Chirp-time and Mass PSO respectively.
- Change the number of PSO iterations ```maxSteps``` and/or number of independent PSO runs ```nruns``` in ```pso.json```.
- If ```signal.json```, ```files.json``` and ```pso.json``` are in different directories, ensure their correct paths are entered in ```allparamfiles.json```

To run PSO in MATLAB on local machine use the following command
```rungwpso <path>/allparamfiles.json```

For example, if ```rungwpso.m``` and ```allparamfiles.json``` are in the same directory,
```rungwpso allparamfiles.json```

The output of the script would include
- PSO estimates of parameters, displayed in command line.
- If requested, plots for the Best Fit Location evolution ```bestlocplot```, Best Fit Fitness Evolution ```bestfitplot``` and PSO-estimated signal ```psoresultplot```.
- Output Struct returned by PSO at ```output_struct_location```.

Sample Command Line Output:
``` 
Original parameters: tau0= 29.6373; tau1p5= 1.1045; m1= 2; m2= 2; A = 8; phi = 0; t_a = 57; FitVal = 314820591.1893
Estimated parameters: tau0=29.6225; tau1p5=1.0967; m1= 1.9868-0.21416i; m2= 1.9868+0.21416i; A = 8.3707; phi = 0.19236; t_a = 57.0112; FitVal = 321155053.9049 
``` 

## Running PSO on Ls6

The PSO script can be run on Ls6 using the sample slurm script ```sample.slurm``` given in ```OptimizedPSOCodes```. 
Before running,
- Change Job Name ```#SBATCH -J <jobname>```
- Change Paths to stdout and stderr files
  ```
  #SBATCH -o /path/sample.o%j       
  #SBATCH -e /path/sample.e%j
  ```
- Provide your email address for job start and end alerts ```#SBATCH --mail-user=username@utrgv.edu```
- Change running time ```#SBATCH -t HH:MM:SS```, typically for 500 PSO iterations on a data length of 4096 seconds and sampling frequency 4096 Hz, PSO takes ~ 7 hours. In this case running time can be set to ```#SBATCH -t 10:00:00```
- Provide paths to SDMBIGDAT19 (on Ls6), ```rungwpso.m``` and ```allparamfiles.json``` as specified in the script

Job to run PSO can then be submitted from the command line as follows
```$> sbatch sample.slurm```

PSO Command line output will now be found in ```/path/sample.o%j```, (```%j``` denotes job number)

## Running Chirp-time Space PSO on GW170817 Data

The ``` examples ``` folder contains the L1 datafile for the GW170817. The following settings could be changed,
- Change the number of PSO runs ``` nruns ``` and iterations  ``` maxSteps ``` in the ``` pso.json ``` file.
- Inject your own CBC signal in the ``` rungwpso.m ``` script

The script can be run in MATLAB by 
``` rungwpso allparamfiles.json ```
