# AccelerAge simulation study
This repository contains the R-code for the simulation study as presented in the paper [The AccelerAge framework: a new statistical approach to predict biological age based on time-to-event data](https://link.springer.com/article/10.1007/s10654-024-01114-8) (2024). Please see the paper for details on the simulation study set-up. For a simpler, more elaborate explanation of how to fit the AccelerAge framework in R, please see [this repository](https://github.com/marije-sluiskes/fitting-accelerage-framework-in-r) instead. 

## Set-up of repository
The folder [scripts](scripts) contains the files that need to be run to recreate the full simulation study results. Please note that this simulation study was run on a cluster computer. The code in script [2_ObtainEstimates_parallelized.R](scripts/2_ObtainEstimates_parallelized.R) is paralellized. This might cause problems when running the code on a personal computer. 

The scripts assume that the necessary repositories to write the results to have been created prior to running the code. Therefore, when reproducing the results, please first recreate the structure of the [output](output) folder as given in this repository, including all subfolders. 

The folder [missellaneous](misscelaneous) contains two scripts for reproducing Figures 1, 3 and 4 of the paper. 


