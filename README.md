# Getting A Step Ahead: Using the Regularized Horseshoe Prior to Select Cross-Loadings in Bayesian Regularized Structural Equation Modeling (SEM)


This repository contains the code of my Master's Thesis in Methodlogy & Statistics at Utrecht University (September '21 - June '22). 

You can clone the repository by running:

`git clone https://github.com/JMBKoch/1vs2StepBayesianRegSEM/`.

- Note that all scripts assume that this repository have been cloned to the home directory of a unix-based system. Hence, if you're on Windows or you want to work from a different path, you will have to adjust the paths in all scripts manually.

- The simulation study can be conducted by sourcing or running `./R/main.R`. 

- Packages should be installed automatically if they are not yet. However, this may not work on all systems/ versions of R. Hence, if the script does not run checking if the packages are installed correctly may be a sensible first step in the debugging process. An overview of the required packages can be found at the top (line 7-l3) of `./R/parameters.R`.

- The model and all other study parameters (e.g. number of iterations) are specified in `.R/parameters.R`. Note that if the model is adjusted, the code in `./stan` needs to be adjusted accordingly as well. 

- `./R/functions` contains all functions that are used in `.R/main.R`.

- `./data` contains the raw datasets that were simulated based on the population conditions. It will be simulated and saved again when running `.R/main.R`.

- `.Rmd/` contains all manuscripts, the most important one being `.Rmd/thesis/thesis.Rmd`.

(c) J.M.B. Koch, 2022
