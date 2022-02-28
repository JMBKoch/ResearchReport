################################################################################
# main.R                                              (c) J.M.B. Koch 2022
################################################################################
# This is the main script running the simulation study's
# Dependencies functions.R; conditions.R; see packages below

# load packages -----------------------------------------------------------
# specify packages that are required for executing the simulation
packages <- c("cmdstanr", # MCMC sampling using stan
              "rstan", # postprocessing of samples
              "tidyverse", # data wrangling, plotting, pipes
              "mvtnorm", # data simulation
              "parallel" # parallelization
              )

# make sure that packages are installed if not present & load packages
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

# source functions --------------------------------------------------------
source('~/1vs2StepBayesianRegSEM/R/functions.R')
# source conditions -------------------------------------------------------
source('~/1vs2StepBayesianRegSEM/R/conditions.R')

# Execute simulation for SVNP ---------------------------------------------
# simulate datasets
datasetsSVNP <- prepareDatasets(condSVNP, nIter, L, Psi, Theta)

# do the sampling
outputFinalSVNP <- sampling(datasetsSVNP, condSVNP, nChain, nWarmup, nSampling)

# write output to .csv ### TBA: adjust to be part of sampling function and append output
write.csv(outputFinalSVNP$results, 
          file = "~/1vs2StepBayesianRegSEM/output/ResultsMiniSimSVNP.csv")
write.csv(outputFinalSVNP$convergence,
          file = "~/1vs2StepBayesianRegSEM/output/ConvergenceMiniSimSVNP.csv")

# Execute simulation for RHSP ---------------------------------------------
# simulate datasets
#datasetsRHSP <- prepareDatasets(condRHSP, nIter, L, Psi, Theta)
## do the sampling
#outputFinalRHSP <- sampling(datasetsRHSP, condRHSP, nIter, L, Psi, Theta)
#  
## write output to .csv
#write.csv(outputFinalSVNP$results, 
#          file = "~/1vs2StepBayesianRegSEM/output/ResultsMiniSimRHSP.csv")
#write.csv(outputFinalSVNP$convergence, 
#          file = "~/1vs2StepBayesianRegSEM/output/ConvergenceMiniSimRHSP.csv")


