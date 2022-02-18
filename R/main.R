################################################################################
# main.R                                              (c) J.M.B. Koch 2022
################################################################################
# This is the main script running the simulation study's
# Dependencies:: cmdstanr; rstan; functions.R; conditions.R; 

# load packages -----------------------------------------------------------
# TBA: make such that packages are installed if not present
library(cmdstanr)
library(rstan)
# source functions --------------------------------------------------------
source('~/1vs2StepBayesianRegSEM/R/functions.R')
# source conditions -------------------------------------------------------
source('~/1vs2StepBayesianRegSEM/R/conditions.R')

# Execute simulation for SVNP ---------------------------------------------
# simulate datasets
datasetsSVNP <- prepareDatasets(condSVNP, nIter, L, Psi, Theta)

# do the sampling
outputFinalSVNP <- sampling(datasetsSVNP, condSVNP, nChain, nWarmup, nSampling)

# write output to .csv
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
