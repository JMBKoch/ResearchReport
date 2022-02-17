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

# write output to .csv
write.csv(file = "~/1vs2StepBayesianRegSEM/output/OutputSVNP.csv")

# Execute simulation for RHSP ---------------------------------------------

# write output to .csv
write.csv(file = "~/1vs2StepBayesianRegSEM/output/OutputRHSP.csv")

