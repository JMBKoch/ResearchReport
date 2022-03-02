################################################################################
# main.R                                              (c) J.M.B. Koch 2022
################################################################################
# This is the main script running the simulation study's
# Dependencies functions.R; conditions.R; see packages below
# set.seed(0704)

# source functions --------------------------------------------------------
source('~/1vs2StepBayesianRegSEM/R/functions.R')
# source conditions -------------------------------------------------------
source('~/1vs2StepBayesianRegSEM/R/conditions.R')
### TBA: Think of how to avoid calling them twice

# load packages -----------------------------------------------------------
# specify packages that are required for executing the simulation
packages <- c("cmdstanr", # MCMC sampling using stan
              "rstan", # postprocessing of samples
              "tidyverse", # data wrangling, plotting, pipes
              "mvtnorm", # data simulation
              "parallel", # parallelization
              "bayesplot" # convergence diagnostics 
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

# Execute simulation for SVNP ---------------------------------------------
# simulate datasets
datasetsSVNP <- prepareDatasets(condSVNP, nIter, L, Psi, Theta)
#### Think about putting this into sampling()

# do the sampling where every available core (nWorkers in condtions.R) does 
#    one unique combination of conditions
clusters <- makePSOCKcluster(nWorkers) # create cluster
# make sure packages are loaded per cluster
clusterCall(clusters, function() library(tidyverse))
clusterCall(clusters, function() library(rstan))
clusterCall(clusters, function() library(cmdstanr))
clusterCall(clusters, function() library(mvtnorm))
clusterCall(clusters, function() library(parallel))
clusterCall(clusters, function() library(bayesplot))

# source functions & conditions within clusters
clusterCall(clusters, function() source('~/1vs2StepBayesianRegSEM/R/functions.R'))
clusterCall(clusters, function() source('~/1vs2StepBayesianRegSEM/R/conditions.R'))
# run functon in clustered way where every set of condition gets it's own core
outputFinalSVNP <- clusterApplyLB(clusters, 
                                  1:nrow(condSVNP), 
                                  sampling,
                                  conditions = condSVNP,
                                  datasets = datasetsSVNP,
                                  nIter = nIter,
                                  nChain = nChain,
                                  nWarmup = nWarmup,
                                  nSampling = nSampling)


# close clusters
stopCluster(clusters) 

# write output to .csv ### TBA: adjust to be part of sampling function and append output
#write.csv(outputFinalSVNP$results, 
#          file = "~/1vs2StepBayesianRegSEM/output/ResultsMiniSimSVNP.csv")
#write.csv(outputFinalSVNP$convergence,
#          file = "~/1vs2StepBayesianRegSEM/output/ConvergenceMiniSimSVNP.csv")
#

# Execute simulation for RHSP ---------------------------------------------
# simulate datasets
#datasetsRHSP <- prepareDatasets(condRHSP, nIter, L, Psi, Theta)
## do the sampling
# do the sampling
#clusters <- makePSOCKcluster(nWorkers) # create cluster
#
## make sure packages are loaded per cluster
#clusterCall(clusters, function() library(tidyverse))
#clusterCall(clusters, function() library(rstan))
#clusterCall(clusters, function() library(cmdstanr))
#clusterCall(clusters, function() library(mvtnorm))
#clusterCall(clusters, function() library(parallel))
#clusterCall(clusters, function() library(bayesplot))
#
## source functions & conditions within clusters
#clusterCall(clusters, function() source('~/1vs2StepBayesianRegSEM/R/functions.R'))
#clusterCall(clusters, function() source('~/1vs2StepBayesianRegSEM/R/conditions.R'))

# run functio in clustered way where every set of condition gets it's own core
#outputFinalRHSP <- clusterApplyLB(clusters, 
#                                  1:nrow(condRHSP), 
#                                  sampling,
#                                  conditions = condRHSP,
#                                  datasets = condRHSP,
#                                  nIter = nIter,
#                                  nChain = nChain,
#                                  nWarmup = nWarmup,
#                                  nSampling = nSampling)
#
#
## close clusters
#stopCluster(clusters) 


## write output to .csv
#write.csv(outputFinalSVNP$results, 
#          file = "~/1vs2StepBayesianRegSEM/output/ResultsMiniSimRHSP.csv")
#write.csv(outputFinalSVNP$convergence, 
#          file = "~/1vs2StepBayesianRegSEM/output/ConvergenceMiniSimRHSP.csv")


