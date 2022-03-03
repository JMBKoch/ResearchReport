################################################################################
# main.R                                              (c) J.M.B. Koch 2022
################################################################################
# This is the main script running the simulation study's
# Dependencies functions.R; conditions.R; see packages below
# set.seed(0704)

# load packages -----------------------------------------------------------
# specify packages that are required for executing the simulation
packages <- c("cmdstanr", # MCMC sampling using stan
              "rstan", # postprocessing of samples
              "tidyverse", # data wrangling, plotting, pipes
              "mvtnorm", # data simulation
              "parallel", # parallelization
              "posterior", # median of posterios
              "bayesplot" # convergence diagnostics 
              )

# make sure that packages are installed if not present
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      #library(x, character.only = TRUE)
    }
  }
)

# source functions and conditions outside clusters ------------------------
source('~/1vs2StepBayesianRegSEM/R/functions.R')
source('~/1vs2StepBayesianRegSEM/R/conditions.R')

# Execute simulation for SVNP ---------------------------------------------
# do the sampling where every available core (nWorkers in condtions.R) does 
#    one unique combination of conditions
clusters <- makePSOCKcluster(6) # create cluster
# make sure packages are loaded per cluster
clusterCall(clusters, function() library(tidyverse)) ### TBA make nicer
clusterCall(clusters, function() library(rstan))
clusterCall(clusters, function() library(cmdstanr))
clusterCall(clusters, function() library(mvtnorm))
clusterCall(clusters, function() library(parallel))
clusterCall(clusters, function() library(bayesplot))
clusterCall(clusters, function() library(posterior))
#clusterCall(clusters, function() lapply(c("cmdstanr", # MCMC sampling using stan
#                                          "rstan", # postprocessing of samples
#                                          "tidyverse", # data wrangling, plotting, pipes
#                                          "mvtnorm", # data simulation
#                                          "parallel", # parallelization
#                                          "bayesplot" # convergence diagnostics 
#                                          ), require, character.only = TRUE))
#

# source functions & conditions within clusters
clusterCall(clusters, function() source('~/1vs2StepBayesianRegSEM/R/functions.R'))
clusterCall(clusters, function() source('~/1vs2StepBayesianRegSEM/R/conditions.R'))
# run functon in clustered way where every set of condition gets it's own core
outputFinalSVNP <- clusterApplyLB(clusters, 
                                  1:nrow(condSVNP), 
                                  sampling,
                                  conditions = condSVNP,
                                  modelPars = modelPars,
                                  nIter = nIter,
                                  samplePars = samplePars)

# close clusters
stopCluster(clusters) 

# Execute simulation for RHSP ---------------------------------------------
## do the sampling
#clusters <- makePSOCKcluster(nWorkers) # create cluster
#
## make sure packages are loaded per cluster
#clusterCall(clusters, function() library(tidyverse))
#clusterCall(clusters, function() library(rstan))
#clusterCall(clusters, function() library(cmdstanr))
#clusterCall(clusters, function() library(mvtnorm))
#clusterCall(clusters, function() library(parallel))
#clusterCall(clusters, function() library(bayesplot))
#clusterCall(clusters, function() library(posterior))
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




