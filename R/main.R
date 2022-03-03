################################################################################
# main.R                                              (c) J.M.B. Koch 2022
################################################################################
# This is the main script running the simulation study's
# Dependencies functions.R; conditions.R; 
# set.seed(0704)
# source functions and conditions outside clusters ------------------------
source('~/1vs2StepBayesianRegSEM/R/functions.R')
source('~/1vs2StepBayesianRegSEM/R/conditions.R')

# Execute simulation for SVNP ---------------------------------------------
# do the sampling where every available core (nWorkers in condtions.R) does 
#    one unique combination of conditions
# create clusters
clusters <- makePSOCKcluster(nWorkers) 
# source functions & conditions within clusters
clusterCall(clusters, 
            function() source('~/1vs2StepBayesianRegSEM/R/functions.R'))
clusterCall(clusters, 
            function() source('~/1vs2StepBayesianRegSEM/R/conditions.R'))
# make sure packages are loaded per cluster
clusterCall(clusters, 
            function() lapply(packages, library, character.only = TRUE))
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
## source functions & conditions within clusters
#clusterCall(clusters, function() source('~/1vs2StepBayesianRegSEM/R/functions.R'))
#clusterCall(clusters, function() source('~/1vs2StepBayesianRegSEM/R/conditions.R'))
# make sure packages are loaded per cluster
# clusterCall(clusters, function() lapply(packages, library, character.only = TRUE))
# run function in clustered way where every set of condition gets it's own core
#outputFinalRHSP <- clusterApplyLB(clusters, 
#                                  1:nrow(condRHSP), 
#                                  sampling,
#                                  conditions = condRHSP,
#                                  datasets = condRHSP,
#                                  nIter = nIter,
#                                  nChain = nChain,
#                                  nWarmup = nWarmup,
#                                  nSampling = nSampling)
## close clusters
#stopCluster(clusters) 
