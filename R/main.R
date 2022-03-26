################################################################################
# main.R                                              (c) J.M.B. Koch 2022
################################################################################
# This is the main script running the simulation 
# Dependencies: functions.R; parameters.R; 
# set.seed(0704)
# source functions and conditions outside clusters ------------------------
source('~/1vs2StepBayesianRegSEM/R/functions.R')
source('~/1vs2StepBayesianRegSEM/R/parameters.R')

## Execute simulation for SVNP ---------------------------------------------
### Make sure output is deleted such that no appendation takes place (BE CAREFUL!)
file.remove(c("~/1vs2StepBayesianRegSEM/output/resultsSVNP.csv",
              "~/1vs2StepBayesianRegSEM/output/convSVNP.csv"))
# do the sampling where every available core (nWorkers in condtions.R) does 
#    one unique combination of conditions
# measure start time
startTimeSVNP <- Sys.time()
# create clusters
clusters <- makePSOCKcluster(nClusters) 
# source functions & conditions within clusters
clusterCall(clusters, 
            function() source('~/1vs2StepBayesianRegSEM/R/functions.R'))
clusterCall(clusters, 
            function() source('~/1vs2StepBayesianRegSEM/R/parameters.R'))
# Load packages per cluster
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
# measure end time
endTimeSVNP <- Sys.time()
# measure elapsed time
endTimeSVNP-startTimeSVNP

# Execute simulation for RHSP ---------------------------------------------
# do the sampling
# Make sure output is deleted such that no appendation takes place (BE CAREFUL!)
#file.remove(c("~/1vs2StepBayesianRegSEM/output/resultsRHSP.csv",
#              "~/1vs2StepBayesianRegSEM/output/convRHSP.csv"))
#startTimeRHSP <- Sys.time()
## create cluster
#clusters <- makePSOCKcluster(nClusters) 
## source functions & conditions within clusters
#clusterCall(clusters, function() source('~/1vs2StepBayesianRegSEM/R/functions.R'))
#clusterCall(clusters, function() source('~/1vs2StepBayesianRegSEM/R/parameters.R'))
##Load packages per cluster
#clusterCall(clusters, function() lapply(packages, library, character.only = TRUE))
##run function in clustered way where every set of condition gets it's own core
#outputFinalRHSP <- clusterApplyLB(clusters, 
#                                  1:nrow(condRHSP), 
#                                  sampling,
#                                  conditions = condRHSP,
#                                  modelPars = modelPars,
#                                  nIter = nIter,
#                                  samplePars = samplePars)
## close clusters
#stopCluster(clusters) 
## measure end time
#endTimeRHSP <- Sys.time()
## measure elapsed time
#endTimeRHSP-startTimeRHSP
