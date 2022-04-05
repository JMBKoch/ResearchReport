################################################################################
# main.R                                              (c) J.M.B. Koch 2022
################################################################################
# This is the main script running the simulation 
# Dependencies: functions.R; parameters.R; 

# setting seed for reproducibility ----------------------------------------
set.seed(0704)

# source functions and conditions outside clusters ------------------------
source('~/1vs2StepBayesianRegSEM/R/functions.R')
source('~/1vs2StepBayesianRegSEM/R/parameters.R')

# Prepare data ------------------------------------------------------------
# simulate data
datasets <- simDatasets(condPop = condPop, modelPars = modelPars, nIter = nIter)
# save raw data
#save(datasets, file = "~/1vs2StepBayesianRegSEM/data/datasets.Rds")
# prepare data for stan
#dataStanSVNP <- prepareDat(datasets, condSVNP, nIter)
# save stan-ready data
#save(dataStanSVNP, file = "~/1vs2StepBayesianRegSEM/data/dataStanSVNP.Rds")
# load stan-ready data generally
load("~/1vs2StepBayesianRegSEM/data/dataStanSVNP.Rds")

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
# source functions & parameters within clusters
clusterCall(clusters, 
            function() source('~/1vs2StepBayesianRegSEM/R/functions.R'))
clusterCall(clusters, 
            function() source('~/1vs2StepBayesianRegSEM/R/parameters.R'))
# Load packages per cluster
clusterCall(clusters, 
            function() lapply(packages, library, character.only = TRUE))
# read in stan-ready data within clusters
clusterCall(clusters,
            function() load("~/1vs2StepBayesianRegSEM/data/dataStanSVNP.Rds"))

# run functon in clustered way where it's clustered over individual combo's of 
#  iteration, condPop and condPrior
outputFinalSVNP <- clusterApplyLB(clusters, 
                                  1:length(dataStanSVNP),
                                  sampling,
                                  dataStan = dataStanSVNP, 
                                  prior = "SVNP",
                                  modelPars = modelPars, 
                                  samplePars = samplePars)
# close clusters
stopCluster(clusters) 
# measure end time
endTimeSVNP <- Sys.time()
# measure elapsed time
endTimeSVNP-startTimeSVNP

# Execute simulation for RHSP ---------------------------------------------
# prepare data for stan
#dataStanRHSP <- prepareDat(datasets, condRHSP, nIter)
# save stan-ready data
#save(dataStanRHSP, file = "~/1vs2StepBayesianRegSEM/data/dataStanRHSP.Rds")

## Execute simulation for RHSP ---------------------------------------------
### Make sure output is deleted such that no appendation takes place (BE CAREFUL!)
file.remove(c("~/1vs2StepBayesianRegSEM/output/resultsRHSP.csv",
              "~/1vs2StepBayesianRegSEM/output/convRHSP.csv"))

# load stan-ready data generally
load("~/1vs2StepBayesianRegSEM/data/dataStanRHSP.Rds")

# do the sampling where every available core (nWorkers in condtions.R) does 
#    one unique combination of conditions
# measure start time
startTimeRHSP <- Sys.time()
# create clusters
clusters <- makePSOCKcluster(nClusters) 
# source functions & parameters within clusters
clusterCall(clusters, 
            function() source('~/1vs2StepBayesianRegSEM/R/functions.R'))
clusterCall(clusters, 
            function() source('~/1vs2StepBayesianRegSEM/R/parameters.R'))
# Load packages per cluster
clusterCall(clusters, 
            function() lapply(packages, library, character.only = TRUE))
# read in stan-ready data within clusters
clusterCall(clusters,
            function() load("~/1vs2StepBayesianRegSEM/data/dataStanRHSP.Rds"))

# run functon in clustered way where it's clustered over individual combo's of 
#  iteration, condPop and condPrior
outputFinalRHSP <- clusterApplyLB(clusters, 
                                  1:length(dataStanRHSP),
                                  sampling,
                                  dataStan = dataStanRHSP, 
                                  prior = "RHSP",
                                  modelPars = modelPars, 
                                  samplePars = samplePars)
# close clusters
stopCluster(clusters) 
# measure end time
endTimeRHSP <- Sys.time()
# measure elapsed time
endTimeRHSP-startTimeRHPS







