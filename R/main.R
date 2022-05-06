################################################################################
# main.R                                              (c) J.M.B. Koch 2022
################################################################################
# This is the main script running the simulation 
# Dependencies: functions.R; parameters.R; 

# setting seed for reproducibility ----------------------------------------
set.seed(0704)

# source functions and conditions outside clusters ------------------------
source('~/ResearchReport/R/functions.R')
source('~/ResearchReport/R/parameters.R')

# Prepare data ------------------------------------------------------------
# simulate data
#datasets <- simDatasets(condPop = condPop, modelPars = modelPars, nIter = nIter)
# save raw data
#save(datasets, file = "~/ResearchReport/data/datasets.Rds")
# prepare data for stan
#dataStanSVNP <- prepareDat(datasets, condSVNP, nIter)
# save stan-ready data
#save(dataStanSVNP, file = "~/ResearchReport/data/dataStanSVNP.Rds")
# load stan-ready data generally
load("~/ResearchReport/data/dataStanSVNP.Rds")

## Execute simulation for SVNP ---------------------------------------------
### Make sure output is deleted such that no appendation takes place (BE CAREFUL!)
#file.remove(c("~/ResearchReport/output/resultsSVNP.csv",
#              "~/ResearchReport/output/convSVNP.csv"))

# do the sampling where every available core (nWorkers in condtions.R) does 
#    one unique combination of conditions
# measure start time
startTimeSVNP <- Sys.time()
# create clusters
clusters <- makePSOCKcluster(nClusters) 
# source functions & parameters within clusters
clusterCall(clusters, 
            function() source('~/ResearchReport/R/functions.R'))
clusterCall(clusters, 
            function() source('~/ResearchReport/R/parameters.R'))
# Load packages per cluster
clusterCall(clusters, 
            function() lapply(packages, library, character.only = TRUE))
# read in stan-ready data within clusters
clusterCall(clusters,
            function() load("~/ResearchReport/data/dataStanSVNP.Rds"))

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
elapsedTimesSVNP <- endTimeSVNP-startTimeSVNP
