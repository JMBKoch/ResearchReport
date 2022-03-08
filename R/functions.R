################################################################################
# functions.R                                             (c) J.M.B. Koch 2022
################################################################################
# All functions used in main.R
# dependencies: tidyverse (magrittr, tidyr, dplyr, ggplot2), mvtnorm, bayesplot

################################################################################
# Part 1: functions for executing simulation study
################################################################################
# prepareDataset() -------------------------------------------------------
# function that prepares a stan ready dataset
prepareDataset <- function(conditions, main, Psi, Theta){
  
  # (inner) helper functions ------------------------------------------------
  # simY() 
  # function to simulate data under desired model sourced from parameters.R
  
  ### TBA: work with different levels of cross-loadings!!!
  simY <- function(main, cross, Psi, Theta, N){
    L <- matrix(NA, nrow = 6, ncol = 2)
    L[1:3, 1] <- main[1:3]
    L[4:6, 2] <- main[4:6]
    L[4:6, 1] <- cross[1:3]
    L[1:3, 2] <- cross[4:6]
    Sigma <- L%*%Psi%*%t(L) + Theta
    Y <- mvtnorm::rmvnorm(N, rep(0, 6), Sigma)
    return(Y)
  }
  
  # prepareDat() 
  # function to prepare stan data object from simdat() based on a unique
  #   combination of hyper-pars
  prepareDat <- function(Y, conditions){ 
    if(conditions$prior == "SVNP"){
      out <-  list(
        N = nrow(Y),
        P = ncol(Y),
        Q = 2,
        Y = Y, 
        sigma = conditions$sigma,
        cross = conditions$cross
      )
    }else if(conditions$prior == "RHSP"){
      out <- list(
        N = nrow(Y),
        P = ncol(Y),
        Q = 2,
        Y = Y, 
        scaleGlobal = conditions$scaleGlobal, # scale omega
        scaleLocal = conditions$scaleLocal, # scale lambda
        dfGlobal = conditions$dfGlobal, # df for half-t prior omega
        dfLocal = conditions$dfLocal, # df for half-t prior tau_j
        nu = conditions$nu, # df IG for c^2
        scaleSlab = conditions$scaleSlab, # scale of slab
        cross = conditions$cross
      )
    } 
    return(out)
  }
  
    # allocate memory for output
    dat <- list()
      
    # specify crossloadings & N based on conditions
    if(conditions$cross == 0.2){
      cross <- cross2
    } else{
      cross <- cross5
    }
    N <- conditions$N
    
    # simulate the data
    Y <- simY(main, cross, Psi, Theta, N)
    dat <-  prepareDat(Y, conditions) 

  # return dataset
  return(dat)
}

# saveResults() ------------------------------------------------------------
# function takes rstan object and saves the estimes 
saveResults <- function(rstanObj,  conditions){

  # estimates Lambda ### median straks wss alleen maar belangrijk voor kruisladingen
  mainEstMean <- apply(as.matrix(rstanObj, pars = "lambdaMainC"), 2, mean)
  mainEstMed <-  apply(as.matrix(rstanObj, pars = "lambdaMainC"), 2, median)
  mainEstVar <-  apply(as.matrix(rstanObj, pars = "lambdaMainC"), 2, var)
  crossMatrix <- as.matrix(rstanObj, pars = "lambdaCrossC") # save because handy for quantiles
  crossEstMean <- apply(crossMatrix, 2, mean)
  crossEstMed <-  apply(crossMatrix, 2, median)
  crossEstVar <-  apply(crossMatrix, 2, var)
  
  # estimates Theta
  thetaEstMean <- apply(as.matrix(rstanObj, pars = "theta"), 2, mean)
  thetaEstMed <- apply(as.matrix(rstanObj, pars = "theta"), 2, median)
  thetaEstVar <- apply(as.matrix(rstanObj, pars = "theta"), 2, var)
  
  # estimate Factor-Corr
  corrEstMean <-  apply(as.matrix(rstanObj, pars = "PsiC[2, 1]"), 2, mean)
  corrEstMed <-   apply(as.matrix(rstanObj, pars = "PsiC[2, 1]"), 2, median)
  corrEstVar <-   apply(as.matrix(rstanObj, pars = "PsiC[2, 1]"), 2, var)

  
  # Parameters for Selection part of CrossLoadins, i.e. different configs of credible intervals
  crossQuantiles <- t(apply(crossMatrix, 2, quantile, seq(0, 1, 0.025)))
  
  # cbind and return output
  out <- cbind(1:6,
               mainEstMean, 
               mainEstMed,
               mainEstVar,
               crossEstMean,
               crossEstMed,
               crossEstVar,
               thetaEstMean,
               thetaEstMed,
               thetaEstVar,
               crossQuantiles) %>% 
        as_tibble()
  colnames(out)[1] <- "item"
  # recode output into wide format and cbind convergence into it
  out <- tidyr::pivot_wider(out, 
                            names_from = item, 
                            values_from = colnames(out[-1])) 
  
  # cbind estimates of corr (only 1 per six items) into output
  out <- cbind(out, corrEstMean, corrEstMed, corrEstVar)
  # cbind conditions into output
  out <- cbind(out, conditions)
  rownames(out) <- NULL
  
  
  # return output
  return(out)
  
}

# convergence() -----------------------------------------------------------
# takes rstan object as input and computes and returns convergence diagnostics
convergence <- function(rstanObj, conditions) {
  
  # start with Rhat and n_eff
  conv <- as.data.frame(
    t(summary(rstanObj, pars = c("lambdaMainC", 
                                 "lambdaCrossC", 
                                 "PsiC[1,2]", 
                                 "theta"))$summary[, 9:10]))
  
  # recode output into a nicer format and including condition config
  conv$parameter <- rownames(conv)
  rownames(conv) <- NULL
  
  # add max treedepth and sum divergent transitions
  tree <- subset(bayesplot::nuts_params(rstanObj), Parameter == "treedepth__")
  conv$maxTree <- max(tree$Value)
  div <- subset(bayesplot::nuts_params(rstanObj), Parameter == "divergent__")
  conv$sumDiv  <- sum(div$Value)
  # save runtime
  time = get_elapsed_time(rstanObj)
  conv$warmupT1 = time["chain:1", "warmup"]
  conv$warmupT2 = time["chain:2", "warmup"]
  conv$sampleT1 = time["chain:1", "sample"]
  conv$sampleT2 = time["chain:2", "sample"]
  # cbind conditions into output
  conv <- cbind(conv, conditions)
  # return output
  return(conv)
}

# sampling() --------------------------------------------------------------
# takes as input the conditions chain-length, warmup, n_chains, n_parallel chains &
#   all hyperparameters sourced from parameters.R
sampling <- function(pos, conditions, modelPars, nIter, samplePars){
  
  # memory allocation final output
  outputFinal <- data.frame()
  convFinal <- data.frame()

  # specify current set of condition based on pos 
  condCurrent <- conditions[pos, ]

  # compile model (if already compiled this will just not be executed)
  if (condCurrent$prior == "SVNP"){
      model <- cmdstan_model("~/1vs2StepBayesianRegSEM/stan/SVNP.stan")
  }else if (condCurrent$prior == "RHSP"){
      model <- cmdstan_model("~/1vs2StepBayesianRegSEM/stan/RHSP.stan")
  }
  
  # Draw the Samples
    for (i in 1:nIter){
      
      # simulate data based on current conditions
      dataset <- prepareDataset(conditions = condCurrent,
                                main = modelPars$main, 
                                Psi = modelPars$Psi, 
                                Theta = modelPars$Theta)
      
      # draw samples
      samples <- model$sample(data = dataset,
                              chains = samplePars$nChain, 
                              iter_warmup = samplePars$nWarmup,
                              iter_sampling = samplePars$nSampling)
      
      # save as rstan object
      rstanObj <- read_stan_csv(samples$output_files())
      
      # save Results
      output <- saveResults(rstanObj, conditions = condCurrent)
      # add iteration to output
      output$iteration <- i
      # rbind output into final output
      outputFinal <- rbind(outputFinal, output)
      
      # save convergence diagnostics
      conv <- convergence(rstanObj, conditions = condCurrent)
      # add iteration
      conv$iteration <- i
      convFinal <- rbind(convFinal, conv)
      
    }
  
  # Write output to disk (per set of conditions in an appending fashion)
  ### THIS ONLY WORKS WHEN files dont already exist, so maybe delete them before?, e.g. in main.R
  resultsName <- ifelse(condCurrent$prior == "SVNP",
                        "~/1vs2StepBayesianRegSEM/output/resultsSVNP.csv",
                        "~/1vs2StepBayesianRegSEM/output/resultsRHSP.csv")
                        
  convName <- ifelse(condCurrent$prior == "SVNP",
                     "~/1vs2StepBayesianRegSEM/output/convSVNP.csv",
                     "~/1vs2StepBayesianRegSEM/output/convRHSP.csv")
  write.table(outputFinal, 
              file = resultsName,
              append = TRUE,
              row.names = FALSE,
              col.names=!file.exists(resultsName))
  write.table(convFinal,
              file = convName,
              append = TRUE,
              row.names = FALSE,
              col.names=!file.exists(convName))
  
  # return list with results and convergence diags
  return(list(results = outputFinal,
              convergence = convFinal))
  
}

# runSim() ----------------------------------------------------------------
# function 
#runSim() <- function()

################################################################################################
# Part 2: Postprocessing Output of Simulation
################################################################################################
# selectConv -------------------------------------------------------------
# function takes whole output and trims dataset such that only converged iterations are included
#selectConv <- function(results){}
### Opsplitsenin Strict en niet zo strict? 

# computeOutcomes ---------------------------------------------------------
# Takes as input the results of a study (minus non converged) and computes all main outcomes
computeOutcomes <- function(resultsTrimmed, modelPars){
  
 ## helper function computeBias()
 computeBias <- function(est, true){
   
   bias <- abs(est-true)
   return(bias)
   
 }
 
  ### subset estimates for more concise code below:
  # subset mean estimates of parameters
  mainMean <- select(resultsTrimmed, mainEstMean_1:mainEstMean_6) 
  crossMean <- select(resultsTrimmed, crossEstMean_1:crossEstMean_6)
  thetaMean <-  select(resultsTrimmed, thetaEstMean_1:thetaEstMean_6)
  corrMean <- select(resultsTrimmed, corrEstMean)
  # subset median estimates of parameters
  mainMed <- select(resultsTrimmed, mainEstMed_1:mainEstMed_6)
  crossMed <- select(resultsTrimmed, crossEstMed_1:crossEstMed_6)
  thetaMed <- select(resultsTrimmed, thetaEstMed_1:thetaEstMed_6)
  corrMed <- select(resultsTrimmed, corrEstMed)
  # subset variances of parameters
  mainVar <- select(resultsTrimmed, mainEstVar_1:mainEstVar_6)
  crossVar <- select(resultsTrimmed, crossEstVar_1:crossEstVar_6)
  #thetaVar <- select(resultsTrimmed, thetaEstVar_1:thetaEstVar_6) ## TBA again later!!!!
  corrVar <- select(resultsTrimmed, corrEstVar)

  ### Bias Mean estimates
  # compute the biases main loadings
  biasMainMean <- t(apply(mainMean, 1, computeBias, true = modelPars$main))
  for(i in 1:6){
    colnames(biasMainMean)[i] <- paste0("biasMainMean", "_", as.character(i))
  }
  # compute biases mean estimates cross-loadings
  crossTrue <- as.data.frame(matrix(NA, nrow = nrow(resultsTrimmed), ncol = 6))
  biasCrossMean <- as.data.frame(matrix(NA, nrow = nrow(resultsTrimmed), ncol = 6))
  for(i in 1:nrow(resultsTrimmed)){
    if(resultsTrimmed[i, ]$cross == 0.5){
      crossTrue[i, ] <- rep(0.5, 6)
    }else{
      crossTrue[i, ] <- rep(0.2, 6)
    }
    biasCrossMean[i, ] <- abs(crossMean[i, ] - crossTrue[i, ])
  }  
  for(i in 1:6){
    colnames(biasCrossMean)[i] <- paste0("biasCrossMean", "_", as.character(i))
  }
  # compute biases mean estimates theta
  biasThetaMean <- t(apply(thetaMean, 1, computeBias, true = diag(modelPars$Theta)))
  for(i in 1:6){
    colnames(biasThetaMean)[i] <- paste0("biasThetaMean", "_", as.character(i))
  }
  # compute biases mean estimate factor Correlation
  biasFactCorrMean <-  abs(corrMean - modelPars$Psi[1, 2])
  names(biasFactCorrMean) <- "biasFactCorrMean"

  ## Bias Median estimates
  # compute the biases main loadings
  biasMainMed <- t(apply(mainMed, 1, computeBias, true = modelPars$main))
  # fix colnames of biases mean estimates main loadings
  for(i in 1:6){
    colnames(biasMainMed)[i] <- paste0("biasMainMed", "_", as.character(i))
  }
  # compute biases median estimates cross-loadings
  biasCrossMed <- as.data.frame(matrix(NA, nrow = nrow(resultsTrimmed), ncol = 6))
  for(i in 1:nrow(resultsTrimmed)){
    biasCrossMed[i, ] <- abs(crossMed[i, ] - crossTrue[i, ])
  }  
  for(i in 1:6){
    colnames(biasCrossMed)[i] <- paste0("biasCrossMed", "_", as.character(i))
  }
  # compute biases median estimates theta
  biasThetaMed <- t(apply(thetaMed, 1, computeBias, true = diag(modelPars$Theta)))
  for(i in 1:6){
    colnames(biasThetaMed)[i] <- paste0("biasThetaMed", "_", as.character(i))
  }
  # compute biases median estimate factor Correlation
  biasFactCorrMed <-  abs(corrMed - modelPars$Psi[1, 2])
  names(biasFactCorrMed) <- "biasFactCorrMed"
  
  ### MSE Mean estimates
  mseMainMean <- biasMainMean + mainVar
  for(i in 1:6){
    colnames(mseMainMean)[i] <- paste0("mseMainMean", "_", as.character(i))
  }
  mseCrossMean <- biasCrossMean + crossVar
  for(i in 1:6){
    colnames(mseCrossMean)[i] <- paste0("mseCrossMean", "_", as.character(i))
  }
  mseFactCorrMean <- biasFactCorrMean + corrVar
  colnames(mseFactCorrMean) <- "mseFactCorrMean"
  #mseThetaMean <- biasThetaMean + thetaVar ### TBA later
  
  ### MSE Median estimates
  mseMainMed <- biasMainMed + mainVar
  for(i in 1:6){
    colnames(mseMainMed)[i] <- paste0("mseMainMed", "_", as.character(i))
  }
  mseCrossMed <- biasCrossMed + crossVar
  for(i in 1:6){
    colnames(mseCrossMed)[i] <- paste0("mseCrossMed", "_", as.character(i))
  }
  mseFactCorrMed <- biasFactCorrMed + corrVar
  colnames(mseFactCorrMed) <- "mseFactCorrMed"
  #mseThetaMed <- biasThetaMean + thetaVar ### TBA later
    

  ## isZero, based on different selection criteria
  ##### TBA: The final outcomes, i.e. Power (p 6 Zhang et al., 2021), type-I error rates, 
  ##### Ratio correct identification to total number identified pars, no established metric
  ## Treshold: 0 if estimate smaller than 0.10
  #isZeroTres10 <- sapply(crossEst, function(x)ifelse(x < 0.10, 0, 1))
  ## Treshold: 0 if estimate is smaller than
  #
  ## 95% Credibility interval containing zero
  #isZeroCred95 <- apply(resultsTrimmed$credInterval, 
  #                      1, 
  #                      function(x) dplyr::between(0, x[1], x[2]))
  ##cred95Power <- 
  ##cred95TypeI <- 
  #
  ## HPD interval containing zero
  
  # cbind everything into a single dataframe
  out <- cbind(biasMainMean, 
               biasMainMed,
               mseMainMean,
               mseMainMed,
               biasCrossMean,
               biasCrossMed,
               mseCrossMean,
               mseCrossMed,
               biasThetaMean,
               biasThetaMed,
               # mseThetaMean, ### TBA
               # mseThetaMed,
               biasFactCorrMean,
               biasFactCorrMed,
               mseFactCorrMean,
               mseFactCorrMed
               )
  
  
  # cbind conditions into output
  out <- cbind(out, select(resultsTrimmed, prior:iteration))
  # return output
  return(out)
  
}

# Plots -----------------------------------------------------------------
# makes all required plots (generally? for AN outcome?) and saves them
#### start with bias for now,ff litertuur weer induiken & Sara spreken over wat handig is
#
#plotsMeanBias <- function(output, parameterName, condition){
#
#      name <- paste0("mean", parameterName)
#  
#      out %>% 
#        group_by(condition) %>% 
#        summarise(name = mean(parameterName)) %>% 
#        ggplot(mapping = aes(x = condition, y = name))+
#        geom_point()
#
#}       
#        
#out <- read.csv("~/1vs2StepBayesianRegSEM/output/ResultsMiniSimSVNP.csv", row.names = 1)
#        
#plotsMeanBias(out, "biasFactCorr", condition= sigma)


