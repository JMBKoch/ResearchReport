################################################################################
# functions.R                                             (c) J.M.B. Koch 2022
################################################################################
# All functions used in main.R
# dependencies: tidyverse (magrittr, tidyr, dplyr, ggplot2), mvtnorm, bayesplot
################################################################################
# Part 1: functions for executing simulation study
################################################################################
# prepareDatasets() -------------------------------------------------------
# function that prepares a list with (nIter X nrow(cond)) datasets 
prepareDatasets <- function(cond, nIter, main, cross, Psi, Theta){

  # (inner) helper functions ------------------------------------------------
  # simY() 
  # function to simulate data under desired model sourced from parameters.R
  
  ### TBA: work with different levels of cross-loadings!!!
  simY <- function(main, cross, Psi, Theta, N){
    L <- matrix(NA, nrow = 6, ncol = 2)
    L[1:3,1] <- main[1:3]
    L[4:6,2] <- main[4:6]
    L[4:6,1] <- cross[1:3]
    L[1:3,2] <- cross[4:6]
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
        sigma = conditions$sigma
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
        scaleSlab = conditions$scaleSlab # scale of slab
      )
    } 
    return(out)
  }
  
    # allocate memory for the final output, a nested list
    datasets <- list()

    # prepare 50 x "# unique combination of conditions" datasets
    for (i in 1:nrow(cond)){
      # simulate & prepare data 50x per set of conditions (row of conditionsRHSP)
      dat <- list()
      for (j in 1:nIter){
        
        if(cond[i, ]$cross == 0.2){
          cross <- 0.2
        } else{
          cross <- 0.5
        }
        
        N <- cond[i, ]$N
        Y <- simY(main, cross, L, Psi, Theta, N)
        conditions <- cond[i, ]
        dat[[j]] <-  prepareDat(Y, conditions) 
      }
      # save data in appropriate element of final output
      datasets[[i]] <- dat
    }
    # return datasets
    return(datasets)
}

# saveResults() ------------------------------------------------------------
# function takes rstan object and saves the estimes 
saveResults <- function(rstanObj, 
                        mainTrue = main, 
                        crossTrue = cross, 
                        PsiTrue = Psi, 
                        ThetaTrue = Theta, 
                        conditions){

  # estimates Lambda
  mainEstMean <- summary(rstanObj,  pars = "lambdaMainC")$summary[, 1]
  mainEstMed <- 
  mainEstSD <- 
  crossEstMean <-  summary(rstanObj, pars = "lambdaCrossC")$summary[, 1]
  crossEstMed <-  summary(rstanObj, pars = "lambdaCrossC")$summary[, 1]
  crossEstSD <- 
  
  # estimate Factor-Corr
  corrEstMean <-  summary(rstanObj, pars = "PsiC[2, 1]")$summary[, 1]
  corrEstMed <- 
  corrEstSD <- 
  # estimates Theta
  thetaEst <- summary(rstanObj, pars = "theta")$summary[, 1]
  thetaEstMed <- 
  
  # Parameters for Selection part of CrossLoadins
  credInterval50
  credInterval90
  credInterval95 <- summary(rstanObj, par = "lambdaCrossC")$summary[, c(4, 8)]
  
  # 
  out <- cbind()
  
  
  
  
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

  # cbind conditions into output
  conv <- cbind(conv, conditions)
  # return output
  return(conv)
}

# sampling() --------------------------------------------------------------
# takes as input the conditions chain-length, warmup, n_chains, n_parallel chains &
#   all hyperparameters sourced from parameters.R
##### VGM kan pos argument gwn weg
sampling <- function(pos, conditions, datasets, nIter, nChain, nWarmup, nSampling){

  # specify condition based on pos 
  cond <- conditions[pos, ]
  
  # memory allocation final output
  outputFinal <- data.frame()
  convFinal <- data.frame()

  # loop over the simulated datasets and save output object per dataset
  # compile model (if already compiled this will just not be executed)
  if (cond$prior == "SVNP"){
      model <- cmdstan_model("~/1vs2StepBayesianRegSEM/stan/SVNP.stan")
  }else if (cond$prior == "RHSP"){
      model <- cmdstan_model("~/1vs2StepBayesianRegSEM/stan/RHSP.stan")
  }
  
  # Draw the Samples
    for (i in 1:nIter){
      
      # draw samples
      samples <- model$sample(data = datasets[[pos]][[i]],
                              chains = nChain, # 2 chains 
                              # parallel_chains = nChain, # in finale setup wss geen parallele chains meer
                              iter_warmup = nWarmup, # 4000 total iterations
                              iter_sampling = nSampling)
      
      # save as rstan object
      rstanObj <- read_stan_csv(samples$output_files())
      
      # save Results
      output <- saveResults(rstanObj, conditions = cond)
      # add iteration to output
      output$iteration <- i
      # rbind output into final output
      outputFinal <- rbind(outputFinal, output)
      
      # save convergence diagnostics
      conv <- convergence(rstanObj, conditions = cond)
      # add iteration
      conv$iteration <- i
      convFinal <- rbind(convFinal, conv)
      
      ## print progress message 
      #print(paste("**************** THIS IS ITERATION ", as.character(i)))
      # does not work within ClusterApplyLB()
      
    }
  
  ### TBA: change such that output is written to disk directly (appending per iteration)
  write.table(outputFinal, 
              file = "~/1vs2StepBayesianRegSEM/output/TestAppendResults.csv",
              append = TRUE,
              row.names = FALSE,
              col.names=!file.exists("~/1vs2StepBayesianRegSEM/output/TestAppendResults.csv"))
  write.table(convFinal,
              file = "~/1vs2StepBayesianRegSEM/output/TestAppendConv.csv",
              append = TRUE,
              row.names = FALSE,
              col.names=!file.exists("~/1vs2StepBayesianRegSEM/output/TestAppendConv.csv"))
  
  # return list with results and convergence diags
  return(list(results = outputFinal,
              convergence = convFinal))
  
}

################################################################################################
# Part 2: Postprocessing Output of Simulation
################################################################################################
# checkIfConv -------------------------------------------------------------
# function takes whole output and trims dataset such that only converged iterations are included

# ComputeOutcomes ---------------------------------------------------------
# Takes as input the results of a study (minus non converged) and computes all main outcomes
ComputeOutcomes(resultsTrimmed){
  
  ## Bias 
  biasMain <- abs(mainEst-mainTrue)
  biasCross <- abs(crossEst-crossTrue)
  biasFactCorr <-  abs(corrEst - PsiTrue[1, 2])
  biasTheta <- abs(thetaEst - diag(ThetaTrue))
  
  # MSE
  mseMain <- biasMain + summary(rstanObj, pars =  "lambdaMainC")$summary[, 3]^2
  mseCross <- biasCross + summary(rstanObj, pars =  "lambdaCrossC")$summary[, 3]^2
  mseFactCorr <- biasFactCorr + summary(rstanObj, pars = "PsiC[2, 1]")$summary[, 3]^2
  mseTheta <- biasTheta + summary(rstanObj, pars = "theta")$summary[, 3]^2
  
  # isZero, based on different selection criteria
  #### TBA: The final outcomes, i.e. Power (p 6 Zhang et al., 2021), type-I error rates, 
  #### Ratio correct identification to total number identified pars, no established metric
  # Treshold: 0 if estimate smaller than 0.10
  isZeroTres10 <- sapply(crossEst, function(x)ifelse(x < 0.10, 0, 1))
  # Treshold: 0 if estimate is smaller than
  
  # 95% Credibility interval containing zero
  isZeroCred95 <- apply(resultsTrimmed$credInterval, 
                        1, 
                        function(x) dplyr::between(0, x[1], x[2]))
  #cred95Power <- 
  #cred95TypeI <- 
  
  # HPD interval containing zero
  
  
  # save output
  out <- cbind(1:6,
               biasMain,
               biasCross,
               biasTheta,
               mseMain,
               mseCross,
               mseTheta,
               isZeroTres10) %>% 
    as_tibble()
  
  # make row and colnames proper
  rownames(out) <- NULL
  colnames(out)[1] <- "item"
  # recode output into wide format and cbind convergence into it
  out <- tidyr::pivot_wider(out, 
                            names_from = item, 
                            values_from = c(biasMain, 
                                            biasCross, 
                                            biasTheta,
                                            mseMain,
                                            mseCross,
                                            mseTheta,
                                            isZeroTres10
                            )) 
  # add factCorr columns
  out$biasFactCorr <- biasFactCorr
  out$mseFactCorr <- mseFactCorr
  
  # cbind conditions into output
  out <- cbind(out, conditions)
  # return output
  return(out)
  

  
}

# Plots -----------------------------------------------------------------
# makes all required plots (generally? for AN outcome?) and saves them
#### start with bias for now,ff litertuur weer induiken & Sara spreken over wat handig is

plotsMeanBias <- function(output, parameterName, condition){

      name <- paste0("mean", parameterName)
  
      out %>% 
        group_by(condition) %>% 
        summarise(name = mean(parameterName)) %>% 
        ggplot(mapping = aes(x = condition, y = name))+
        geom_point()

}       
        
out <- read.csv("~/1vs2StepBayesianRegSEM/output/ResultsMiniSimSVNP.csv", row.names = 1)
        
plotsMeanBias(out, "biasFactCorr", condition= sigma)


