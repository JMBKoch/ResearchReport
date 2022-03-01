################################################################################
# functions.R                                             (c) J.M.B. Koch 2022
################################################################################
# All functions used in main.R
#  dependencies: tidyverse (magrittr, tidyr, dplyr, ggplot2)

# prepareDatasets() -------------------------------------------------------
# function that prepares a list with (nIter X nrow(cond)) datasets 
prepareDatasets <- function(cond, nIter, L, Psi, Theta){

  # (inner) helper functions ------------------------------------------------
  # simY() 
  # function to simulate data under desired model sourced from parameters.R
  simY <- function(L, Psi, Theta, N){
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
        N <- cond[i, ]$N
        Y <- simY(L, Psi, Theta, N)
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
# function takes rstan object and computes outcomes, and saves them
saveResults <- function(rstanObj, 
                        mainTrue = main, 
                        crossTrue = cross, 
                        PsiTrue = Psi, 
                        ThetaTrue = Theta, 
                        conditions){

  # estimates Lambda #####?? Change into output summary(rstanObj)$summary???
  mainEst <- colMeans(as.matrix(rstanObj, pars = c("lambdaMainC[1]",
                                               "lambdaMainC[2]",
                                               "lambdaMainC[3]",
                                               "lambdaMainC[4]",
                                               "lambdaMainC[5]",
                                               "lambdaMainC[6]")))
  
  crossEst <- colMeans(as.matrix(rstanObj, pars = c("lambdaCrossC[1]",
                                             "lambdaCrossC[2]",
                                             "lambdaCrossC[3]",
                                             "lambdaCrossC[4]",
                                             "lambdaCrossC[5]",
                                             "lambdaCrossC[6]")))
  # estimate Factor-Corr
  corrEst <- colMeans(as.matrix(rstanObj, pars = c("Psi[2, 1]"))) # adjust to summary
  # estimates Theta
  thetaEst <- colMeans(as.matrix(rstanObj, pars = "theta"))
  
  
  ## Bias 
  # Bias Lambda
  biasMain <- abs(mainEst-mainTrue)
  biasCross <- abs(crossEst-crossTrue)
  # Bias Factor Correlation
  biasFactCorr <-  abs(corrEst - PsiTrue[1, 2])
  # Bias Theta
  biasTheta <- abs(thetaEst - diag(ThetaTrue))
  
  # MSE
  #mseMain <- biasMain + var()
  #mseCross <- 
  #mseFactCorr <- 
  #mseTheta <- 
  # Just compute as bias plus variance
  # TBUse all selection criteria 

  # save output
  out <- cbind(1:6,
               biasMain,
               biasCross,
               biasTheta) %>% 
         as_tibble()
            
  
  # make row and colnames proper
  rownames(out) <- NULL
  colnames(out)[1] <- "item"
  
  # recode output into wide format and cbind convergence into it
  out <- tidyr::pivot_wider(out, 
                            names_from = item, 
                            values_from = c(biasMain, 
                                            biasCross, 
                                            biasTheta
                                            )) 
  # add factCorr column
  out$biasFactCorr <- biasFactCorr
  
  
  # cbind conditions into output
  out <- cbind(out, conditions)
  # return output
  return(out)

}


# convergence() -----------------------------------------------------------
# takes rstan object as input and computes and returns convergence diagnostics
convergence <- function(rstanObj, conditions) {
  
  # save convergence diagnostics
  conv <- as.data.frame(
    t(summary(rstanObj, pars = c("lambdaMainC", 
                                 "lambdaCrossC", 
                                 "PsiC[1,2]", 
                                 "theta"))$summary[, 9:10]))
  # recode output into a nicer format and including condition config
  conv$parameter <- rownames(conv)
  rownames(conv) <- NULL
  
  # cbind conditions into output
  conv <- cbind(conv, conditions)
  # return output
  return(conv)
}



# sampling() --------------------------------------------------------------
# takes as input the conditions chain-length, warmup, n_chains, n_parallel chains &
#   all hyperparameters sourced from parameters.R
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
  # start  loop 
    for (i in 1:nIter){
      # do the sampling
      samples <- model$sample(data = datasets[[pos]][[i]],
                              chains = nChain, # 2 chains 
                              # parallel_chains = nChain, # in finale setup wss geen parallele chains meer
                              iter_warmup = nWarmup, # 4000 total iterations
                              iter_sampling = nSampling)
      
      # save as rstan object
      rstanObj <- read_stan_csv(samples$output_files())
      
      # save Results
      output <- saveResults(rstanObj, conditions = cond)
      #output <- cbind(output, )
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
      # does not work within ClusterApplyLB
      
    }
  
  ### TBA: change such that output is written to disk directly (appending per iteration)
  #### Does this work with clusteraplY
  # return list with results and convergence diags
  return(list(results = outputFinal,
              convergence = convFinal))
  
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


