################################################################################
# functions.R                                             (c) J.M.B. Koch 2022
################################################################################
# All functions used in main.R
# Dependencies:: mvtnorm ; tidyr; magrittr ; conditions.R, cmdstanr, rstan

# prepareDatasets() -------------------------------------------------------

# helper functions --------------------------------------------------------
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
      scaleSlab = conditions$scaleSlab
    )
  } 
  return(out)
}

# function that prepares a list with (nIter X nrow(cond)) datasets 
#  has two inner functions, one to simulate Y and one to prepare a single 
#  dataset (list) ready to be processed by stan 
prepareDatasets <- function(conditions, nIter = nIter, L = L, Psi = Psi, Theta = Theta){
  
    # allocate memory for the final output, a nested list
    datasets <- list()

    # prepare 50 x "# unique combination of conditions" datasets
    for (i in 1:nrow(conditions)){
      # simulate & prepare data 50x per set of conditions (row of conditionsRHSP)
      dat <- list()
      for (j in 1:nIter){
        Y <- simY(L, Psi, Theta, N = conditions[i, ]$N)
        dat[[j]] <-  prepareDat(Y, conditions = conditions[i, ]) 
      }
      # save data in appropriate element of final output
      datasets[[i]] <- dat
    }
    # return datasets
    return(datasets)
    
}

# saveOutput() ------------------------------------------------------------
# function takes rstan object and computes outcomes, and saves them
saveOutput <- function(rstanObj, 
                       mainTrue = main, 
                       crossTrue = cross, 
                       PsiTrue = Psi, 
                       ThetaTrue = Theta, 
                       conditions = conditions){

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
  corrEst <- colMeans(as.matrix(rstanObj, pars = c("Psi[2, 1]")))
  # estimates Theta
  thetaEst <- colMeans(as.matrix(rstanObj, pars = "theta"))
  
  # Bias Lambda
  biasMain <- abs(mainEst-mainTrue)
  biasCross <- abs(crossEst-crossTrue)
  # Bias Factor Correlation
  biasFactCorr <-  abs(corrEst - PsiTrue[1, 2])
  # Bias Theta
  biasTheta <- abs(thetaEst - diag(ThetaTrue))
  
  # TBA: MSE
  # TBA: True & False positives in estimating truly non-0 as non-0
  #   THINK WELL OF SELECTION CRITERIA

  # save output
  out <- as.data.frame(
          cbind(
               1:6,
               biasMain,
               biasCross,
               biasFactCorr,
               biasTheta
               )
                        )
  
  # make row and colnames proper
  rownames(out) <- NULL
  colnames(out)[1] <- "item"
  
  # recode output into wide format and cbind convergence into it
  out <- tidyr::pivot_wider(out, 
                            names_from = item, 
                            values_from = c(biasMain, 
                                            biasCross, 
                                            biasFactCorr, 
                                            biasTheta
                                            ))
  
  
  # cbind conditions into output
  out <- cbind(out, conditions)
  # return output
  return(out)

}

# sampling() --------------------------------------------------------------
# takes as input the conditions chain-length, warmup, n_chains, n_parallel chains &
#   all hyperparameters sourced from parameters.R
sampling <- function(datasets, cond, nChain = nChain, nWarmup = nWarmup, nSampling = nSampling){

# loop over the simulated datasets and save output object per dataset
# compile model (if already compiled this will just not be executed)
if (cond$prior == "SVNP"){
    model <- cmdstan_model("~/1vs2StepBayesianRegSEM/stan/SmallVarNormal.stan")
}else if (cond$prior == "RHSP"){
    model <- cmdstan_model("~/1vs2StepBayesianRegSEM/stan/SmallVarNormal.stan")
}

# Draw the Samples
# allocate memory for final output
outputFinal <- data.frame()

# start nested loop (i = conditions config, j = iteration)
for (i in 1:nrow(cond)){
  for (j in 1:nIter){
    # do the sampling
    samples <- model$sample(data = datasets[[i]][[j]],
                            chains = nChain,
                            parallel_chains = nChain,
                            iter_warmup = nWarmup, # 4000 total iterations
                            iter_sampling = nSampling)
    
    # save as rstan object
    rstanObj <- read_stan_csv(samples$output_files())
    
    # save desired ouput
    output <- saveOutput(rstanObj, conditions = cond[i, ])
    #output <- cbind(output, )
    # add iteration to output
    output$iteration <- j
    # rbind output into final output
    outputFinal <- rbind(outputFinal, output)
    
    ## print progress message 
    print(paste("**************** THIS IS ITERATION ", 
                as.character(j), 
                "OF ROW ", 
                as.character(i),
                "****************"))
    
  }
 }
}
# Plots -----------------------------------------------------------------
# makes all required plots (generally? for AN outcome?) and saves them
# plotsBias <- ()

# convergence() -----------------------------------------------------------
# takes rstan object as input and computes and returns convergence diagnostics
convergence <- function(rstanObj) {
  
  as.data.frame(
    summary(rstanObj, pars = c("lambdaMainC", 
                               "lambdaCrossC", 
                               "PsiC[1,2]", 
                               "theta"))$summary[, 9:10]
              )
}




