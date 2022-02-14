################################################################################
# functions.R
################################################################################
# All functions used in main.R

# simDat() ----------------------------------------------------------------
# function to simulate data under desired model sourced from parameters.R
# function
simDat <- function(L, Psi, Theta, N){
      Sigma <- L%*%Psi%*%t(L) + Theta
      Y <- mvtnorm::rmvnorm(N, rep(0, 6), Sigma)
      return(Y)
}

# prepareDat() ------------------------------------------------------------
# function to prepare stan data object from simdat() based on a unique
#   combination of hyper-pars
# This functions helps to clearly separate population 
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

# sampling() --------------------------------------------------------------
# takes as input the chain-length, warmup, n_chains, n_parallel chains &
#   all hyperparameters sourced from parameters.R
sampling <- function()

# output() ----------------------------------------------------------------
# function takes rstan object and computes outcomes, and saves them

output <- function(rstanObj, 
                   mainTrue = main, 
                   crossTrue = cross, 
                   PsiTrue = Psi, 
                   ThetaTrue = Theta
                   #, 
                   #cond = conditions[i, ]
                   ){
  

  # estimates Lambda
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
  # estimates Factor-Corr
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
  
  # TBA: save output (in list?)
  
  
  # Output
  
  out <- as.data.frame(cbind(
    
    #cond, # save results
               biasMain,
               biasCross,
               biasFactCorr,
               biasTheta

  ))
  
  rownames(out) <- NULL
  
  # recode output into wide format
  
  # cbind conditions into output
 # outFinal <-cbind(
 #   #cond, 
 #             biasMain,
 #             biasCross,
 #             biasFactCorr,
 #             biasTheta
 #             )
  
  # return output
  #return(outFinal)
  return(out)
  
}

# plots() -----------------------------------------------------------------
# makes all required plots and saves them

# convergence() -----------------------------------------------------------
# takes rstan object as input and computes and returns convergence diagnostics

  #



