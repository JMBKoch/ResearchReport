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
        scaleGlobal = conditions$scaleGlobal,
        scaleLocal = conditions$scaleLocal,
        dfGlobal = conditions$dfGlobal, # df for half-t prior omega
        dfLocal = conditions$dfLocal, # df for half-t prior tau_j
        omegaSquZero = conditions$omegaSquZero, # omega^2_0 
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

output <- function(rstanObj, L, Psi, Theta){
  
  # estimates Lambda
  main <- colMeans(as.matrix(FitA, pars = c("lambdaMainC[1]",
                                    "lambdaMainC[2]",
                                    "lambdaMainC[3]",
                                    "lambdaMainC[4]",
                                    "lambdaMainC[5]",
                                    "lambdaMainC[6]")))
  
  cross <- colMeans(as.matrix(FitA, pars = c("lambdaCrossC[1]",
                                             "lambdaCrossC[2]",
                                             "lambdaCrossC[3]",
                                             "lambdaCrossC[4]",
                                             "lambdaCrossC[5]",
                                             "lambdaCrossC[6]")))
  # estimates Factor-Corr
  corr <- colMeans(as.matrix(FitA, pars = c("Psi[2, 1]")))
  # estimates Theta
  theta <- colMeans(as.matrix(FitA, pars = "theta"))
  
  # Bias Lambda
  biasMain <- abs(main-c(L1, L2))
  biasCross <- abs(cross-c(CL1, CL2))
  # Bias Factor Correlation
  biasFactCorr <-  abs(corr-0.5)
  # Bias Thetay
  biasTheta <- abs(theta - rep(0.3, 6))
  
  # TBA: MSE
  
  # TBA: True & False positives in estimating truly non-0 as non-0
  
  # TBA: save output (in list?)
  
  # TBA: save output (format?)
}


# plots() -----------------------------------------------------------------
# makes all required plots and saves them

# convergence() -----------------------------------------------------------
# takes rstan object as input and computes and returns convergence diagnostics

  #

# main() ------------------------------------------------------------------
# runs the main simulation, including the setup of parralell computing etc.

