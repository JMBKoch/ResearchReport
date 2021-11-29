################################################################################
# functions.R: All functions used in main.R, J.M.B. Koch, 2021
################################################################################
# simdat() ----------------------------------------------------------------
# function to simulate data under desired model
# function takes the four objects below as input; adjust as desired 
# TBA: Extract the matrix from Hyper-parameters.R instead in final setup
# 1. Lambda
L1 <- matrix(c(.5,.75,.25), 3, 1)
CL1 <- matrix(c(.25, 0, 0), 3, 1)
L2 <- matrix(c(.1, .5, .25), 3 , 1)
CL2 <- matrix(c(0, 0, .25), 3, 1)
L <-  cbind(rbind(L1, CL1), rbind(CL2, L2))
# 2. Psi
Psi <- matrix(rep(NA, 4), ncol = 2)
diag(Psi) <- 1
Psi[1, 2] <- Psi[2, 1] <- 0.5
# 3. Theta
Theta <- diag(rep(0.3, 6))
# 4. N
N <- 200

# function
simdat <- function(L, Psi, Theta, N){
      Sigma <- L%*%Psi%*%t(L) + Theta
      Y <- mvtnorm::rmvnorm(N, rep(0, 6), Sigma)
      return(Y)
}


# sampling() --------------------------------------------------------------
# takes as output the chain-length, warmup, n_chains, n_parallel chains &
#   all hyperparameters

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
  
  # TBA: save output in list
  
}


# convergence() -----------------------------------------------------------
# takes rstan object as input and computes and returns convergence diagnostics

# main() ------------------------------------------------------------------
# runs the main simulation, including the setup of parralell computing etc.

