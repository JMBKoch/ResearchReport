################################################################################
# parameters.R
################################################################################
# This file contains the specification of all relevant study parameters
# model--------------------------------------------------------------------
# Lambda
main <- c(.5, .75, .25, .1, .5, .25)
cross <- c(.8, 0, 0, 0,  0, .8)
L <- matrix(NA, nrow = 6, ncol = 2)
L[1:3,1] <- main[1:3]
L[4:6,2] <- main[4:6]
L[4:6,1] <- cross[1:3]
L[1:3,2] <- cross[4:6]
# Psi
Psi <- matrix(rep(NA, 4), ncol = 2)
diag(Psi) <- 1
Psi[1, 2] <- Psi[2, 1] <- 0.5
# Theta
Theta <- diag(rep(0.3, 6))

# Hyper-Parameters: -------------------------------------------------------
# Small Variance Normal Prior ---------------------------------------------
sigma <- c(0.1^2, 0.01^2, 0.001^2)

# Regularized Horseshoe Prior ---------------------------------------------
scaleGlobal <- c(0.1, 0.5, 1) # scale for half-t prior omega
scaleLocal <- c(0.1, 0.5, 1) # scale for half-t prior tau_j
dfGlobal <- c(1, 2, 3) # df for half-t prior omega
dfLocal <- c(1, 2, 3) # df for half-t prior tau_j
omegaSquZero <- c(0.1, 0.5, 1) # omega^2_0 
nu <- c(1, 2, 3) # df IG for c^2 (slab)
scaleSlab <- c(0.1, 0.5, 1) # scale of slab

# Population conditions ----------------------------------------------------
n <- c(100, 200, 300)

# specify combinations of population conditions & hyperpars per prior
condRHSP <- 
      expand.grid(
         prior = prior,
         sigma = sigma, 
         scaleGlobal = scaleGlobal, 
         scaleLocal = scaleLocal,
         dfGlobal = dfGlobal,
         dfLocal = dfLocal,
         omegaSquZero = omegaSquZero,
         nu = nu,
         scaleSlab = scaleSlab, 
         n
         )
condSVNP <- 
  expand.grid(
    simga = sigma,
    n = n
  )



