################################################################################
# conditions.R                                            (c) J.M.B. Koch 2022
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
sigma <- c(0.1, 0.01, 0.001)

# Regularized Horseshoe Prior ---------------------------------------------
scaleGlobal <- c(0.1, 1) # scale for half-t prior omega
scaleLocal <- c(0.1, 1) # scale for half-t prior tau_j
dfGlobal <- c(1, 3) # df for half-t prior omega
dfLocal <- c(1, 3) # df for half-t prior tau_j
nu <- c(1, 3) # df IG for c^2 (slab)
scaleSlab <- c(0.1, 1, 5) # scale of slab

# Population conditions ----------------------------------------------------
#N <- c(100, 200, 300)
N <- 200 # nu ff alleen maar dit for MiniMiniSim

# specify combinations of population conditions & hyperpars per prior
condSVNP <- 
  expand.grid(
    prior = "SVNP",
    sigma = sigma,
    N = N
  )

condRHSP <- 
  expand.grid(
    prior = "RHSP",
    scaleGlobal = scaleGlobal, 
    scaleLocal = scaleLocal,
    dfGlobal = dfGlobal,
    dfLocal = dfLocal,
    nu = nu,
    scaleSlab = scaleSlab, 
    N = N
  )



# Sampling parameters -----------------------------------------------------
nChain <- 2
nWarmup <- 2000
nSampling <- 2000


# other study parameters --------------------------------------------------
nIter <- 50






