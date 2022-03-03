################################################################################
# conditions.R                                            (c) J.M.B. Koch 2022
################################################################################
# This file contains the specification of all relevant study parameters
# model--------------------------------------------------------------------
# Lambda
main <- c(.5, .75, .25, .1, .5, .25)
cross2 <- c(.2, 0, 0, 0, 0, .2)
cross5 <- c(.5, 0, 0, 0, 0, .5)

# Psi
Psi <- matrix(rep(NA, 4), ncol = 2)
diag(Psi) <- 1
Psi[1, 2] <- Psi[2, 1] <- 0.5
# Theta
Theta <- diag(rep(0.3, 6))

# Hyper-Parameters: -------------------------------------------------------
# Small Variance Normal Prior ---------------------------------------------
sigma <- c(sqrt(0.1), # specified such that sigma^2 > sigma
           sqrt(0.01), 
           sqrt(0.001))

# Regularized Horseshoe Prior ---------------------------------------------
scaleGlobal <- c(0.1, 1) # scale for half-t prior omega
scaleLocal <- c(0.1, 1) # scale for half-t prior tau_j
dfGlobal <- c(1, 3) # df for half-t prior omega
dfLocal <- c(1, 3) # df for half-t prior tau_j
nu <- c(1, 3) # df IG for c^2 (slab)
scaleSlab <- c(0.1, 1, 5) # scale of slab

# Population conditions ----------------------------------------------------
#N <- c(100, 200, 300)
N <- c(100, 200) # nu ff alleen maar dit for MiniMiniSim
cross <- c(0.2, 0.5)

# specify combinations of population conditions & hyperpars per prior
condSVNP <- 
  expand.grid(
    prior = "SVNP",
    sigma = sigma,
    N = N,
    cross = cross
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
    N = N,
    cross = cross
  )

# Sampling parameters -----------------------------------------------------
nChain <- 2
nWarmup <- 50
nSampling <- 50


# parallelization parameters ----------------------------------------------
nWorkers <- nrow(condSVNP) # 16 total virtual coress

# other study parameters --------------------------------------------------
nIter <- 5 # ff dit proberen






