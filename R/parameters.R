################################################################################
# parameters.R                                            (c) J.M.B. Koch 2022
################################################################################
# This file contains the specification of all relevant study parameters
# Packages ----------------------------------------------------------------
# specify packages that are required for executing the simulation
packages <- c("cmdstanr", # MCMC sampling using stan
              "rstan", # postprocessing of samples
              "tidyverse", # data wrangling, plotting, pipes
              "mvtnorm", # data simulation
              "parallel", # parallelization
              "bayesplot" # convergence diagnostics 
              )
# make sure that packages are installed if not present
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
    }
  }
)

# Model--------------------------------------------------------------------
# Lambda
main <- c(.75, .75, .75, .75, .75, .75)
cross2 <- c(.2, 0, 0, 0, 0, .2) # TBA: make this more generalizable?
cross5 <- c(.5, 0, 0, 0, 0, .5)
# Psi
Psi <- matrix(rep(NA, 4), ncol = 2)
diag(Psi) <- 1
Psi[1, 2] <- Psi[2, 1] <- 0.5
# Theta
Theta <- diag(rep(0.3, 6))
# save all in one object for easier passing to functions
modelPars <- list(
                main = main,
                cross2 = cross2,
                cross5 = cross5,
                Psi = Psi,
                Theta = Theta
                  )

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
N <- c(100, 200)
#N <- 200
cross <- c(0.2, 0.5)
#cross <- 0.5

# Making condition objects ------------------------------------------------
## gwn opsplitsen in condPop en de overige twee
condPop   <- 
  expand.grid(
    N = N,
    cross = cross
  )

condSVNP <- 
  expand.grid(
    prior = "SVNP",
    sigma = sigma
  )

condRHSP <- 
  expand.grid(
    prior = "RHSP",
    scaleGlobal = scaleGlobal, 
    scaleLocal = scaleLocal,
    dfGlobal = dfGlobal,
    dfLocal = dfLocal,
    nu = nu,
    scaleSlab = scaleSlab
  )

# Sampling parameters -----------------------------------------------------
# save in one list for easier passing to functions
samplePars <- list(
                nChain = 2,
                nWarmup = 2000,
                nSampling = 4000
                )

# Parallelization Parameter -----------------------------------------------
nClusters <- 46

# other study parameters --------------------------------------------------
nIter <- 200
