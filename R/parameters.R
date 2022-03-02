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
hyperSVNP <- list()
hyperSVNP$sigma <- sigma

# Regularized Horseshoe Prior ---------------------------------------------
hyperRHSP <- list(
  dfGlobal = c(1), # df for half-t prior omega
  dfLocal = c(1), # df for half-t prior tau_j
  omegaSquZero = c(1), # omega^2_0 
  nu = c(1), # df IG for c^2
  s2 = c(1)
)
# or as data frame?

# Other conditions --------------------------------------------------------
# Sample Sizes
n <- c(100, 200)
# 
prior <- c("S", "SV")

expand.grid(n, prior)