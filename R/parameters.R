################################################################################
# parameters.R
################################################################################
# This file contains the specification of all relevant study parameters

# data --------------------------------------------------------------------
# Lambda
main <- c(.5, .75,.25,.1, .5, .25)
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
sigma <- c(sqrt(0.1), sqrt(0.01), sqrt(0.001))

# Regularized Horseshoe Prior ---------------------------------------------


# Conditions --------------------------------------------------------------
# Sample Sizes
n <- c(100, 200)

# 