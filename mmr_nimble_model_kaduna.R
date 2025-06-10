
##############################################
##     MMR Model with Under-Reporting       ##
##        Kaduna State, Nigeria (2022–24)   ##
##############################################

library(nimble)

# Constants
n <- 23  # number of LGAs
R <- 23  # number of spatial regions (same here)

# Adjacency structure (pre-processed)
adj <- c(2, 3, 4, 1, 4, 8, 9, 10, 14, 1, 4, 16, 19, 23, 1, 2, 3, 9, 10, 12, 14, 21, 23,
         15, 18, 21, 7, 8, 11, 22, 6, 13, 20, 22, 2, 6, 11, 14, 22, 2, 4, 10, 2, 4, 9,
         6, 8, 4, 13, 14, 15, 17, 21, 22, 7, 12, 22, 2, 4, 8, 12, 22, 5, 12, 17, 21,
         3, 18, 19, 21, 12, 15, 5, 16, 21, 3, 16, 21, 23, 7, 4, 5, 12, 15, 16, 18,
         19, 23, 6, 7, 8, 12, 13, 14, 3, 4, 19, 21)
num <- c(3, 6, 5, 9, 3, 4, 4, 5, 3, 3, 2, 7, 3, 5, 4, 4, 2, 3, 4, 1, 8, 6, 4)
weights <- rep(1, length(adj))
l_adj <- length(adj)

# Covariates (you must create or load these)
# These should be n×1 or n×p vectors/matrices
# education[n], distance[n], conflict[n], community[n], live_births[n]

# Model code
mmr_code <- nimbleCode({
  for (i in 1:n) {
    pi[i] <- ilogit(b[1] + conflict[i]*b[2] + community[i]*b[3] + gamma[i])
    lambda[i] <- exp(log(live_births[i]) + a[1] + education[i]*a[2] + distance[i]*a[3] + 
                     phi[i] + theta[i])
    z[i] ~ dpois(pi[i] * lambda[i])
    gamma[i] ~ dnorm(0, sd=epsilon)
  }

  for (j in 1:R) {
    theta[j] ~ dnorm(0, sd=sigma)
  }

  phi[1:R] ~ dcar_normal(adj=adj[1:l_adj], num=num[1:R], tau=tau, zero_mean=1)

  # Priors for regression parameters
  a[1] ~ dnorm(-8, sd=1)
  a[2] ~ dnorm(0, sd=10)
  a[3] ~ dnorm(0, sd=10)

  b[1] ~ dnorm(2, sd=0.6)     # informative prior for reporting rate
  b[2] ~ dnorm(0, sd=10)
  b[3] ~ dnorm(0, sd=10)

  # Priors for random effect variances
  sigma ~ T(dnorm(0, 1), 0, )
  epsilon ~ T(dnorm(0, 1), 0, )
  nu ~ T(dnorm(0, 1), 0, )
  tau <- 1 / (nu^2)
})

# Constants and data list
mmr_constants <- list(
  n = n,
  R = R,
  adj = adj,
  num = num,
  weights = weights,
  l_adj = l_adj
)

# Example placeholders for covariate data
# Replace with your actual data
mmr_data <- list(
  z = rep(1, n),                      # observed maternal deaths (placeholder)
  education = rnorm(n),              # example covariate
  distance = rnorm(n),
  conflict = rnorm(n),
  community = rnorm(n),
  live_births = rep(1000, n)         # offset (placeholder)
)

mmr_inits <- list(
  a = c(-8, 0, 0),
  b = c(2, 0, 0),
  sigma = 0.5,
  epsilon = 0.5,
  nu = 0.5,
  phi = rnorm(R, 0, 1),
  theta = rnorm(R, 0, 0.5),
  gamma = rnorm(n, 0, 0.5)
)
