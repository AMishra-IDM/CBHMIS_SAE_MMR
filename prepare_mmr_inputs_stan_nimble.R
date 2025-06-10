
##############################################
##   Prepare Inputs for STAN & NIMBLE       ##
##     Using Mock MMR Kaduna Dataset        ##
##############################################

# Load dataset
mock_data <- read.csv("mock_mmr_kaduna_dataset.csv")

# ----- COMMON -----

# Basic dimensions
N <- nrow(mock_data)
K_x <- 2  # education, distance
K_w <- 2  # conflict, community

# Covariates
x <- as.matrix(mock_data[, c("education", "distance")])
w <- as.matrix(mock_data[, c("conflict", "community")])

# Offset
live_births <- mock_data$live_births

# Observed counts
z <- mock_data$z_obs


# ----- FOR STAN -----

# Example placeholder for adjacency structure
# You will replace these with real adjacency data
# neighbors = 2 x N_edges array of adjacent region indices
neighbors <- matrix(c(1,2, 2,3, 3,4, 4,5, 5,1), nrow = 2)  # Dummy edges
N_edges <- ncol(neighbors)

# Placeholder for tau (can be 1.0 if unsure)
tau <- 1.0

stan_data <- list(
  N = N,
  z = z,
  live_births = live_births,
  K_x = K_x,
  x = x,
  K_w = K_w,
  w = w,
  N_edges = N_edges,
  neighbors = neighbors,
  tau = tau
)


# ----- FOR NIMBLE -----

nimble_constants <- list(
  n = N,
  R = N,  # assuming 1 region per row
  adj = neighbors[1,],  # needs proper adjacency formatting
  num = rep(1, N),       # dummy: each has 1 neighbor (replace!)
  weights = rep(1, length(neighbors[1,])),
  l_adj = length(neighbors[1,])
)

nimble_data <- list(
  z = z,
  education = mock_data$education,
  distance = mock_data$distance,
  conflict = mock_data$conflict,
  community = mock_data$community,
  live_births = live_births
)

nimble_inits <- list(
  a = c(-8, 0, 0),
  b = c(2, 0, 0),
  sigma = 0.5,
  epsilon = 0.5,
  nu = 0.5,
  phi = rnorm(N, 0, 1),
  theta = rnorm(N, 0, 0.5),
  gamma = rnorm(N, 0, 0.5)
)

# Save for STAN or NIMBLE usage
save(stan_data, file = "stan_data.RData")
save(nimble_constants, nimble_data, nimble_inits, file = "nimble_inputs.RData")
