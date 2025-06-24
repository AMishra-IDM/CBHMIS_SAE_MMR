#################################################################
##   LGA estimates of community MMR                            ##
##   2. RUN_NIMBLE                                             ##
##   Purpose: Run model using NIMBLE                           ##
##   Author: Anu Mishra                                        ##  
##   Created 6/9/25                                            ##
#################################################################

library(nimble)

rm(list=ls())

set.seed(105)
setwd("C:/Users/anumi/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub/CBHMIS_SAE_MMR")
outdir <- "C:/Users/anumi/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub/CBHMIS_SAE_MMR/output/"

#### Load NIMBLE data
load("nimble_inputs_final.RData")


# --- Model code block ---
mmr_code <- nimbleCode({
  for (i in 1:n) {
    # Under-reporting model (logit scale)
    pi[i] <- ilogit(b[1] + ancRef[i] * b[2] + gamma[i])
    
    # Incidence model (log scale)
    lambda[i] <- exp(log(live_births[i]) +
                       a[1] + educ[i] * a[2] +
                       travel[i] * a[3] + anc4[i] * a[4] +
                       phi[i] + theta[i])
    
    # Likelihood for observed deaths
    z[i] ~ dpois(pi[i] * lambda[i])
    
    # Observation-level noise in reporting probability
    gamma[i] ~ dnorm(0, sd = epsilon)
  }
  
  # Unstructured spatial effect (iid)
  for (j in 1:R) {
    theta[j] ~ dnorm(0, sd = sigma)
  }
  
  # Structured spatial effect (ICAR)
  phi[1:R] ~ dcar_normal(adj = adj[1:l_adj], num = num[1:R], tau = tau, zero_mean = 1)
  
  # Priors: incidence coefficients
  for (k in 1:4) {
    a[k] ~ dnorm(0, sd = 10)  # vague priors for all
  }
  
  # Priors: reporting coefficients
  b[1] ~ dnorm(2, sd = 0.6)   # informative prior on reporting rate
  b[2] ~ dnorm(0, sd = 10)
  
  # Priors: variance parameters
  sigma ~ T(dnorm(0, 1), 0, )
  epsilon ~ T(dnorm(0, 1), 0, )
  nu ~ T(dnorm(0, 1), 0, )
  tau <- 1 / (nu^2)
})



# --- Build and compile model ---
mmr_model <- nimbleModel(mmr_code, constants = nimble_constants, data = nimble_data, inits = nimble_inits)
compiled_model <- compileNimble(mmr_model)

# --- Configure MCMC ---
### DECIDE IF WANT WAIC here
mmr_conf <- configureMCMC(mmr_model,
                          monitors = c("a", "b", "sigma", "epsilon", "nu", "phi", "theta", "pi", "lambda"),
                          useConjugacy = TRUE)

# Customize samplers (e.g., slice sampling for highly correlated terms) -- Not sure if I need this, but following STONER paper
mmr_conf$removeSamplers(c("a[1]", "b[1]", "sigma", "nu", "epsilon"))
mmr_conf$addSampler(target = c("a[1]", "b[1]", "epsilon"), type = "AF_slice")
mmr_conf$addSampler(target = c("sigma", "nu"), type = "AF_slice")

# --- Compile and run MCMC ---
mmr_mcmc <- buildMCMC(mmr_conf)
compiled_mcmc <- compileNimble(mmr_mcmc, project = mmr_model)

# --- Run model (example: 4 chains × 100,000 iterations each) ---
samples <- runMCMC(compiled_mcmc,
                   nchains = 4,
                   niter = 50000,
                   nburnin = 25000,
                   thin = 1,
                   summary = TRUE,
                   WAIC = FALSE)

saveRDS(samples, file = paste0(outdir,"/NIMBLE output/samples_full_output.rds"))
