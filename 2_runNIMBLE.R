#################################################################
##   LGA estimates of community MMR                            ##
##   2. RUN_NIMBLE                                             ##
##   Purpose: Run model using NIMBLE                           ##
##   Author: Anu Mishra                                        ##  
##   Created 6/9/25                                            ##
#################################################################

library(nimble)

rm(list=ls())

start_time <- Sys.time()

set.seed(105)
setwd("C:/Users/anumi/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub/CBHMIS_SAE_MMR")

model_run_date <- Sys.Date()
model_output <- paste0("model_out_",model_run_date,".rds")  ## which version of the model you want to read in
model_out_dir <- paste0("model runs/model_",model_run_date,"/")
dir.create(model_out_dir)

#### Load NIMBLE data
load("nimble_inputs_final.RData")
modData <- read.csv("cbhmis_data_for_model.csv")


# --- Model code block ---
mmr_code <- nimbleCode({
  for (i in 1:n) {
    # Under-reporting model (logit scale)                   # pi is probability of death getting reported -- one for each region
    pi[i] <- ilogit(b[1] + ancRef[i] * b[2] + gamma[i])     # b[1] baseline log reporting rate
                                                            # b[2] is effect of ANC references on reporting
                                                            # gamma is LGA-specific unexplained noise for reporting
    
    # Incidence model (log scale)
    lambda[i] <- exp(log(live_births[i]) +                  # lambda is true maternal deaths count -- one for each region
                       a[1] + educ[i] * a[2] +              #a[1] baseline log MMR, 
                       travel[i] * a[3] + anc4[i] * a[4] +  #a[2]-a[4] effect of XX predictor on MMR
                       phi[i] + theta[i])                   #phi is the spatial random effect  
                                                            #theta is LGA-specific unexplained noise for incidence
    # Likelihood for observed deaths
    z[i] ~ dpois(pi[i] * lambda[i])                         #z is the under-reported observed maternal deaths  -- one for each region                 
                                                            
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
    a[k] ~ dnorm(0, sd = 2.5)  # vague priors for all -- this is less than Stoner et all, but I think that was going to give implausible values based on prior predictive checks        
  }
  
  # Priors: reporting coefficients
  b[1] ~ dnorm(2, sd = 0.6)   # informative prior on reporting rate --- NEEDS FOLLOW-UP UPDATE TO BE BASED ON OBSERVED REPORTING RATE
  b[2] ~ dnorm(0, sd = 2.5)
  
  # Priors: variance parameters   --- these are truncated normal distributions truncated at 0
  sigma ~ T(dnorm(0, 1), 0, )                              #sigma SD for theta (unstructred RE for incidence)
  epsilon ~ T(dnorm(0, 1), 0, )                            #epsilon SD for gamma (unstructred RE for reporting)
  #nu ~ T(dnorm(0, 1), 0, )
  #nu ~ T(dnorm(0, 2), 0, )                                #Update 7/2, trying a slightly more prior
  tau <- 1 / (nu^2)                                        #tau is precision (1/var) for spatial RE
})



# --- Prior predictive checks --
## temporarily resetting z so I can do prior predictive check
nimble_data_prior <- nimble_data
nimble_inits_prior <- nimble_inits
nimble_data_prior$z <- NULL
nimble_inits_prior$z <- rep(NA, nrow(modData))

mmr_model <- nimbleModel(mmr_code, constants = nimble_constants, data = nimble_data_prior, inits = nimble_inits_prior)
compiled_model <- compileNimble(mmr_model)

sim_results <- list()
n_sim <- 100

for (i in 1:n_sim) {
  simulate(compiled_model, includeData = FALSE)
  lambda <- compiled_model$lambda
  y_sim <- rpois(length(lambda), lambda)  
  
  sim_results[[i]] <- list(
    z = compiled_model$z,
    pi = compiled_model$pi, 
    y = y_sim,
    lambda = compiled_model$lambda,  
    MMR = y_sim / mmr_model$live_births * 100000,
    total_z = sum(compiled_model$z),
    total_y = sum(y_sim)
  )
}


### checking priors compared to real data
total_z <- sapply(sim_results, `[[`, "total_z")
total_y <- sapply(sim_results, `[[`, "total_y")
hist(total_z, main = "Prior Predictive: Total Observed Deaths (z)", xlab = "Sum(z_i)")
abline(v=sum(nimble_data$z),col="red")
hist(total_y, main = "Prior Predictive: Total True Deaths (y)", xlab = "Sum(y_i)")

mmr_all <- do.call(rbind, lapply(sim_results, `[[`, "MMR"))  # rows = simulations, cols = LGAs
boxplot(mmr_all, main = "Prior Predictive MMR by LGA", ylab = "MMR",ylim=c(0,1e6))
points(modData$MMR, col = "red", pch = 16)  # observed overlay

pi_all <- do.call(rbind, lapply(sim_results, `[[`, "pi"))
hist(pi_all, breaks = 30, main = "Prior Predictive: Reporting Probability (π)", xlab = "π_i")

lambda_all <- do.call(rbind, lapply(sim_results, `[[`, "lambda"))
hist(lambda_all, breaks = 30, main = "Prior Predictive: Poisson λ_i", xlab = "λ_i")
### NEEDS FOLLOW-UP: Currently the prior values are implausible for data generating. We could tighten everything

# --- Build and compile model ---
mmr_model <- nimbleModel(mmr_code, constants = nimble_constants, data = nimble_data, inits = nimble_inits)
compiled_model <- compileNimble(mmr_model)


# --- Configure MCMC ---
### DECIDE IF WANT WAIC here
mmr_conf <- configureMCMC(mmr_model,
                          monitors = c("a", "b", "sigma", "epsilon", "tau","nu","phi", "theta", "pi", "lambda"),
                          useConjugacy = TRUE,enableWAIC = TRUE)

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
                   niter = 100000,
                   nburnin = 50000,
                   thin = 1,
                   summary = FALSE,
                   samplesAsCodaMCMC = TRUE,
                   WAIC = TRUE)

end_time <- Sys.time()

elapsed_time <- end_time - start_time
print(elapsed_time)

saveRDS(samples, file = paste0(model_out_dir,"/",model_output))

#### NEEDED FOR SAVING TO CLUSTER ###
# https://adb-8293166442374260.0.azuredatabricks.net/explore/data/volumes/idm_general/cbhmis_sae_mmr/cbhmis_sae_mmr?o=8293166442374260
# file_name = "samples_test_output.rds"
# rds_save_path = paste0("/Volumes/idm_general/cbhmis_sae_mmr/cbhmis_sae_mmr/", file_name)
# dir.create(dirname(rds_save_path), showWarnings = F, recursive = T)
# 
# saveRDS(samples, file = rds_save_path)
