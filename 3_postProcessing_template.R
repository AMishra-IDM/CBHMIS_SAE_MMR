#################################################################
##   LGA estimates of community MMR                            ##
##   4. Post-processing                                        ##
##   Purpose: Fit checks and post-processing of raw model      ##
##            model output data                                ##
##   Author: Anu Mishra                                        ##  
##   Created 6/9/25                                            ##
#################################################################
rm(list=ls())

library(tidyverse)
library(ggplot2)
library(nimble)
library(coda) # For manipulation of MCMC results.
library(sf)
library(gridExtra)
library(patchwork)
library(ggrepel)


source("~/GitHub/CBHMIS_SAE_MMR/0_functions.R", echo=TRUE)
setwd("C:/Users/anumi/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub/CBHMIS_SAE_MMR")

model_date <- "2025-07-02_no_ANCRef"
model_output <- paste0("model_out_",model_date,".rds")  ## which version of the model you want to read in
model_out_dir <- paste0("model runs/model_",model_date,"/")

#### Loading data ####
## model output
full_out <- readRDS(paste0(model_out_dir,"/",model_output))
samples <- full_out$samples

## shape file 
nga_shp <- st_read("raw data/gadm41_NGA_shp/gadm41_NGA_2.shp")
kdn_shp <- nga_shp[nga_shp$NAME_1=="Kaduna",]
kdn_shp$NAME_2 <- gsub(pattern = " ",replacement = "_", kdn_shp$NAME_2)
kdn_shp$NAME_2 <- gsub(pattern = "'",replacement = "", kdn_shp$NAME_2)

## CBHMIS data
data <- read.csv("cbhmis_data_for_model.csv")
load("nimble_inputs_final.RData")
LGA_names <- data$LGA
LGA_names2 <- data.frame("LGA"=1:23,"name"=LGA_names)




#### 1: Model Diagnostics and Convergence ####

# traceplots and convergence 
pdf(file = paste0(model_out_dir,"/trace_plots.pdf"),height=8,width=12)
traceplot(samples[,c('a[1]','a[2]','a[3]','a[4]',   
                         'b[1]','b[2]',                  
                         'epsilon','sigma','tau','nu',
                         'phi[1]','phi[4]','phi[8]','phi[10]','phi[12]','phi[17]','phi[23]')])    
dev.off()
### Adding a couple of phis because tau looks bad but nu looks ok --> phi's seem to be mixing well

gelman.diag(samples[,c('a[1]','a[2]','a[3]','a[4]',   
                       'b[1]','b[2]',                  
                       'epsilon','sigma','tau','nu',
                       'phi[1]','phi[4]','phi[8]','phi[10]','phi[12]','phi[17]','phi[23]')])

effectiveSize(samples[,c('a[1]','a[2]','a[3]','a[4]',   
                         'b[1]','b[2]',                  
                         'epsilon','sigma','tau','nu',
                         'phi[1]','phi[4]','phi[8]','phi[10]','phi[12]','phi[17]','phi[23]')])



#### 2: Prior checks ####
## prior vs posterior
param_names <- c('a[1]', 'a[2]', 'a[3]', 'a[4]', 
                 'b[1]', 'b[2]', 'epsilon', 'sigma', 'nu')

prior_defs <- list(
  "a[1]" = list(type = "norm", mean = log(0.01), sd = 0.5),
  "a[2]" = list(type = "norm", mean = 0, sd = 1),
  "a[3]" = list(type = "norm", mean = 0, sd = 1),
  "a[4]" = list(type = "norm", mean = 0, sd = 1),
  "b[1]" = list(type = "norm", mean = 2, sd = 0.6),
  "b[2]" = list(type = "norm", mean = 0, sd = 1),
  "epsilon" = list(type = "half-normal", mean = 0, sd = 1),  # tighter prior, often used in under-reporting
  "sigma" = list(type = "half-normal", mean = 0, sd = 0.3),
  "nu" = list(type = "half-normal", mean = 0, sd = 2)
)

plot_list <- lapply(param_names, function(param) {
  plot_prior_vs_posterior(samples = samples, param_name = param, prior_defs = prior_defs)
})

pdf(file = paste0(model_out_dir, "/prior_post_plots.pdf"), height = 8, width = 12)
print(marrangeGrob(grobs = plot_list, nrow = 2, ncol = 2) )
dev.off()


#### 3: Posterior predictive checks ### 
## checking to see if  my model were true, would I be likely to see data like the data I actually observed
post_mat <- as.matrix(samples)
n_draws <- nrow(post_mat) 
n_areas <- nrow(data)

sim_y <- matrix(NA, nrow = n_draws, ncol = n_areas)
sim_z <- matrix(NA, nrow = n_draws, ncol = n_areas)
sim_mmr <- matrix(NA, nrow = n_draws, ncol = n_areas)

## simulating reported (z) vs true (y) deaths and corresponding MMR 
for (d in 1:n_draws) {
  lambda <- post_mat[d, paste0("lambda[", 1:n_areas, "]")]
  pi <- post_mat[d, paste0("pi[", 1:n_areas, "]")]
  
  # Simulate y_i
  y <- rpois(n_areas, lambda = lambda)
  sim_y[d, ] <- y
  
  # Simulate z_i
  sim_z[d, ] <- rbinom(n_areas, size = y, prob = pi)
  
  # e. Simulate MMR
  sim_mmr[d, ] <- y / nimble_data$live_births * 1e5
}

pdf(file = paste0(model_out_dir, "/posterior_pred_checks.pdf"), height = 8, width = 12)
par(mfrow=c(2,3))

### Observed vs simulated deaths (not corrected deaths )
hist(rowSums(sim_z), col = "skyblue", main = "Simulated vs Observed Total Deaths",
     xlab = "Simulated total z_i")
abline(v = sum(nimble_data$z), col = "red", lwd = 2, lty = 2)

plot(nimble_data$z, colMeans(sim_z), pch = 16, xlab = "Observed z_i", ylab = "Predicted z_i",
     main = "Observed vs Expected Reported Deaths by LGA")
abline(0, 1, col = "red", lty = 2)

## coverage of observed deaths 
ci_lower <- apply(sim_z, 2, quantile, 0.025)
ci_upper <- apply(sim_z, 2, quantile, 0.975)
ci_width <- ci_upper - ci_lower

coverage <- nimble_data$z >= ci_lower & nimble_data$z <= ci_upper
mean(coverage)  # Proportion covered by 95% CrI
plot(ci_width, type = "h", main = "Width of 95% CrIs for z_i", xlab = "LGA Index", ylab = "CrI Width",xaxt = "n")
axis(1, at = 1:length(LGA_names), labels = LGA_names, las = 2, cex.axis = 0.7)
legend("topleft",c(paste0("Mean coverage =",mean(coverage))))

## summaries of reporting
y_mean <- colMeans(sim_y)
plot(nimble_data$z, y_mean,
     xlab = "Observed z_i (Reported Deaths)", ylab = "Posterior mean of y_i (True Deaths)",
     main = "Observed vs. Posterior Expected True Deaths")
abline(0, 1, col = "red", lty = 2)
sum(nimble_data$z > y_mean)

## Reporting adjustment factor
ratio <- y_mean / nimble_data$z
barplot(ratio, main = "Implied Under-reporting Adjustment: y_i / z_i",
        ylab = "Correction Factor", names.arg = LGA_names, las = 2)
abline(h = 1, col = "red", lty = 2)
dev.off()
par(mfrow=c(1,1))


#### Log mean squared error
rate_z <- nimble_data$z / nimble_data$live_births
rate_z_hat <- colMeans(sim_z) / nimble_data$live_births

lmse <- log(mean((rate_z - rate_z_hat)^2))
print(paste("LMSE =", round(lmse, 4)))

#### p-value for posterior predictive check
observed_total <- sum(nimble_data$z)
simulated_totals <- rowSums(sim_z)

p_ppc <- mean(simulated_totals > observed_total)
print(paste("Posterior Predictive p-value =", round(p_ppc, 3)))

##### 3: Summary Estimates #######
samp_matrix <- as.matrix(samples)
lambda_cols <- grep("^lambda\\[[0-9]+\\]$", colnames(samp_matrix), value = TRUE)
pi_cols     <- grep("^pi\\[[0-9]+\\]$", colnames(samp_matrix), value = TRUE)

lambda_samples <- samp_matrix[, lambda_cols]
pi_samples     <- samp_matrix[, pi_cols]
correction_samples <- 1 / pi_samples
z_samples <- lambda_samples * pi_samples

# Extract summaries of key outcomes
y_summary   <- summarize_draws(samp_matrix, "y")
z_summary   <- summarize_draws(samp_matrix, "z")
pi_summary  <- summarize_draws(samp_matrix, "pi")

# get births in matrix form
b_i <- nimble_data$live_births
b_matrix <- matrix(rep(b_i, each = nrow(lambda_samples)), ncol = ncol(lambda_samples), byrow = FALSE)
MMR_samples <- (lambda_samples / b_matrix) * 100000

## Summaries of main parameters
lambda_df     <- summarize_draws(lambda_samples, "lambda")
pi_df         <- summarize_draws(pi_samples, "pi")
z_df          <- summarize_draws(z_samples, "z")
MMR_df        <- summarize_draws(MMR_samples, "MMR")
correction_df <- summarize_draws(correction_samples, "corr")

## combining to one output
summary_df <- reduce(
  list(lambda_df, z_df, pi_df, MMR_df, correction_df),
  full_join,
  by = "LGA"
) 
summary_df$LGA <- LGA_names
write.csv(summary_df,file = paste0(model_out_dir,"/summary_stats.csv"),row.names = F)

# a full dataset that can be used to compare observed vs output vs covariates etc.
combined_df <- left_join(summary_df,data[,c("LGA","repRate","ANC_ref","lb","deaths","MMR","anyEd",
                                           "tt_mean_unweighted","ANC.1st")])
combined_df$ANC_ref_scaled <- nimble_data$ancRef
combined_df$ANC_4_scaled <- nimble_data$anc
combined_df$tt_scaled <- nimble_data$travel
combined_df$educ_scaled <- nimble_data$educ

#some other statistics that would be helpful
combined_df$MMR_diff <- abs(combined_df$MMR_median - combined_df$MMR)
# observed 95% CI
combined_df <- combined_df %>%
  mutate(
    MMR_obs_SE = sqrt(deaths) / lb * 100000,
    MMR_obs_lower_95 = MMR - 1.96 * MMR_obs_SE,
    MMR_obs_upper_95 = MMR + 1.96 * MMR_obs_SE
  )

#joining with shapefile
map_out <- left_join(kdn_shp,combined_df,by=c("NAME_2"="LGA"))
combined_df$MMR_obs_CiW <- combined_df$MMR_obs_upper_95 - combined_df$MMR_obs_lower_95
combined_df$MMR_CiW <- combined_df$MMR_upper_95 - combined_df$MMR_lower_95

###### Plots of summary statistics ## 
dir.create(paste0(model_out_dir,"/plots/"))

pZobsVExpec <- ggplot(combined_df, aes(x = reorder(factor(LGA), lambda_median))) +
                  geom_point(aes(y = lambda_median, color = "Estimated")) +
                  geom_errorbar(aes(ymin = lambda_lower_95, ymax = lambda_upper_95), width = 0.2) +
                  geom_point(aes(y = deaths, color = "Observed")) +
                  labs(x = "LGA", y = "Expected Maternal Deaths (λ[i])",
                       title = "Posterior Estimates of Adjusted vs. Observed Maternal Deaths by LGA",
                       color = "Type") +
                  theme_minimal() +
                  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                  scale_color_manual(values = c("Estimated" = "black", "Observed" = "red"),name = NULL)

png(filename = paste0(model_out_dir,"/plots/PostDeaths vs Observed.png"),height=900,width=1200)
pZobsVExpec
dev.off()

pMMRobsVExpec <- ggplot(combined_df, aes(x = reorder(factor(LGA), MMR_diff))) +
  geom_point(aes(y = MMR_median, color = "Estimated")) +
  geom_errorbar(aes(ymin = MMR_lower_95, ymax = MMR_upper_95), width = 0.2) +
  geom_point(aes(y = MMR, color = "Observed")) +
  labs(x = "LGA", y = "MMR",
       title = "Posterior Estimates of Adjusted vs. Observed MMR by LGA",
       color = "Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c("Estimated" = "black", "Observed" = "red"),name = NULL)

png(filename = paste0(model_out_dir,"/plots/PostMMR vs Observed.png"),height=900,width=1200)
pMMRobsVExpec
dev.off()

pMMRCrIWidth <- ggplot(combined_df,aes(y=MMR_CiW, x = MMR_obs_CiW)) +
  geom_point() +
  labs(x = "Observed MMR 95% CI Width", y = "Posterior MMR CrI Width",
       title = "Comparison of uncertainty post modeling") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  geom_text_repel(data = filter(combined_df, MMR_CiW > 2000), aes(label = LGA)) +
  theme_minimal() 
png(filename = paste0(model_out_dir,"/plots/Uncertainity.png"),height=900,width=1200)
pMMRCrIWidth
dev.off()

pRepAdjust <- ggplot(combined_df, aes(x = reorder(factor(LGA), pi_median), y = pi_median)) +
  geom_point() +
  geom_errorbar(aes(ymin = pi_lower_95, ymax = pi_upper_95), width = 0.2) +
  labs(x = "LGA", y = "Reporting Rate (pi)",
       title = "Estimated Median Posterior Reporting Rate by LGA") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

png(filename = paste0(model_out_dir,"/plots/ReportingPost.png"),height=900,width=1200)
pRepAdjust
dev.off()

###### Maps of summary statistics #### 
dir.create(paste0(model_out_dir,"/maps/"))

#High burden LGAs
pLambdamap <- ggplot(map_out) +
  geom_sf(aes(fill = lambda_median )) +
  scale_fill_gradientn(
    colors = c("#7FB77E", "#F7DC6F", "#E74C3C"),
    name = "Deaths") +
  geom_sf_text(aes(label = NAME_2 ), size = 4, check_overlap = TRUE,color="white",fontface ="bold") +
  theme_void() +
  labs(fill = "Deaths",
       title = "Posterior Median Maternal Deaths (Adjusted for Under-reporting)")

pMMR <- ggplot(map_out) +
  geom_sf(aes(fill = MMR_median )) +
  scale_fill_gradientn(
    colors = c("#7FB77E", "#F7DC6F", "#E74C3C"),
    name = "MMR") +
  geom_sf_text(aes(label = NAME_2 ), size = 4, check_overlap = TRUE,color="white",fontface ="bold") +
  theme_void() +
  labs(fill = "MMR (per 100k LB)",
       title = "Posterior median MMR (Adjusted for Under-reporting)")

pMMRObs <- ggplot(map_out) +
  geom_sf(aes(fill = MMR )) +
  scale_fill_gradientn(
    colors = c("#7FB77E", "#F7DC6F", "#E74C3C"),
    name = "MMR") +
  geom_sf_text(aes(label = NAME_2 ), size = 4, check_overlap = TRUE,color="white",fontface ="bold") +
  theme_void() +
  labs(fill = "MMR (per 100k LB)",
       title = "Observed MMR")

pRiskBurdenMap <- pLambdamap + pMMR + pMMRObs +plot_layout(ncol = 2)

png(filename = paste0(model_out_dir,"/maps/Map_MMR_Deaths.png"),height=900,width=1200)
pRiskBurdenMap
dev.off()


pRepMap <- ggplot(map_out) +
  geom_sf(aes(fill = pi_median )) +
  scale_fill_viridis_c() +
  geom_sf_text(aes(label = NAME_2 ), size = 4, check_overlap = TRUE,color="white",fontface ="bold") +
  theme_void() +
  labs(fill = "Proportion",
       title = "Posterior Median Reporting Probability")
png(filename = paste0(model_out_dir,"/maps/Map_RepRate.png"),height=900,width=1200)
pRepMap
dev.off()


##### Comparison of incidence/reporting to predictors 
pRepPred <- ggplot(combined_df, aes(x = ANC_ref_scaled, y = pi_median)) +
              geom_point() +
              geom_smooth(method = "lm", se = FALSE, color = "blue") +
              labs(x = "ANC Reporting Covariate (Scaled)",
                   y = "Posterior Median of π[i]",
                   title = "Correlation between Reporting Predictors and Covariates") +
              geom_text_repel(data = filter(combined_df, pi_median < 0.6), aes(label = LGA)) +
              theme_minimal()

png(filename = paste0(model_out_dir,"/plots/ANC Ref vs Reporting.png"),height=900,width=1200)
pRepPred
dev.off()




#### Spatial vs unstructured effects for incidence
# Get phi and theta column names
phi_cols <- grep("^phi\\[[0-9]+\\]$", colnames(post_mat), value = TRUE)
theta_cols <- grep("^theta\\[[0-9]+\\]$", colnames(post_mat), value = TRUE)

# Extract posterior samples & summarize
phi_samples <- post_mat[, phi_cols]
theta_samples <- post_mat[, theta_cols]
tau_samples <- post_mat[, "tau"]

phi_df <- summarize_draws(phi_samples, "phi")
phi_df$LGA <- LGA_names
theta_df <- summarize_draws(theta_samples, "theta")
theta_df$LGA <- LGA_names

## add to output data
combined_df <- combined_df %>%
  left_join(phi_df, by = "LGA") %>%
  left_join(theta_df, by = "LGA")

map_out <- map_out  %>%
  left_join(phi_df, by = c("NAME_2"="LGA")) %>%
  left_join(theta_df, by = c("NAME_2"="LGA")) 


### Proportion of variation explained by spatial random effect
Var_phi <- var(combined_df$phi_mean)
Var_theta <- var(combined_df$theta_mean)
Var_spatial <- Var_phi + Var_theta
p_phi <- Var_phi / Var_spatial

## posteriors of phi & theta
phi_df <- as.data.frame(phi_samples)
phi_df$draw <- seq_len(nrow(phi_df))

theta_df <- as.data.frame(theta_samples)
theta_df$draw <- seq_len(nrow(theta_df))

theta_df <- as.data.frame(theta_samples)
theta_df$draw <- seq_len(nrow(theta_df))

# Reshape to long format for both phi and theta
phi_long <- phi_df %>%
  pivot_longer(
    cols = -draw,  # exclude 'draw' from pivoting
    names_to = "param",
    values_to = "phi_draw"
  ) %>%
  mutate(LGA = as.integer(gsub("phi\\[|\\]", "", param)))
phi_long <- left_join(phi_long,LGA_names2,by="LGA")

theta_long <- theta_df %>%
  pivot_longer(-draw, names_to = "param", values_to = "theta_draw") %>%
  mutate(LGA = as.integer(gsub("theta\\[|\\]", "", param)))
theta_long <- left_join(theta_long,LGA_names2,by="LGA")


pPhiPOst <- ggplot(phi_long, aes(x = phi_draw)) +
              geom_density(fill = "steelblue", alpha = 0.6) +
              facet_wrap(~ name, scales = "free", ncol = 6) +
              geom_vline(xintercept = 0, linetype = "dashed") +
              labs(title = "Posterior Distributions of Structured Spatial Effects (φ[i])",
                   x = "φ[i]", y = "Density") +
              theme_minimal()

png(filename = paste0(model_out_dir,"/plots/Spatial_RE_post.png"),height=900,width=1200)
pPhiPOst
dev.off()

pThetaPOst <- ggplot(theta_long, aes(x = theta_draw)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
  facet_wrap(~ name, scales = "free", ncol = 6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Posterior Distributions of Unstructured Spatial Effects (Theta[i])",
       x = "theta[i]", y = "Density") +
  theme_minimal()
png(filename = paste0(model_out_dir,"/plots/UnstructuredSpatial_RE_post.png"),height=900,width=1200)
pThetaPOst
dev.off()

pTauPost <- ggplot(data.frame(tau = tau_samples), aes(x = tau)) +
                geom_density(fill = "tomato", alpha = 0.6) +
                labs(title = "Posterior Distribution of τ (ICAR Precision)",
                     x = "τ", y = "Density") +
                theme_minimal()
png(filename = paste0(model_out_dir,"/plots/tau_post.png"),height=900,width=1200)
pTauPost
dev.off()


pPhiMap <- ggplot(map_out) +
  geom_sf(aes(fill = phi_median), color = "gray30") +
  scale_fill_gradient2(low = "darkred", mid = "white", high = "darkgreen", midpoint = 0) +
  labs(title = "Structured Spatial Effect (φ[i])") +
  theme_void()
 ## dark red is where MMR is higher than expected
 ## dark green is where MMR is lower than expected because of sptial
png(filename = paste0(model_out_dir,"/maps/Map_phi.png"),height=900,width=1200)
pPhiMap
dev.off()

pThetaMap <- ggplot(map_out) +
  geom_sf(aes(fill = theta_median), color = "gray30") +
  scale_fill_gradient2(
    low = "steelblue",
    mid = "white",
    high = "firebrick",
    midpoint = 0,
    name = "θ[i]"
  ) +  labs(title = "Unstructured Spatial Effect (theta[i])") +
  theme_void()
png(filename = paste0(model_out_dir,"/maps/Theta_phi.png"),height=900,width=1200)
pThetaMap
dev.off()
#Negative residuals (blue): model pulled mortality down
# neutral/no residual (white): model fit as expected
# Positive residuals (red): model pushed mortality up


