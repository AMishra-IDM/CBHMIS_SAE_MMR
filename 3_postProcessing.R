#################################################################
##   LGA estimates of community MMR                            ##
##   4. Post-processing                                        ##
##   Purpose: Fit checks and post-processing of raw model      ##
##            model output data                                ##
##   Author: Anu Mishra                                        ##  
##   Created 6/9/25                                            ##
#################################################################

library(tidyverse)
library(ggplot2)
library(nimble)
library(coda) # For manipulation of MCMC results.
library(sf)

source("~/GitHub/CBHMIS_SAE_MMR/0_functions.R", echo=TRUE)
setwd("C:/Users/anumi/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub/CBHMIS_SAE_MMR")
outdir <- "C:/Users/anumi/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub/CBHMIS_SAE_MMR/output/"

#### Loading data ####
## model output
full_out <- readRDS("output/model output/samples_full_output.RDS")
coda_samples <- mcmc.list(as.mcmc(full_out$samples$chain1), as.mcmc(full_out$samples$chain2),
                          as.mcmc(full_out$samples$chain3), as.mcmc(full_out$samples$chain4))  # for easy plotting

## shape file 
nga_shp <- st_read("raw data/gadm41_NGA_shp/gadm41_NGA_2.shp")
kdn_shp <- nga_shp[nga_shp$NAME_1=="Kaduna",]
kdn_shp$NAME_2 <- gsub(pattern = " ",replacement = "_", kdn_shp$NAME_2)
kdn_shp$NAME_2 <- gsub(pattern = "'",replacement = "", kdn_shp$NAME_2)

## CBHMIS data
data <- read.csv("cbhmis_data_for_model.csv")
LGA_names <- data$LGA

#### Checking convergence ####
# Check chains for convergence

pdf(file = "output/model output/convergence checks/trace_plots.pdf",height=8,width=12)
traceplot(coda_samples[,c('a[1]','a[2]','a[3]','a[4]',   ## education parameters
                         'b[1]','b[2]',                 ## under-reporting parameters
                         'epsilon','sigma','nu')])     ## epsilon is unstructured noise (s.d) for the reporting
                                                        ## sigma is unstructured noise (s.d.) for the incidence
                                                        ## precision for the structured spatial effect
dev.off()
                  
gelman.diag(coda_samples[,c('a[1]','a[2]','a[3]','a[4]',   
                                 'b[1]','b[2]',                 
                                 'epsilon','sigma',"nu")])



#### Producing summary estimates ####
mcmc_combined <- as.mcmc(do.call(rbind, full_out$samples))

# Calculate posterior summary
lambda_draws <- mcmc_combined[, grep("^lambda\\[", colnames(mcmc_combined), value = TRUE)]
pi_draws <- mcmc_combined[, grep("^pi\\[", colnames(mcmc_combined), value = TRUE)]
mmr_samples <- sweep(lambda_draws, 2, data$lb, "/") * 1e5

lambda_summary <- summarize_matrix(lambda_draws)
names(lambda_summary) <- paste0("lambda.",names(lambda_summary))

mmr_summary <- summarize_matrix(mmr_samples)
names(mmr_summary) <- paste0("mmr.",names(mmr_summary))

pi_summary <- summarize_matrix(pi_draws)
names(pi_summary) <- paste0("rep.",names(pi_summary))

all_summary <- cbind(mmr_summary,pi_summary,lambda_summary)
all_summary$LGA <- data$LGA
names(all_summary) <- gsub(pattern = ".50%|.2.5%|.97.5%",replacement = "",names(all_summary))
row.names(all_summary) <- NULL

all_summary_long <- all_summary %>%
  pivot_longer(
    cols = -LGA,
    names_to = c("metric", ".value"),
    names_pattern = "(.*)\\.(mean|sd|median|l95|u95)"
  )


### Overall summary 
pSummary <- ggplot(all_summary_long, aes(x =LGA, y = median)) + 
              facet_wrap(~metric,scales="free_y") +
              geom_point() +
              geom_errorbar(aes(ymin = l95, ymax = u95)) + theme_bw() +
          theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

### Maps 
kdn_shp_out <- left_join(kdn_shp,all_summary,by=c("NAME_2"="LGA"))
pEstMMR <- ggplot(kdn_shp_out) +
                geom_sf(aes(fill = mmr.median)) +
                scale_fill_viridis_c() +
                geom_sf_text(aes(label = NAME_2 ), size = 4, check_overlap = TRUE,color="white",fontface ="bold") +
                theme_void() +
                labs(fill = "MMR (per 100k LB)",
                     title = "Estimated Median MMR")

pEstRep <- ggplot(kdn_shp_out) +
  geom_sf(aes(fill = rep.median)) +
  scale_fill_viridis_c() +
  geom_sf_text(aes(label = NAME_2 ), size = 4, check_overlap = TRUE,color="white",fontface ="bold") +
  theme_void() +
  labs(fill = "Reporting Rate (%)",
       title = "Estimated Median Rep Rate")


#### Fitted vs predictive checks ####