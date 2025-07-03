rm(list=ls())

library(tidyverse)
library(ggplot2)
library(nimble)
library(coda) # For manipulation of MCMC results.
library(sf)
library(gridExtra)
library(patchwork)
library(ggrepel)
library(corrplot)


source("~/GitHub/CBHMIS_SAE_MMR/0_functions.R", echo=TRUE)
setwd("C:/Users/anumi/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub/CBHMIS_SAE_MMR")

model_date1 <- "2025-07-02_no_ANCRef"
model_output <- paste0("model_out_",model_date1,".rds")  ## which version of the model you want to read in
model_out_dir <- paste0("model runs/model_",model_date1,"/")

model_org <-  "2025-07-02"
model_output_orig <- paste0("model_out_",model_org,".rds")  ## which version of the model you want to read in
model_out_dir_orig <- paste0("model runs/model_",model_org,"/")

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

load("nimble_inputs_final.RData")
nimble_data_LGA <- nimble_data[sapply(nimble_data, length) == 23]
covariate_data <- as.data.frame(nimble_data_LGA)
covariate_data$LGA <- LGA_names

#### Loading data ####
## model output
full_out <- readRDS(paste0(model_out_dir,"/",model_output))
samples <- full_out$samples

full_out_orig <- readRDS("C:/Users/anumi/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub/CBHMIS_SAE_MMR/model runs/model_2025-07-02/model_out_2025-07-02.rds")
samples_orig <- full_out_orig$samples


post_mat <- as.matrix(samples)
post_mat_orig <- as.matrix(samples_orig)


##### 
pi_names <- grep("^pi\\[[0-9]+\\]$", colnames(mat_with), value = TRUE)

pi_with <- colMeans(post_mat_orig[, pi_names])
pi_without <- colMeans(post_mat[, pi_names])


pi_df <- tibble(
  LGA = LGA_names,
  pi_with_covariate = pi_with,
  pi_without_covariate = pi_without,
  delta = pi_without_covariate - pi_with_covariate
)

ggplot(pi_df, aes(x = reorder(LGA, delta), y = delta)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Change in πᵢ After Removing ANC_ref_scaled",
       x = "LGA", y = "πᵢ (without) − πᵢ (with)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

### correlations between the coviarates
cor_df <- covariate_data %>%
  select(ancRef, educ, travel, anc, mmr) %>%
  cor()


corrplot(cor_df, method = "color", addCoef.col = "black", number.cex = 0.8)
ggplot(mapping = aes(x=scale(data$ANC_ref_scaled,center = TRUE,scale = TRUE),y=data$MMR)) + geom_point() + geom_smooth(method="lm")

       