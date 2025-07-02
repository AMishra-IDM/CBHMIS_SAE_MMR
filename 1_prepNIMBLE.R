#################################################################
##   LGA estimates of community MMR                            ##
##   1. Prep_NIMBLE                                            ##
##   Purpose: Using cleaned data for model prepare for NIMBLE  ##
##            input                                            ##
##   Author: Anu Mishra                                        ##  
##   Created 6/9/25                                            ##
#################################################################
library(sf)
library(spdep)
library(ggrepel)

rm(list=ls())
set.seed(105)
setwd("C:/Users/anumi/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub/CBHMIS_SAE_MMR")
outdir <- "C:/Users/anumi/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub/CBHMIS_SAE_MMR/output/"


#### 1: Read in data ####
modData <- read.csv("cbhmis_data_for_model.csv")

## define Outcome and offset
# ------------------------
z <- modData$deaths
live_births <- modData$lb
N <- length(z)  # number of LGAs


#### 2: Center and covariates  ####
incidence_covars <- scale(modData[, c("anyEd", "tt_mean_unweighted", "ANC.1st")],
                          center = TRUE, scale = FALSE)
ggplot(modData, aes(x = anyEd, y = incidence_covars[,1], label = LGA)) +
  geom_point() +
  geom_text_repel(size = 3) +
  theme_minimal()

ggplot(modData, aes(x = tt_mean_unweighted, y = incidence_covars[,2], label = LGA)) +
  geom_point() +
  geom_text_repel(size = 3) +
  theme_minimal()

ggplot(modData, aes(x = ANC.1st, y = incidence_covars[,3], label = LGA)) +
  geom_point() +
  geom_text_repel(size = 3) +
  theme_minimal()


# Reporting model covariate: center
reporting_covars <- scale(modData[, "ANC_ref_scaled"],
                          center = TRUE, scale = FALSE)
ggplot(modData, aes(x = ANC_ref_scaled, y = reporting_covars[,1], label = LGA)) +
  geom_point() +
  geom_text_repel(size = 3) +
  theme_minimal()

#### 3: Create spatial inputs  ####
nga_shp <- st_read("raw data/gadm41_NGA_shp/gadm41_NGA_2.shp")
kdn_shp <- nga_shp[nga_shp$NAME_1=="Kaduna",]
kdn_shp$NAME_2 <- gsub(pattern = " ",replacement = "_", kdn_shp$NAME_2)
kdn_shp$NAME_2 <- gsub(pattern = "'",replacement = "", kdn_shp$NAME_2)

nb <- poly2nb(kdn_shp)
lw <- nb2listw(nb, style = "B", zero.policy = TRUE)  # binary style

# Extract NIMBLE structures
adj <- unlist(nb)
num <- sapply(nb, length)
weights <- rep(1, length(adj))
l_adj <- length(adj)


##### 5: Assemble NIMBLE inputs ##### 

nimble_constants <- list(
  n = N,
  R = N,
  adj = adj,
  num = num,
  weights = weights,
  l_adj = l_adj
)

nimble_data <- list(
  z = z,
  educ = incidence_covars[, 1],
  travel = incidence_covars[, 2],
  anc4 = incidence_covars[, 3],
  ancRef = reporting_covars[, 1],
  live_births = live_births
)

### Intial values -- REALLY SHOULD DO ONE PER CHAIN
nimble_inits <- list(
  a = c(0, 0, 0, 0),  # intercept + 3 incidence covariates
  b = c(2, 0),         # intercept + 1 reporting covariate
  sigma = 0.5,
  epsilon = 0.5,
  nu = 0.5,
  phi = rnorm(N, 0, 1),
  theta = rnorm(N, 0, 0.5),
  gamma = rnorm(N, 0, 0.5)
)


#### 3: Save NIMBLE data  ####
save(nimble_constants, nimble_data, nimble_inits, file = "nimble_inputs_final.RData")
