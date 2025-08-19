#########################################
##     Correcting Under-Reporting      ##
##            Master Script            ##  
#########################################

# Welcome to the supplementary materials accompanying
# 'A Hierarchical Framework for Correcting Under-Reporting in Count Data'.
# These scripts are intended to be used to reproduce results found in the paper
# and to illustrate how the model can be used in practice.

# If you would just like to reproduce results, just source this script and the
# figures will be saved as PDF files in your working directory. Alternatively,
# use the relelvant commands below.

# R devtools may need to be installed to use NIMBLE.
# Please ensure all packages listed below are installed and up to date.

thin_multiplier <- 1 # Thinning less improves effective sample size but
# this code uses approximately 20GB of system memory when this is set to 1. If you 
# don't have enough memory increase this to 10, which will use approximately
# 8GB of system memory. This only affects the simulation experiments.

# To further reduce memory usage, you may need to reduce the number of chains
# and increase the thinning level in the TB model.

library(ggplot2) # For reproducing plots seen in the paper.
library(nimble) # For MCMC computation using NIMBLE.
library(coda) # For manipulation of MCMC results.
library(mgcv)
library(dplyr)
library(ngspatial)
library(sp) # For plotting the micro-regions.
library(spdep) # For computing the neighbourhood and adjancency objects.
library(maps) # For adding a map to the plots of Brazil.
library(viridis)
library(mapproj)

# Load in some functions.
source('Functions.R')

vp=viridis_pal()(20) # Colour pallette for plots.

seed <- 794637 # Seed for R and NIMBLE, 794637 is used for results in the article.

load('Data.RData') # Load in the data.

# Simulate some data.
source('Simulation.R')

# Run the sensitivity analysis.
source('Experiments.R')

# Run the tuberculosis analysis.
source('Tuberculosis.R')