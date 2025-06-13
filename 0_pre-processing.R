#################################################################
##   LGA estimates of community MMR                            ##
##   0. Pre-processing data                                    ##
##   Purpose: To get all covariate, outcome, and spatial data  ##
##            cleaned and prepped for model                    ##
##   Author: Anu Mishra                                        ##  
##   Created 6/9/25                                            ##
#################################################################
library(sf)
library(spdep)
library(mapview)

outdir <- "C:/Users/anumi/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub/CBHMIS_SAE_MMR/"
cbhmis <- readxl::read_xlsx("C:/Users/anumi/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub/CBHMIS_SAE_MMR/Kaduna TL dataset _ Jul2024.xlsx",sheet=2)

#### A: Spatial data for adjacency matrix in model #### 
#Nigeria shape file from GADM
nga_shp <- st_read("C:/Users/anumi/OneDrive - Bill & Melinda Gates Foundation/Data/shape_files/NGA/nga_admbnda_adm2_osgof/nga_admbnda_adm2_osgof_20170222.shp")
kdn_shp <- nga_shp[nga_shp$admin1Name=="Kaduna",]

st_crs(kdn_shp) ## checking projection of shape file to match GPS coordinates
mapview(kdn_shp) ## view map of Kaduna

# create adjacency matrix  -- Don't actually need the matrix for stan but just the nb file
nb <- poly2nb(kdn_shp)
adj_matrix <- nb2mat(nb, style = "B", zero.policy = TRUE)
id_col <- kdn_shp$admin2Name
rownames(adj_matrix) <- id_col
colnames(adj_matrix) <- id_col
write.csv(adj_matrix,paste0(outdir,"/adj_mat.csv"))


#### B: Prep covariate data #### 

##### B1: Referral data for under-reporting #####
