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
library(tidyverse)
library(ggplot2)
library(patchwork)


rm(list=ls())

source()

outdir <- "C:/Users/anumi/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub/CBHMIS_SAE_MMR/output/"
cbhmis <- readxl::read_xlsx("C:/Users/anumi/OneDrive - Bill & Melinda Gates Foundation/HDS - Health and Data Systems/CBHMIS CHAI Nigeria/Source data/Kaduna_Final Data for Analysis.xlsx",sheet="Merged Data")
names(cbhmis) <- gsub(pattern = " ",replacement = ".",x=names(cbhmis))

#light data cleaning 
cbhmis <- cbhmis[cbhmis$LGA!="0",]
cbhmis$LGA <- ifelse(cbhmis$LGA=="Birnin_gwari","Birnin_Gwari",cbhmis$LGA)
cbhmis$LGA <- ifelse(cbhmis$LGA=="Jema'a","Jemaa",cbhmis$LGA)
cbhmis$LGA <- ifelse(cbhmis$LGA=="Kaduna_north","Kaduna_North",cbhmis$LGA)
cbhmis$LGA <- ifelse(cbhmis$LGA=="Kaduna_south","Kaduna_South",cbhmis$LGA)
cbhmis$LGA <- ifelse(cbhmis$LGA=="Sabon_gari","Sabon_Gari",cbhmis$LGA)
cbhmis$LGA <- ifelse(cbhmis$LGA=="Zangon_kataf","Zango_Kataf",cbhmis$LGA)
cbhmis$LGA <- ifelse(cbhmis$LGA=="Makarfi","Markafi",cbhmis$LGA)


#subsetting to inclusion period
cbhmis$time <- as.Date(paste0(cbhmis$Month.of.reporting, cbhmis$Year.of.reporting,"1"), format="%B%Y%d")
cbhmis <- cbhmis[cbhmis$time >= as.Date("2023-10-01") & cbhmis$time < as.Date("2025-04-01"),]
cbhmis$submission_time <- cbhmis$`_submission_time`
cbhmis <- cbhmis %>% dplyr::select("LGA","Ward","time","Total.number.of.CVs.in.the.ward","Number.of.CVs.that.reported",
                        "No.of.ANC.referral.in.the.last.one.month","Total.number.of.live.births","Number.of.maternal.deaths",
                        "submission_time")
cbhmis$id <- 1:nrow(cbhmis)

#number of wards per LGA -- seems to be some duplicates 
ward_counts <- aggregate(Ward ~ LGA, data = cbhmis, FUN = function(x) length(unique(x)))
ward_time_counts <- aggregate(id ~ LGA  + Ward + time, data = cbhmis, FUN = length)  ## number of duplicate ward/month
names(ward_time_counts)[4] <- "count" 
ward_time_counts[ward_time_counts$count>1,]
write.csv(ward_time_counts[ward_time_counts$count>1,],paste0(outdir,"data checks/duplicates.csv"))


#only keeping last submission
cbhmis_clean <-  cbhmis %>%
  group_by(LGA, Ward, time)  %>%
  slice_max(order_by = submission_time, n = 1, with_ties = FALSE) %>%
  ungroup()

#check for duplicates again
ward_time_counts_clean <- aggregate(id ~ LGA  + Ward + time, data = cbhmis_clean, FUN = length)  ## number of duplicate ward/month
names(ward_time_counts_clean)[4] <- "count" 
ward_counts_clean <- aggregate(id ~ LGA  + Ward , data = cbhmis_clean, FUN = length)  ## number of duplicate ward/month
names(ward_counts_clean)[3] <- "count" 
hist(ward_counts_clean$count) ## shouldn't have more than 18 obs per ward given the months we've included
write.csv(ward_counts_clean,paste0(outdir,"data checks/num_obs_per_ward.csv"),row.names = F)

# some manual cleaning -- total number of CVs and rep CVs clearly getting enetered in one cell here
cbhmis_clean$Total.number.of.CVs.in.the.ward <- ifelse(cbhmis_clean$Total.number.of.CVs.in.the.ward>1000, 
                                                 as.numeric(substr(as.character(cbhmis_clean$Number.of.CVs.that.reported), 1, 2)),
                                                 cbhmis_clean$Total.number.of.CVs.in.the.ward)

# ASSUMPTION: If the number of CVs per ward is 3 digits, we're going to take the first two & one ward seems kind of off
cbhmis_clean$Total.number.of.CVs.in.the.ward <- ifelse(cbhmis_clean$Total.number.of.CVs.in.the.ward>100, 
                                                       as.numeric(substr(as.character(cbhmis_clean$Number.of.CVs.that.reported), 1, 2)),
                                                       cbhmis_clean$Total.number.of.CVs.in.the.ward)
cbhmis_clean$Total.number.of.CVs.in.the.ward[cbhmis_clean$Ward=="Turawa" & cbhmis_clean$time=="2024-11-01"] <- 9

names(cbhmis_clean) <- c("LGA", "Ward","time","totalCVs","repCVs","ANC_ref","lb","deaths","sub_time","id")



#### A: Spatial data for adjacency matrix in model #### 
#Nigeria shape file from GADM
nga_shp <- st_read("C:/Users/anumi/OneDrive - Bill & Melinda Gates Foundation/Data/shape_files/NGA/nga_admbnda_adm2_osgof/nga_admbnda_adm2_osgof_20170222.shp")
kdn_shp <- nga_shp[nga_shp$admin1Name=="Kaduna",]
kdn_shp$admin2Name <- gsub(pattern = "-",replacement = "_", kdn_shp$admin2Name)
kdn_shp$admin2Name <- gsub(pattern = " ",replacement = "_", kdn_shp$admin2Name)
kdn_shp$admin2Name <- gsub(pattern = "'",replacement = "", kdn_shp$admin2Name)

st_crs(kdn_shp) ## checking projection of shape file to match GPS coordinates
mapview(kdn_shp) ## view map of Kaduna

# create adjacency matrix  -- Don't actually need the matrix for stan but just the nb file
nb <- poly2nb(kdn_shp)
adj_matrix <- nb2mat(nb, style = "B", zero.policy = TRUE)
id_col <- kdn_shp$admin2Name
rownames(adj_matrix) <- id_col
colnames(adj_matrix) <- id_col
write.csv(adj_matrix,paste0(outdir,"/adj_mat.csv"))


#### B: Analysis of Under-reporting ####
## Ward-level reporting rates over time
cbhmis_clean$repRate <- cbhmis_clean$repCVs/cbhmis_clean$totalCVs
table(cbhmis_clean$repRate,exclude = NULL)

cbhmis_clean <- cbhmis_clean %>% group_by(LGA) %>%
  mutate(Ward_color_group = as.factor((as.numeric(factor(Ward)) - 1) %% 10 + 1))

p1 <- ggplot(cbhmis_clean,aes(x=time,y=repRate,group=Ward,color=Ward_color_group)) + facet_wrap(~LGA) + 
    geom_point() + geom_line() + theme_bw() +   theme(legend.position = "none")
p1


## aggregating to LGA level
cbhmis_lga <- cbhmis_clean[,!names(cbhmis_clean) %in% c("Ward","repRate","Ward_color_group")] %>% 
                  group_by(time, LGA) %>% summarise(across(everything(), sum))

cbhmis_lga$repRate <- 100*cbhmis_lga$repCVs/cbhmis_lga$totalCVs
p.RR <- ggplot(cbhmis_lga, aes(x = time, y = repRate)) +
  geom_line() +
  geom_point() +
  facet_wrap(~LGA) +
  ylab("Reporting Rate") +
  scale_x_date(date_labels = "%b %y") +
  theme_bw() +
  ylim(0, 100) +
  theme(
    strip.text = element_text(size = 14),      # Facet labels
    axis.title = element_text(size = 14),      # Axis labels
    axis.text = element_text(size = 12),       # Tick labels
    plot.title = element_text(size = 16),      # Title (if you add one)
    legend.text = element_text(size = 12),     # Legend text (if present)
    legend.title = element_text(size = 14)     # Legend title (if present)
  )
png(filename = paste0(outdir,"/plots/repRateMonthly.png"),height=800,width=1400)
p.RR
dev.off()

cbhmis_lga_rep <- cbhmis_lga[,names(cbhmis_lga) %in% c("LGA","repRate")] %>% 
  group_by(LGA) %>% summarise(across(everything(), mean)) %>%
  rename("avgRegRate"="repRate")
  

repRate <- left_join(kdn_shp,cbhmis_lga_rep,by=c("admin2Name"="LGA"))
pRRmap <- ggplot(repRate) +
  geom_sf(aes(fill = avgRegRate)) +
  scale_fill_viridis_c() +
  geom_sf_text(aes(label = admin2Name ), size = 4, check_overlap = TRUE,color="white",fontface ="bold") +
  theme_void() +
  labs(fill = "Average reporting rate",
       title = "Reporting Rate for LGA averaged over monthly rates")
png(filename = paste0(outdir,"/plots/repRateMap.png"),height=900,width=1200)
pRRmap
dev.off()



#### C: Prep covariate data #### 

##### C1: Referral data for under-reporting #####
## We will use ANC referral rates as our proxy of under-reporing




