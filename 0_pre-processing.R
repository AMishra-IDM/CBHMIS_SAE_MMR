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
library(terra)
library(lubridate)


rm(list=ls())
set.seed(105)

setwd("C:/Users/anumi/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub/CBHMIS_SAE_MMR")

outdir <- "C:/Users/anumi/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub/CBHMIS_SAE_MMR/output/"
cbhmis <- readxl::read_xlsx("raw data/Kaduna_Final Data for Analysis.xlsx",sheet="Merged Data")
names(cbhmis) <- gsub(pattern = " ",replacement = ".",x=names(cbhmis))

#light data cleaning 
cbhmis <- cbhmis[cbhmis$LGA!="0",]
cbhmis$LGA <- ifelse(cbhmis$LGA=="Birnin_gwari","Birnin_Gwari",cbhmis$LGA)
cbhmis$LGA <- ifelse(cbhmis$LGA=="Jema'a","Jemaa",cbhmis$LGA)
cbhmis$LGA <- ifelse(cbhmis$LGA=="Kaduna_north","Kaduna_North",cbhmis$LGA)
cbhmis$LGA <- ifelse(cbhmis$LGA=="Kaduna_south","Kaduna_South",cbhmis$LGA)
cbhmis$LGA <- ifelse(cbhmis$LGA=="Sabon_gari","Sabon_Gari",cbhmis$LGA)
cbhmis$LGA <- ifelse(cbhmis$LGA=="Zangon_kataf","Zangon_Kataf",cbhmis$LGA)


#subsetting to inclusion period
cbhmis$time <- as.Date(paste0(cbhmis$Month.of.reporting, cbhmis$Year.of.reporting,"1"), format="%B%Y%d")
cbhmis <- cbhmis[cbhmis$time >= as.Date("2023-10-01") & cbhmis$time < as.Date("2025-04-01"),]
cbhmis$submission_time <- cbhmis$`_submission_time`

#subsetting to a smaller set of variable
cbhmis <- cbhmis %>% dplyr::select("LGA","Ward","time","Total.number.of.CVs.in.the.ward","Number.of.CVs.that.reported",
                        "No.of.ANC.referral.in.the.last.one.month","No.of.FP.referral.in.the.last.one.month",
                        "No.of.Labour.&.delivery.referrals.in.the.last.one.month","No.of.PNC.referral.in.the.last.one.month",
                        "No.of.PPFP.referral.in.the.last.one.month","No.of.immunization.referral","No.of.Nutrition.referral",
                        "No.of.Integrated.Management.of.Childhood.Illness.(IMCI).referral",
                        "Total.number.of.live.births","Number.of.live.births.-.Males","Number.of.live.births.-.Females",
                        "Number.of.maternal.deaths",
                        "submission_time")
cbhmis$id <- 1:nrow(cbhmis)

## checking number of missing obs for a variable
colSums(is.na(cbhmis))


#number of wards per LGA -- seems to be some duplicates 
ward_counts <- aggregate(Ward ~ LGA, data = cbhmis, FUN = function(x) length(unique(x)))
ward_time_counts <- aggregate(id ~ LGA  + Ward + time, data = cbhmis, FUN = length)  ## number of duplicate ward/month
names(ward_time_counts)[4] <- "count" 
ward_time_counts[ward_time_counts$count>1,]
write.csv(ward_time_counts[ward_time_counts$count>1,],paste0(outdir,"raw data checks/duplicates.csv"))


#only keeping last submission
cbhmis_clean <-  cbhmis %>%
  group_by(LGA, Ward, time)  %>%
  slice_max(order_by = submission_time, n = 1, with_ties = FALSE) %>%
  ungroup()

## checking number of missing obs for a variable
colSums(is.na(cbhmis_clean))


#check for duplicates again
ward_time_counts_clean <- aggregate(id ~ LGA  + Ward + time, data = cbhmis_clean, FUN = length)  ## number of duplicate ward/month
names(ward_time_counts_clean)[4] <- "count" 
ward_counts_clean <- aggregate(id ~ LGA  + Ward , data = cbhmis_clean, FUN = length)  ## number of duplicate ward/month
names(ward_counts_clean)[3] <- "count" 
hist(ward_counts_clean$count) ## shouldn't have more than 18 obs per ward given the months we've included
write.csv(ward_counts_clean,paste0(outdir,"raw data checks/num_obs_per_ward.csv"),row.names = F)

# some manual cleaning -- total number of CVs and rep CVs clearly getting enetered in one cell here
cbhmis_clean$Total.number.of.CVs.in.the.ward <- ifelse(cbhmis_clean$Total.number.of.CVs.in.the.ward>1000, 
                                                 as.numeric(substr(as.character(cbhmis_clean$Number.of.CVs.that.reported), 1, 2)),
                                                 cbhmis_clean$Total.number.of.CVs.in.the.ward)

# ASSUMPTION: If the number of CVs per ward is 3 digits, we're going to take the first two & one ward seems kind of off
cbhmis_clean$Total.number.of.CVs.in.the.ward <- ifelse(cbhmis_clean$Total.number.of.CVs.in.the.ward>100, 
                                                       as.numeric(substr(as.character(cbhmis_clean$Number.of.CVs.that.reported), 1, 2)),
                                                       cbhmis_clean$Total.number.of.CVs.in.the.ward)
cbhmis_clean$Total.number.of.CVs.in.the.ward[cbhmis_clean$Ward=="Turawa" & cbhmis_clean$time=="2024-11-01"] <- 9

names(cbhmis_clean) <- c("LGA", "Ward","time","totalCVs","repCVs",
                         "ANC_ref","FP_ref","LD_ref","PNC_ref","PPFP_ref","Imm_ref","Nut_ref","ICMI_ref",
                         "lb","lb_male","lb_female",
                         "deaths","sub_time","id")



#### A: Spatial data for adjacency matrix in model #### 
#Nigeria shape file from GADM
nga_shp <- st_read("raw data/gadm41_NGA_shp/gadm41_NGA_2.shp")
kdn_shp <- nga_shp[nga_shp$NAME_1=="Kaduna",]
kdn_shp$NAME_2 <- gsub(pattern = " ",replacement = "_", kdn_shp$NAME_2)
kdn_shp$NAME_2 <- gsub(pattern = "'",replacement = "", kdn_shp$NAME_2)

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
#table(cbhmis_clean$repRate,exclude = NULL)

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
  

repRate <- left_join(kdn_shp,cbhmis_lga_rep,by=c("NAME_2"="LGA"))
pRRmap <- ggplot(repRate) +
  geom_sf(aes(fill = avgRegRate)) +
  scale_fill_viridis_c() +
  geom_sf_text(aes(label = NAME_2 ), size = 4, check_overlap = TRUE,color="white",fontface ="bold") +
  theme_void() +
  labs(fill = "Average reporting rate",
       title = "Reporting Rate for LGA averaged over monthly rates")
png(filename = paste0(outdir,"/plots/repRateMap.png"),height=900,width=1200)
pRRmap
dev.off()


####  C: Checking CBHMIS Referrral and MMR data      ####
##### C1: Cleaning referral data for under-reporting #####
## cleaning will happen on the monthly level rather than aggregate level

referral_vars <- c("ANC_ref","FP_ref","LD_ref","PNC_ref","PPFP_ref","Imm_ref","Nut_ref","ICMI_ref")

### NEEDS FOLLOW_UP: Use consistent population data once decided

## We load in external population data to scale the referral data
pop_wra <- rast("C:/Users/anumi/OneDrive - Bill & Melinda Gates Foundation/Data/Population data/NGA_population_v2_1_agesex/NGA_population_v2_1_agesex/NGA_population_v2_1_agesex_f15_49.tif")
boundary <- kdn_shp
boundary$wra_pop <- exact_extract(pop_wra, boundary, 'sum')

pop_dat <- boundary %>% 
  st_drop_geometry() %>% 
  as.data.frame() %>% dplyr::select(NAME_2, wra_pop)

## Checking referral data to see how it looks against reporting rate and for outliers
palette23 <- createPalette(23, seedcolors = c("#000000"))
names(palette23) <- unique(cbhmis_lga$LGA)  # Named vector for manual scale

scat.plots <- map(referral_vars, function(var) {
  ggplot(cbhmis_lga, aes(x = repRate, y = .data[[var]])) +
    geom_point(aes(color = LGA)) +
    geom_smooth(
      method = "loess",
      se = FALSE,
      color = "black",
      linetype = "dashed") +
    scale_color_manual(values = palette23) +
    theme_minimal() +
    labs(
      title = paste(var),
      x = "Reporting Rate",
      y = "Num. Referrals"
    )
})
ref_plot <- wrap_plots(scat.plots, ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")
png(filename = paste0(outdir,"/plots/referralplots.png"),height=900,width=1200)
ref_plot
dev.off()

## We will use ANC referral rates as our proxy of under-reporing. Code from here on
##  out will only focus on ANC

## lots of outliers again for Sabon Gari -- will look at exclusion based on (1) 3SD and
## (2)assuming 5% of women of WRA are pregnant and 70% are attending ANC (gives an upper bound)
outlier.sd <-  cbhmis_lga %>%
  group_by(LGA)  %>%
  summarise(threeSD = 3*sd(ANC_ref) + mean(ANC_ref, na.rm=T))

outlier.wra <- pop_dat
outlier.wra$preg_women <- outlier.wra$wra_pop*0.05
outlier.wra$anc_max <- outlier.wra$preg_women*0.7
outlier <- left_join(outlier.sd,outlier.wra,by=c("LGA"="NAME_2"))

pANC <- ggplot(cbhmis_lga) + 
         facet_wrap(~LGA,scales = 'free_y') + 
          geom_point(aes(x = repRate, y = ANC_ref)) + 
          geom_hline(data=outlier, aes(yintercept=threeSD, color = "3 SD Threshold")) + 
          geom_hline(data=outlier, aes(yintercept=anc_max, color = "Max ANC")) + 
         theme_bw() + xlim(c(0,100)) +   
         scale_color_manual(values = c("3 SD Threshold" = "blue", "Max ANC" = "green"),name="Outlier Rule") 
png(filename = paste0(outdir,"/plots/ANC_ref_outlier.png"),height=900,width=1200)
pANC
dev.off()

### There is not a best rule, 3SD looks generally less conservative, but the maximum logical
##    ANC rule seems to overly exclude the Kauru and Jaba data. Will go with this one for now

## join with CBHMIS data
cbhmis_lga <- dplyr::left_join(cbhmis_lga,outlier,by=c("LGA"))

## surpress those LGA obs where the obs are above threshold
cbhmis_lga$ANC_ref <- ifelse(cbhmis_lga$ANC_ref > cbhmis_lga$anc_max,NA,cbhmis_lga$ANC_ref)

pANC_remove <- ggplot(cbhmis_lga,aes(x = repRate, y = ANC_ref)) + 
  geom_point(aes(color = LGA)) + 
  scale_color_manual(values = palette23) +
  geom_smooth(method = "loess", se = FALSE, color = "black",linetype = "dashed") + 
  theme_bw() + xlim(c(50,100)) 
png(filename = paste0(outdir,"/plots/ANC_ref_remove.png"),height=900,width=1200)
pANC_remove
dev.off()

### NEEDS FOLLOW_UP: This is not a super strong relationship, so perhaps we need to look at FP -- for now use this one

##### C2: Checking monthly MMR data #####
sd_lines <- cbhmis_lga %>%
  group_by(LGA) %>%
  summarise(deaths_3sd = 3*sd(deaths, na.rm = TRUE),
            lb_3sd = 3*sd(lb, na.rm = TRUE) + mean(lb, na.rm=T),
            lb_4sd = 4*sd(lb, na.rm = TRUE) + mean(lb, na.rm=T),)

pMD <- ggplot(data=cbhmis_lga) + facet_wrap(~LGA) + ylab("Maternal Deaths") +
        geom_line(aes(x=time,y=deaths)) + geom_point(aes(x=time,y=deaths)) + theme_bw() + 
       geom_hline(data = sd_lines, aes(yintercept = deaths_3sd,color="3 SD"), 
             linetype = "dashed") #+ scale_color_manual(values = c("3 SD" = "red"),name="Outlier Rule") 


pLB <- ggplot(data=cbhmis_lga,aes(x=time,y=lb)) + facet_wrap(~LGA) + ylab("Live Births") +
  geom_line() + geom_point() + theme_bw() +
  geom_hline(data = sd_lines, aes(yintercept = lb_3sd,color="3 SD"), 
             linetype = "dashed") +
  geom_hline(data = sd_lines, aes(yintercept = lb_4sd,color="4 SD"), 
             linetype = "dashed") + 
    scale_color_manual(values = c("3 SD" = "red","4 SD" = "blue"),name="Outlier Rule") 

MMR_line <- pMD + pLB + plot_layout(ncol = 2)

png(filename = paste0(outdir,"/plots/MMR_line.png"),height=900,width=1200)
MMR_line
dev.off()

### NEEDS FOLLOW_UP: There are some outliers that we need to check but leaving them for now


##### C3: Aggregating data across time for final CBHMIS data #####
cbhmis_agg <- cbhmis_lga %>% 
  select(LGA,repRate,ANC_ref,lb,deaths,wra_pop) %>%
  group_by(LGA) %>% 
  summarise(repRate = mean(repRate,na.rm=TRUE),
            ANC_ref = mean(ANC_ref, na.rm=TRUE),
            lb = sum(lb,na.rm = T),
            deaths=sum(deaths, na.rm = T),
            wra = max(wra_pop))

cbhmis_agg$ANC_ref_scaled <- 1000*cbhmis_agg$ANC_ref/cbhmis_agg$wra
cbhmis_agg$MMR <- 100000*cbhmis_agg$deaths/cbhmis_agg$lb
kdn_shp <- left_join(kdn_shp,cbhmis_agg,by=c("NAME_2"="LGA"))


pMMRmap <- ggplot(kdn_shp) +
  geom_sf(aes(fill = MMR)) +
  scale_fill_viridis_c() +
  geom_sf_text(aes(label = NAME_2 ), size = 4, check_overlap = TRUE,color="white",fontface ="bold") +
  theme_void() +
  labs(fill = "MMR (per 100k LB)",
       title = "MMR for period of Oct '23-Mar '25")

png(filename = paste0(outdir,"/plots/MMR_map.png"),height=1000,width=1200)
pMMRmap
dev.off()


pANCrefmap <- ggplot(kdn_shp) +
  geom_sf(aes(fill = ANC_ref_scaled)) +
  scale_fill_viridis_c() +
  geom_sf_text(aes(label = NAME_2 ), size = 4, check_overlap = TRUE,color="white",fontface ="bold") +
  theme_void() +
  labs(fill = "ANC referrals (per 1000 WRA)",
       title = "Average ANC referrals over Oct '23-Mar '25 period")

png(filename = paste0(outdir,"/plots/ANCref_map.png"),height=1000,width=1200)
pANCrefmap
dev.off()


#### D: Pulling in External Predictor data ####
##### D1: Education levels from Kaduna HHS #####
educ <- read.csv("raw data/education_khhs.csv")
educ$LGA <- ifelse(educ$LGA=="MAKARFI","Makarfi",educ$LGA)
educ$LGA <- ifelse(educ$LGA=="ZANGON_KATAF","Zangon_Kataf",educ$LGA)

educ$LGA <- gsub(pattern = "_",replacement = " ",x = educ$LGA)
educ <- educ %>% select(LGA,ever_school,second.plus) %>%
         mutate(ever_school=100*ever_school,
                second.plus=100*second.plus,
                LGA = str_to_title(LGA, locale = "en"))
educ$LGA <- gsub(pattern = " ",replacement = "_",x = educ$LGA)
kdn_shp <- left_join(kdn_shp,educ,by=c("NAME_2"="LGA"))



pAnyEdmap <- ggplot(kdn_shp) +
  geom_sf(aes(fill = ever_school )) +
  scale_fill_viridis_c() +
  geom_sf_text(aes(label = NAME_2 ), size = 4, check_overlap = TRUE,color="white",fontface ="bold") +
  theme_void() +
  labs(fill = "Proportion",
       title = "Women who have attended any school")


pSecondmap <- ggplot(kdn_shp) +
  geom_sf(aes(fill = second.plus )) +
  scale_fill_viridis_c() +
  geom_sf_text(aes(label = NAME_2 ), size = 4, check_overlap = TRUE,color="white",fontface ="bold") +
  theme_void() +
  labs(fill = "Proportion",
       title = "Women who attended secondary or higher")

educ_map <- pAnyEdmap + pSecondmap + plot_layout(ncol = 2)
png(filename = paste0(outdir,"/plots/educ_map.png"),height=900,width=1200)
educ_map
dev.off()



##### D2: Travel time from Malaria Atlas project #####
## We want to calculate average or median travel time using the MAP data 
## cite: https://www.nature.com/articles/s41591-020-1059-1
## To aggregate to LGA level we need to appropriately weight by population 
## Using GRID3 population data to do this
kaduna_vect <- vect(kdn_shp)

# travel time to health facility by motor
tt_raster <- rast("raw data/travel time/202001_Global_Motorized_Travel_Time_to_Healthcare_NGA.tiff")
kdn_shp$tt_mean_unweighted <- exact_extract(tt_raster, kdn_shp, 'mean')

### NEEDS FOLLOW_UP: NEED TO FIGURE OUT HOW TO COMPUTE THE WEIGHTED TRAVEL TIME GIVEN THE SPARSITY OF POPULATION DATA

# population data 
# pop_raster <- rast("raw data/NGA_population_v2_1_gridded.tif")
# 
# #GRID3 data is 100m x 100m so need to aggregate to 1km x 1km to match the MAP data
# pop_1km <- aggregate(pop_raster, fact = 10, fun = sum)
# summary(values(pop_1km))
# 
# #check to see all rasters are the same CRS
# crs(pop_kdn); crs(pop_1km); st_crs(kdn_shp)
# 
# ## Align the travel time and population data 
# ## We need to make sure the grids between the two raster datasets align
# pop_1km_aligned <- resample(pop_1km, tt_raster, method = "near")
# summary(values(pop_1km_aligned))
# 
# ## some checks on coordinates
# crs(tt_raster); crs(pop_1km_aligned); st_crs(kdn_shp) #should all be same coordinate system
# 
# plot(tt_raster); plot(st_geometry(kdn_shp), add = TRUE, col = "red")
# plot(pop_1km_aligned); plot(st_geometry(kdn_shp), add = TRUE, col = "red")
# 
# 
# ## Stack rasters 
# stacked <- c(tt_raster, pop_1km_aligned)
# names(stacked) <- c("travel_time", "population")
# 
# ## computed weighted MEAN travel times
# # this function computes the weighted mean trvel time for LGA by LGA
# tt_mean <- exact_extract(stacked, kdn_shp[1,], function(df) {
#   weighted.mean(df$travel_time, df$population, na.rm = TRUE)
# }, summarize_df = TRUE)


pTTmap <- ggplot(kdn_shp) +
  geom_sf(aes(fill = tt_mean_unweighted )) +
  scale_fill_viridis_c() +
  geom_sf_text(aes(label = NAME_2 ), size = 4, check_overlap = TRUE,color="white",fontface ="bold") +
  theme_void() +
  labs(fill = "Travel Time (min)",
       title = "Unweighted travel time to nearest facility")
png(filename = paste0(outdir,"/plots/traveltime_map.png"),height=900,width=1200)
pTTmap
dev.off()



##### D3: ANC volumes #####
## HMIS data sent by Hadiza
hmis_full <- readxl::read_excel("raw data/Modeldata-IDMCHAI_02032025.xls")
### NEEDS FOLLOW_UP: Get updated data to match inclusion period for MMR

names(hmis_full) <- gsub(pattern = " ",replacement = ".",x=names(hmis_full))
anc <- hmis_full[,c("organisationunitname",grep(pattern = "ANC",x = names(hmis_full),value = TRUE))]

## Fixing LGA names
anc$LGA <- anc$organisationunitname
anc$LGA <- trimws(gsub("kd | Local Government Area","",anc$LGA))
anc$LGA <- gsub(" ",replacement = "_",anc$LGA)
anc$LGA <- gsub("'",replacement = "",anc$LGA)

## Converting from wide to long for easy average
anc_long <- anc %>%
  pivot_longer(
    cols = matches("^ANC\\."), 
    names_to = c("visit_type", "month"),
    names_pattern = "^(ANC\\.[^\\.]+)\\.*[Vv]isit\\.*(.*)$",
    values_to = "value"
  ) %>%
  mutate(
    # Normalize and parse month names
    month = str_replace_all(month, "\\.+", "."),
    date = parse_date_time(month, orders = "my", locale = "en_US")
  )



pANCvis <- ggplot(anc_long, aes(x = date, y = value, color = visit_type)) +
              geom_line(na.rm = TRUE) +
               geom_point(na.rm = TRUE) +
  
              facet_wrap(~ LGA, scales = "free_y") +
              labs(
                title = "ANC Visit Trends by LGA",
                x = "Date",
                y = "Number of Visits",
                color = "Visit Type"
              ) +
              theme_minimal() +
              theme(
                axis.text.x = element_text(angle = 45, hjust = 1),
                legend.position = "bottom"
              )
png(filename = paste0(outdir,"/plots/ANCvol_map.png"),height=900,width=1200)
pANCvis
dev.off()


## Calculate average over the observation period
ancAvg <- anc_long %>%
  group_by(visit_type, LGA) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

## convert to wide
ancAvg <- ancAvg %>%
  pivot_wider(
    names_from = visit_type,
    values_from = value
  )



#### E: Putting all together in one dataset for modeling input ####
pred1 <- st_drop_geometry(kdn_shp[,c("NAME_2","ever_school","second.plus","tt_mean_unweighted")])
cbhmis_final <- left_join(cbhmis_agg,pred1,by=c("LGA"="NAME_2"))
cbhmis_final <- left_join(cbhmis_final,ancAvg,by=c("LGA"))

write.csv(cbhmis_final,"cbhmis_data_for_model.csv",row.names = F)                           
