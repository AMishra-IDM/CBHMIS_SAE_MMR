#### Plots for SAE Conference ###

rm(list=ls())

library(tidyverse)
library(ggplot2)
library(sf)


orig <- read.csv("C:/Users/anumi/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub/CBHMIS_SAE_MMR/model runs/model_2025-07-02/summary_stats.csv")
noANC <- read.csv("C:/Users/anumi/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub/CBHMIS_SAE_MMR/model runs/model_2025-07-02_no_ANCRef/summary_stats.csv")

data <- read.csv("cbhmis_data_for_model.csv")

orig <- left_join(orig,data[,c("LGA","repRate","ANC_ref","lb","deaths","MMR","anyEd",
                                            "tt_mean_unweighted","ANC.1st")])


### Kaduna Shape file ###
## shape file 
nga_shp <- st_read("raw data/gadm41_NGA_shp/gadm41_NGA_2.shp")
kdn_shp <- nga_shp[nga_shp$NAME_1=="Kaduna",]
kdn_shp$NAME_2 <- gsub(pattern = " ",replacement = "_", kdn_shp$NAME_2)
kdn_shp$NAME_2 <- gsub(pattern = "'",replacement = "", kdn_shp$NAME_2)


######## MMR vs interval width ########
orig$crWidth <- orig$MMR_upper_95 - orig$MMR_lower_95
map_out <- left_join(kdn_shp,orig,by=c("NAME_2"="LGA"))

p1 <- ggplot(orig, aes(x = MMR_median, y = crWidth)) +
  geom_point(aes(size = lb), color = "lightblue", alpha = 0.7) +
  scale_size_continuous(
    name = "Live Births",
    range = c(2, 12),
    breaks = c(2000, 4000, 6000, 8000, 10000),
    labels = paste0(seq(2, 10, 2), "k")) +
  labs(
    x = "Posterior Median MMR (per 100k LB)",
    y = "95% Credible Interval Width"
  ) + theme_bw() + 
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  ) + theme(legend.position = "bottom")


png(filename = "C:/Users/anumi/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub/CBHMIS_SAE_MMR/presentation plots/mean_vs_width.png",
    height=600,width=600)
p1
dev.off()

p2 <- ggplot(map_out) +
  geom_sf(aes(fill = MMR_median)) +
  scale_fill_gradientn(
    colors = c("#7FB77E", "#F7DC6F", "#E74C3C"),
    name = "MMR",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      label.position = "bottom",
      label.theme = element_text(angle = 45, size = 12)
    )
  ) +
  geom_sf_text(aes(label = NAME_2), size = 4, check_overlap = TRUE, color = "white", fontface = "bold") +
  theme_void() +
  labs(
    fill = "MMR (per 100k LB)",
    title = "Posterior median MMR"
  ) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.position = "bottom"
  )

png(filename = "C:/Users/anumi/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub/CBHMIS_SAE_MMR/presentation plots/MMRmap.png",
    height=600,width=600)
p2
dev.off()

summary(orig$MMR_median)



p3 <- ggplot(map_out) +
  geom_sf(aes(fill = pi_median )) +
  scale_fill_viridis_c() +
  geom_sf_text(aes(label = NAME_2 ), size = 4, check_overlap = TRUE,color="white",fontface ="bold") +
  theme_void() +
  labs(fill = "Proportion",
       title = "Posterior Median Reporting Probability") + 
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.position = "bottom"
  )

png(filename = "C:/Users/anumi/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub/CBHMIS_SAE_MMR/presentation plots/Pimap.png",
    height=600,width=600)
p3
dev.off()

p4 <- ggplot(orig, aes(x = reorder(factor(LGA), lambda_median))) +
  geom_point(aes(y = lambda_median, color = "Estimated (lambda)"),size=3) +
  geom_errorbar(aes(ymin = lambda_lower_95, ymax = lambda_upper_95), width = 0.5) +
  geom_point(aes(y = deaths, color = "Observed (z)"),size=3) +
  labs(x = "LGA", y = "Maternal Deaths",
       color = "Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c("Estimated (lambda)" = "black", "Observed (z)" = "red"),name = NULL) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  ) + theme(legend.position = "bottom")

png(filename = "C:/Users/anumi/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub/CBHMIS_SAE_MMR/presentation plots/deaths_obs_expected.png",
    height=600,width=600)
p4
dev.off()

range(orig$pi_median)


p5 <- ggplot(mapping = aes(x=scale(orig$ANC_ref,center = TRUE,scale = TRUE),y=orig$pi_median)) + 
          geom_point(size=3) + geom_smooth(method="lm",linetype=2) + 
        theme_bw() +
        xlab("ANC Referrals (scaled and centered)") + ylab("Posterior median reporting rate (pi)") +
theme(
  plot.title = element_text(size = 20, face = "bold"),
  axis.title = element_text(size = 18),
  axis.text = element_text(size = 16),
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 14)
) + theme(legend.position = "bottom")

png(filename = "C:/Users/anumi/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub/CBHMIS_SAE_MMR/presentation plots/anc_ref_coef.png",
    height=600,width=600)
p5
dev.off()


comp <- left_join(orig[,c("LGA","pi_median")],noANC[,c("LGA","pi_median")],by = "LGA",suffix = c(".orig",".noANC"))

long_pi <- comp %>%
  pivot_longer(
    cols = c(pi_median.orig, pi_median.noANC),
    names_to = "model",
    values_to = "pi"
  ) %>%
  mutate(
    model = recode(model,
                   "pi_median.orig" = "With ANC covariate",
                   "pi_median.noANC" = "Without ANC covariate"
    )
  )

p6 <- ggplot(long_pi, aes(x = model, y = pi)) +
  geom_boxplot() +
  labs(
    x = "",
    y = "Posterior Median Pi (Reporting Probability)",
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  ) + theme(legend.position = "none")

png(filename = "C:/Users/anumi/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub/CBHMIS_SAE_MMR/presentation plots/model_comp.png",
    height=400,width=600)
p6
dev.off()

