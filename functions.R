#################################################################
##   LGA estimates of community MMR                            ##
##   functions.R                                               ##
##   Purpose:  Assorted functions for plotting, model etc.     ##
##   Author: Anu Mishra                                        ##  
##   Created 6/9/25                                            ##
#################################################################



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
