#################################################################
##   LGA estimates of community MMR                            ##
##   functions.R                                               ##
##   Purpose:  Assorted functions for plotting, model etc.     ##
##   Author: Anu Mishra                                        ##  
##   Created 6/9/25                                            ##
#################################################################


summarize_matrix <- function(mat) {
  apply(mat, 2, function(x) {
    c(mean = mean(x), median = quantile(x, 0.5), sd = sd(x),
      l95 = quantile(x, 0.025),
      u95 = quantile(x, 0.975))
  }) |> t() |> as.data.frame()
}