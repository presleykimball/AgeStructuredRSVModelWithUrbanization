# Auxiliary function used in Error_Bar_Plots. Calculates means and confidence intervals.
mean_sd_df <- function(data_frm,top_num,conf_lvl){
  data_frm <- data_frm %>%
    mutate(mean = apply(data_frm[ , 2:2:(top_num+1)], 1, mean)) %>%
    mutate(std_dev = apply(data_frm[ , 2:(top_num+1)], 1, sd)) %>%
    mutate(upper_bd = mean+std_dev*qnorm(1-(1-conf_lvl)/2)) %>%
    mutate(lower_bd = mean+std_dev*qnorm((1-conf_lvl)/2)) 
  return(data_frm)
}