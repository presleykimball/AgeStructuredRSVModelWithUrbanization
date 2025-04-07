Urban_Rural_ErrorBarPlots <- function(parms,urban_pop_size,rural_pop_size,
                                      top_num,conf_lvl,
                                      urban_age_cat,rural_age_cat,
                                      urban_tau,rural_tau,State_name){
  
  parms <- parms[order(parms$Loss, decreasing = FALSE),]
  parms <- parms[1:top_num,]
  parms$tau_1 <- 0
  
  # set burn period
  burn_pd = 20
  
  # create data frames for storage
  
  
  urban_prop <- as.data.frame(age_cats)
  rural_prop <- as.data.frame(age_cats)
  
  
  for(i in 1:top_num){
    # Get urban simulated solution and hospitalizations
    parms$tau_2 <- urban_tau
    urban_solution <- ODE_solve_SIR(parms[i,],States_IC_vec,burn_pd,522)
    urban_hospitalizations <- Track_Hospitalizations(urban_solution,parms[i,])
    urban_hospitalizations[,2:ncol(urban_hospitalizations)] <-
      urban_hospitalizations[,2:ncol(urban_hospitalizations)]*urban_pop_size
    
    # Get rural simulated solution and hospitalizations
    parms$tau_2 <- rural_tau
    rural_solution <- ODE_solve_SIR(parms[i,],States_IC_vec,burn_pd,522)
    rural_hospitalizations <- Track_Hospitalizations(rural_solution,parms[i,])
    rural_hospitalizations[,2:ncol(rural_hospitalizations)] <-
      rural_hospitalizations[,2:ncol(rural_hospitalizations)]*rural_pop_size
    
    # calculate RSV fraction
    urban_frac <- colSums(urban_hospitalizations[,2:ncol(urban_hospitalizations)])
    urban_frac <- urban_frac/sum(urban_frac)
    
    rural_frac <- colSums(rural_hospitalizations[,2:ncol(rural_hospitalizations)])
    rural_frac <- rural_frac/sum(rural_frac)
    
    # store the appropriate hospitalizations in each data frame
    urban_prop[[paste0("Iteration_", i)]] <- urban_frac
    rural_prop[[paste0("Iteration_", i)]] <- rural_frac
  }
  
  # calculate means and lower and upper bounds of 95% confidence interval
  urban_prop<- mean_sd_df(urban_prop,top_num,conf_lvl)
  rural_prop <- mean_sd_df(rural_prop,top_num,conf_lvl)
  
  # calculate error bar limits for each age
  urban_prop <- urban_prop %>%
    mutate(upper_bd = mean+std_dev*qnorm(1-(1-conf_lvl)/2)) %>%
    mutate(lower_bd = mean+std_dev*qnorm((1-conf_lvl)/2)) 
  rural_prop <- rural_prop %>%
    mutate(upper_bd = mean+std_dev*qnorm(1-(1-conf_lvl)/2)) %>%
    mutate(lower_bd = mean+std_dev*qnorm((1-conf_lvl)/2)) 
  
  # calculate the fractions of RSV per age cat for the data
  urban_data_prop <- urban_age_cat %>%
    group_by(age_cat)%>%
    summarise(data = sum(cases))
  urban_data_prop$data <- urban_data_prop$data/sum(urban_data_prop$data)
  
  rural_data_prop <- rural_age_cat %>%
    group_by(age_cat)%>%
    summarise(data = sum(cases))
  rural_data_prop$data <- rural_data_prop$data/sum(rural_data_prop$data)
  
  # select only under 10 data
  urban_prop <- urban_prop[urban_prop$age_cat %in% c("0-1","1-2","2-3","3-4","4-5","5-10"),]
  rural_prop <- rural_prop[rural_prop$age_cat %in% c("0-1","1-2","2-3","3-4","4-5","5-10"),]
  urban_data_prop <- urban_data_prop[urban_data_prop$age_cat %in%
                                       c("0-1","1-2","2-3","3-4","4-5","5-10"),] 
  rural_data_prop <- rural_data_prop[rural_data_prop$age_cat %in%
                                       c("0-1","1-2","2-3","3-4","4-5","5-10"),] 
  
  #select only the columns needed from our predicted
  urban_prop <- urban_prop[,c("age_cats","mean","upper_bd","lower_bd")]
  rural_prop <- rural_prop[,c("age_cats","mean","upper_bd","lower_bd")]
  
  # start creating data frames for bar plot
  comparison <- merge(urban_prop,rural_prop,by= "age_cats")
  colnames(comparison) <- c("age_cat","urban_mean","urban_upper_bd","urban_lower_bd",
                            "rural_mean","rural_upper_bd","rural_lower_bd")
  data_comparison <- merge(urban_data_prop,rural_data_prop,by = "age_cat")
  colnames(data_comparison) <- c("age_cat","urban_data","rural_data")
  data_comparison <- merge(comparison,data_comparison,by="age_cat")
  
  row.names(comparison) <- comparison[,1]
  row.names(data_comparison) <- data_comparison[,1]
  comparison <- comparison[,-1]
  data_comparison <- data_comparison[,-1]
  
  main_comparison <- comparison[,c("urban_mean","rural_mean")]
  main_data_comparison <- data_comparison[,c("urban_mean","rural_mean","urban_data","rural_data")]
  
  main_comparison <- as.matrix(main_comparison)
  main_data_comparison <- as.matrix(main_data_comparison)
  main_comparison<- apply(t(main_comparison),2,rev)
  main_data_comparison<- apply(t(main_data_comparison),2,rev)
  
  main_comparison <- main_comparison[c("urban_mean","rural_mean"),]
  row.names(main_comparison) <- c("Urban (mean)","Rural (mean)") 
  
  main_data_comparison <- main_data_comparison[c("urban_data","urban_mean",
                                                 "rural_data","rural_mean"),]
  row.names(main_data_comparison) <- c("Urban Data","Predicted Urban (mean)",
                                       "Rural Data","Predicted Rural (mean)") 
  
  base_bar <- barplot(main_comparison,beside = TRUE,legend.text = TRUE,
                      ylim = c(0,ceiling(max(comparison)*10)/10),
                      main = paste("Proportion of Reported RSV by Age Group in", State_name),
                      col = c(palette_4[1],palette_4[2]),
                      xlab = "Age Category",
                      ylab = "Proportion of Cases")
  arrows(x0 = base_bar[1, ], 
         y0 = comparison$urban_lower_bd, 
         x1 = base_bar[1,], 
         y1 = comparison$urban_upper_bd, 
         angle = 90, 
         code = 3, 
         length = 0.1, 
         col = 'black', 
         lwd = 3)
  arrows(x0 = base_bar[2,], 
         y0 = comparison$rural_lower_bd, 
         x1 = base_bar[2,], 
         y1 = comparison$rural_upper_bd, 
         angle = 90, 
         code = 3, 
         length = 0.1, 
         col = 'black', 
         lwd = 3)
  
  base_data_bar <- barplot(main_data_comparison,beside = TRUE,legend.text = TRUE,
                           ylim = c(0,ceiling(max(comparison)*10)/10),
                           main = paste("Proportion of Reported RSV by Age Group in", State_name),
                           col = c(palette_4[1],palette_4[2],palette_4[3],palette_4[4]),
                           xlab = "Age Category",
                           ylab = "Proportion of Cases")
  arrows(x0 = base_data_bar[2, ], 
         y0 = comparison$urban_lower_bd, 
         x1 = base_data_bar[2,], 
         y1 = comparison$urban_upper_bd, 
         angle = 90, 
         code = 3, 
         length = 0.1, 
         col = 'black', 
         lwd = 3)
  arrows(x0 = base_data_bar[4,], 
         y0 = comparison$rural_lower_bd, 
         x1 = base_data_bar[4,], 
         y1 = comparison$rural_upper_bd, 
         angle = 90, 
         code = 3, 
         length = 0.1, 
         col = 'black', 
         lwd = 3)
}