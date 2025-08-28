# This function creates same plots as Make_Plots.R except with error bars
ErrorBarPlots <- function(parms,pop_size,top_num,conf_lvl,data_sum,data_age_cat,State_name){
  
  parms <- parms[order(parms$Loss, decreasing = FALSE),]
  parms <- parms[1:top_num,]
  
  # set burn period
  burn_pd = 20
  
  # create data frames for storage 
  time <- 0:521
  total_hosp <- as.data.frame(time)
  hosp_0_1 <- as.data.frame(time)
  hosp_1_2 <- as.data.frame(time)
  hosp_2_3 <- as.data.frame(time)
  hosp_3_4 <- as.data.frame(time)
  hosp_4_5 <- as.data.frame(time)
  hosp_5_10 <- as.data.frame(time)
  for(i in 1:top_num){
    # Get simulated solution
    solution <- ODE_solve_SIR(parms[i,],States_IC_vec,burn_pd,522)
    
    # Get simulated hospitalizations
    hospitalizations <- Track_Hospitalizations(solution,parms[i,])
    
    # Reset hospitalizations time-scale to match the timescale of the data
    hospitalizations[,"time"] <- hospitalizations[,"time"]- min(hospitalizations[,"time"])
    
    # Re-scale hospitalizations to match the data
    hospitalizations[,2:ncol(hospitalizations)] <- hospitalizations[,2:ncol(hospitalizations)]*pop_size
    
    # store the appropriate hospitalizations in each data frame
    total_hosp[[paste0("Iteration_", i)]] <- rowSums(hospitalizations[,2:ncol(hospitalizations)])
    hosp_0_1[[paste0("Iteration_", i)]] <- hospitalizations[,2]
    hosp_1_2[[paste0("Iteration_", i)]] <- hospitalizations[,3]
    hosp_2_3[[paste0("Iteration_", i)]] <- hospitalizations[,4]
    hosp_3_4[[paste0("Iteration_", i)]] <- hospitalizations[,5]
    hosp_4_5[[paste0("Iteration_", i)]] <- hospitalizations[,6]
    hosp_5_10[[paste0("Iteration_", i)]] <- hospitalizations[,7]
  }
  
  # calculate means and lower and upper bounds of 95% confidence interval
  total_hosp <- mean_sd_df(total_hosp,top_num,conf_lvl)
  hosp_0_1 <- mean_sd_df(hosp_0_1,top_num,conf_lvl)
  hosp_1_2 <- mean_sd_df(hosp_1_2,top_num,conf_lvl)
  hosp_2_3 <- mean_sd_df(hosp_2_3,top_num,conf_lvl)
  hosp_3_4 <- mean_sd_df(hosp_3_4,top_num,conf_lvl)
  hosp_4_5 <- mean_sd_df(hosp_4_5,top_num,conf_lvl)
  hosp_5_10 <- mean_sd_df(hosp_5_10,top_num,conf_lvl)
  
  data_sum <- data_sum[order(data_sum$AWEEK1,decreasing = FALSE),]
  
  #### Plot hospitalizations against data for the 6 youngest age categories
  pdf("figures/Errorbars_totalhosp.pdf")
    plot(data_sum$AWEEK1,data_sum$cases,type = 'l', 
         main = paste("Total Hospitalizations in", State_name),xlab = "Weeks",
         ylab = "Count",lwd = 2,col =palette_4[1])
    lines(total_hosp$time,total_hosp$upper_bd,col= palette_4[4],lty = 2,lwd = 2)
    lines(total_hosp$time,total_hosp$lower_bd,col= palette_4[4],lty = 2,lwd = 2)
    lines(total_hosp$time,total_hosp$mean,col=palette_4[3],lwd = 2)
    lines(total_hosp$time,total_hosp[,2],col=palette_4[2],lty = 2,lwd = 2)
    legend("bottomright",legend = c("data","best case","mean predicted","95% CI"),
           col = c(palette_4[1],palette_4[2],palette_4[3],palette_4[4]),
           lty=c(1,2,1,2),cex=0.8,lwd = c(2,2,2,2))
  dev.off()
  plot(data_sum$AWEEK1,data_sum$cases,type = 'l', 
       main = paste("Total Hospitalizations in", State_name),xlab = "Weeks",
       ylab = "Count",lwd = 2,col =palette_4[1])
  lines(total_hosp$time,total_hosp$upper_bd,col= palette_4[4],lty = 2,lwd = 2)
  lines(total_hosp$time,total_hosp$lower_bd,col= palette_4[4],lty = 2,lwd = 2)
  lines(total_hosp$time,total_hosp$mean,col=palette_4[3],lwd = 2)
  lines(total_hosp$time,total_hosp[,2],col=palette_4[2],lty = 2,lwd = 2)
  legend("bottomright",legend = c("data","best case","mean predicted","95% CI"),
         col = c(palette_4[1],palette_4[2],palette_4[3],palette_4[4]),
         lty=c(1,2,1,2),cex=0.8,lwd = c(2,2,2,2))
  
  sub <- data_age_cat[data_age_cat$age_cat=="0-1",]
  pdf("figures/Errorbars_0_1.pdf")
    plot(sub$AWEEK1,sub$cases,col = palette_4[1],type = 'l',
         main =paste("Total Hospitalizations for 0-1 year olds in", State_name), 
         xlab = "Weeks", ylab = "Count",lwd = 2)
    lines(hosp_0_1$time,hosp_0_1$upper_bd,col= palette_4[4],lty = 2,lwd = 2)
    lines(hosp_0_1$time,hosp_0_1$lower_bd,col= palette_4[4],lty = 2,lwd = 2)
    lines(hosp_0_1$time,hosp_0_1$mean,col=palette_4[3],lwd = 2)
    lines(hosp_0_1$time,hosp_0_1[,2],col=palette_4[2],lty = 2,lwd = 2)
    legend("bottomright",legend = c("data","best case","mean predicted","95% CI"),
           col = c(palette_4[1],palette_4[2],palette_4[3],palette_4[4]),
           lty=c(1,2,1,2),cex=0.8,lwd = c(2,2,2,2))
  dev.off()
  plot(sub$AWEEK1,sub$cases,col = palette_4[1],type = 'l',
       main =paste("Total Hospitalizations for 0-1 year olds in", State_name), 
       xlab = "Weeks", ylab = "Count",lwd = 2)
  lines(hosp_0_1$time,hosp_0_1$upper_bd,col= palette_4[4],lty = 2,lwd = 2)
  lines(hosp_0_1$time,hosp_0_1$lower_bd,col= palette_4[4],lty = 2,lwd = 2)
  lines(hosp_0_1$time,hosp_0_1$mean,col=palette_4[3],lwd = 2)
  lines(hosp_0_1$time,hosp_0_1[,2],col=palette_4[2],lty = 2,lwd = 2)
  legend("bottomright",legend = c("data","best case","mean predicted","95% CI"),
         col = c(palette_4[1],palette_4[2],palette_4[3],palette_4[4]),
         lty=c(1,2,1,2),cex=0.8,lwd = c(2,2,2,2))
  
  sub <- data_age_cat[data_age_cat$age_cat=="1-2",]
  pdf("figures/Errorbars_1_2.pdf")
    plot(sub$AWEEK1,sub$cases,col = palette_4[1],type = 'l',
         main =paste("Total Hospitalizations for 1-2 year olds in", State_name), 
         xlab = "Weeks", ylab = "Count",lwd = 2)
    lines(hosp_1_2$time,hosp_1_2$upper_bd,col= palette_4[4],lty = 2,lwd = 2)
    lines(hosp_1_2$time,hosp_1_2$lower_bd,col= palette_4[4],lty = 2,lwd = 2)
    lines(hosp_1_2$time,hosp_1_2$mean,col=palette_4[3],lwd = 2)
    lines(hosp_1_2$time,hosp_1_2[,2],col=palette_4[2],lty = 2,lwd = 2)
    legend("bottomright",legend = c("data","best case","mean predicted","95% CI"),
           col = c(palette_4[1],palette_4[2],palette_4[3],palette_4[4]),
           lty=c(1,2,1,2),cex=0.8,lwd = c(2,2,2,2))
  dev.off()
  plot(sub$AWEEK1,sub$cases,col = palette_4[1],type = 'l',
       main =paste("Total Hospitalizations for 1-2 year olds in", State_name), 
       xlab = "Weeks", ylab = "Count",lwd = 2)
  lines(hosp_1_2$time,hosp_1_2$upper_bd,col= palette_4[4],lty = 2,lwd = 2)
  lines(hosp_1_2$time,hosp_1_2$lower_bd,col= palette_4[4],lty = 2,lwd = 2)
  lines(hosp_1_2$time,hosp_1_2$mean,col=palette_4[3],lwd = 2)
  lines(hosp_1_2$time,hosp_1_2[,2],col=palette_4[2],lty = 2,lwd = 2)
  legend("bottomright",legend = c("data","best case","mean predicted","95% CI"),
         col = c(palette_4[1],palette_4[2],palette_4[3],palette_4[4]),
         lty=c(1,2,1,2),cex=0.8,lwd = c(2,2,2,2))
  
  sub <- data_age_cat[data_age_cat$age_cat=="2-3",]
  pdf("figures/Errorbars_2_3.pdf")
    plot(sub$AWEEK1,sub$cases,col = palette_4[1],type = 'l',
         main =paste("Total Hospitalizations for 2-3 year olds in", State_name), 
         xlab = "Weeks", ylab = "Count")
    lines(hosp_2_3$time,hosp_2_3$upper_bd,col= palette_4[4],lty = 2,lwd = 2)
    lines(hosp_2_3$time,hosp_2_3$lower_bd,col= palette_4[4],lty = 2,lwd = 2)
    lines(hosp_2_3$time,hosp_2_3$mean,col=palette_4[3],lwd = 2)
    lines(hosp_2_3$time,hosp_2_3[,2],col=palette_4[2],lty = 2,lwd = 2)
    legend("bottomright",legend = c("data","best case","mean predicted","95% CI"),
           col = c(palette_4[1],palette_4[2],palette_4[3],palette_4[4]),
           lty=c(1,2,1,2),cex=0.8,lwd = c(2,2,2,2))
  dev.off()
  plot(sub$AWEEK1,sub$cases,col = palette_4[1],type = 'l',
       main =paste("Total Hospitalizations for 2-3 year olds in", State_name), 
       xlab = "Weeks", ylab = "Count")
  lines(hosp_2_3$time,hosp_2_3$upper_bd,col= palette_4[4],lty = 2,lwd = 2)
  lines(hosp_2_3$time,hosp_2_3$lower_bd,col= palette_4[4],lty = 2,lwd = 2)
  lines(hosp_2_3$time,hosp_2_3$mean,col=palette_4[3],lwd = 2)
  lines(hosp_2_3$time,hosp_2_3[,2],col=palette_4[2],lty = 2,lwd = 2)
  legend("bottomright",legend = c("data","best case","mean predicted","95% CI"),
         col = c(palette_4[1],palette_4[2],palette_4[3],palette_4[4]),
         lty=c(1,2,1,2),cex=0.8,lwd = c(2,2,2,2))
  
  sub <- data_age_cat[data_age_cat$age_cat=="3-4",]
  pdf("figures/Errorbars_3_4.pdf")
    plot(sub$AWEEK1,sub$cases,col = palette_4[1],type = 'l',
         main =paste("Total Hospitalizations for 3-4 year olds in", State_name), 
         xlab = "Weeks", ylab = "Count")
    lines(hosp_3_4$time,hosp_3_4$upper_bd,col= palette_4[4],lty = 2,lwd = 2)
    lines(hosp_3_4$time,hosp_3_4$lower_bd,col= palette_4[4],lty = 2,lwd = 2)
    lines(hosp_3_4$time,hosp_3_4$mean,col=palette_4[3],lwd = 2)
    lines(hosp_3_4$time,hosp_3_4[,2],col=palette_4[2],lty = 2,lwd = 2)
    legend("bottomright",legend = c("data","best case","mean predicted","95% CI"),
           col = c(palette_4[1],palette_4[2],palette_4[3],palette_4[4]),
           lty=c(1,2,1,2),cex=0.8,lwd = c(2,2,2,2))
  dev.off()
  plot(sub$AWEEK1,sub$cases,col = palette_4[1],type = 'l',
       main =paste("Total Hospitalizations for 3-4 year olds in", State_name), 
       xlab = "Weeks", ylab = "Count")
  lines(hosp_3_4$time,hosp_3_4$upper_bd,col= palette_4[4],lty = 2,lwd = 2)
  lines(hosp_3_4$time,hosp_3_4$lower_bd,col= palette_4[4],lty = 2,lwd = 2)
  lines(hosp_3_4$time,hosp_3_4$mean,col=palette_4[3],lwd = 2)
  lines(hosp_3_4$time,hosp_3_4[,2],col=palette_4[2],lty = 2,lwd = 2)
  legend("bottomright",legend = c("data","best case","mean predicted","95% CI"),
         col = c(palette_4[1],palette_4[2],palette_4[3],palette_4[4]),
         lty=c(1,2,1,2),cex=0.8,lwd = c(2,2,2,2))
  
  sub <- data_age_cat[data_age_cat$age_cat=="4-5",]
  pdf("figures/Errorbars_4_5.pdf")
    plot(sub$AWEEK1,sub$cases,col = palette_4[1],type = 'l',
         main =paste("Total Hospitalizations for 4-5 year olds in", State_name), 
         xlab = "Weeks", ylab = "Count")
    lines(hosp_4_5$time,hosp_4_5$upper_bd,col= palette_4[4],lty = 2,lwd = 2)
    lines(hosp_4_5$time,hosp_4_5$lower_bd,col= palette_4[4],lty = 2,lwd = 2)
    lines(hosp_4_5$time,hosp_4_5$mean,col=palette_4[3],lwd = 2)
    lines(hosp_4_5$time,hosp_4_5[,2],col=palette_4[2],lty = 2,lwd = 2)
    legend("bottomright",legend = c("data","best case","mean predicted","95% CI"),
           col = c(palette_4[1],palette_4[2],palette_4[3],palette_4[4]),
           lty=c(1,2,1,2),cex=0.8,lwd = c(2,2,2,2))
  dev.off()
  plot(sub$AWEEK1,sub$cases,col = palette_4[1],type = 'l',
       main =paste("Total Hospitalizations for 4-5 year olds in", State_name), 
       xlab = "Weeks", ylab = "Count")
  lines(hosp_4_5$time,hosp_4_5$upper_bd,col= palette_4[4],lty = 2,lwd = 2)
  lines(hosp_4_5$time,hosp_4_5$lower_bd,col= palette_4[4],lty = 2,lwd = 2)
  lines(hosp_4_5$time,hosp_4_5$mean,col=palette_4[3],lwd = 2)
  lines(hosp_4_5$time,hosp_4_5[,2],col=palette_4[2],lty = 2,lwd = 2)
  legend("bottomright",legend = c("data","best case","mean predicted","95% CI"),
         col = c(palette_4[1],palette_4[2],palette_4[3],palette_4[4]),
         lty=c(1,2,1,2),cex=0.8,lwd = c(2,2,2,2))
  
  sub <- data_age_cat[data_age_cat$age_cat=="5-10",]
  pdf("figures/Errorbars_4_5.pdf")
    plot(sub$AWEEK1,sub$cases,col = palette_4[1],type = 'l',
         main =paste("Total Hospitalizations for 5-10 year olds in", State_name), 
         xlab = "Weeks", ylab = "Count")
    lines(hosp_5_10$time,hosp_5_10$upper_bd,col= palette_4[4],lty = 2,lwd = 2)
    lines(hosp_5_10$time,hosp_5_10$lower_bd,col= palette_4[4],lty = 2,lwd = 2)
    lines(hosp_5_10$time,hosp_5_10$mean,col=palette_4[3],lwd = 2)
    lines(hosp_5_10$time,hosp_5_10[,2],col=palette_4[2],lty = 2,lwd = 2)
    legend("bottomright",legend = c("data","best case","mean predicted","95% CI"),
           col = c(palette_4[1],palette_4[2],palette_4[3],palette_4[4]),
           lty=c(1,2,1,2),cex=0.8,lwd = c(2,2,2,2))
  dev.off()
  plot(sub$AWEEK1,sub$cases,col = palette_4[1],type = 'l',
       main =paste("Total Hospitalizations for 5-10 year olds in", State_name), 
       xlab = "Weeks", ylab = "Count")
  lines(hosp_5_10$time,hosp_5_10$upper_bd,col= palette_4[4],lty = 2,lwd = 2)
  lines(hosp_5_10$time,hosp_5_10$lower_bd,col= palette_4[4],lty = 2,lwd = 2)
  lines(hosp_5_10$time,hosp_5_10$mean,col=palette_4[3],lwd = 2)
  lines(hosp_5_10$time,hosp_5_10[,2],col=palette_4[2],lty = 2,lwd = 2)
  legend("bottomright",legend = c("data","best case","mean predicted","95% CI"),
         col = c(palette_4[1],palette_4[2],palette_4[3],palette_4[4]),
         lty=c(1,2,1,2),cex=0.8,lwd = c(2,2,2,2))
  
  
  #### Creating bar plot
  # create data frame of counts by age cats
  hosp_age_cat <- as.data.frame(c("0-1","1-2","2-3","3-4","4-5","5-10"))
  colnames(hosp_age_cat) <-"age_cat"
  
  # calculate ratio of hospitalizations for each iteration among all subsets
  tot_hosp_iter <- colSums(total_hosp[,2:(top_num+1)])
  hosp_iter_0_1 <- colSums(hosp_0_1[,2:(top_num+1)])/tot_hosp_iter
  hosp_iter_1_2 <- colSums(hosp_1_2[,2:(top_num+1)])/tot_hosp_iter
  hosp_iter_2_3 <- colSums(hosp_2_3[,2:(top_num+1)])/tot_hosp_iter
  hosp_iter_3_4 <- colSums(hosp_3_4[,2:(top_num+1)])/tot_hosp_iter
  hosp_iter_4_5 <- colSums(hosp_4_5[,2:(top_num+1)])/tot_hosp_iter
  hosp_iter_5_10 <- colSums(hosp_5_10[,2:(top_num+1)])/tot_hosp_iter
  
  # calculate mean and standard deviation for each age
  hosp_age_cat$mean <- c(mean(hosp_iter_0_1),mean(hosp_iter_1_2),mean(hosp_iter_2_3),
                         mean(hosp_iter_3_4),mean(hosp_iter_4_5),mean(hosp_iter_5_10))
  hosp_age_cat$std_dev <- c(sd(hosp_iter_0_1),sd(hosp_iter_1_2),sd(hosp_iter_2_3),
                            sd(hosp_iter_3_4),sd(hosp_iter_4_5),sd(hosp_iter_5_10))
  # calculate error bar limits for each age
  hosp_age_cat <- hosp_age_cat %>%
    mutate(upper_bd = mean+std_dev*qnorm(1-(1-conf_lvl)/2)) %>%
    mutate(lower_bd = mean+std_dev*qnorm((1-conf_lvl)/2)) 
  
  # calculate the fractions of RSV per age cat for the data
  data_prop <- data_age_cat %>%
    group_by(age_cat)%>%
    summarise(data = sum(cases))
  data_prop$data <- data_prop$data/sum(data_prop$data)
  data_prop <- data_prop[data_prop$age_cat %in% c("0-1","1-2","2-3","3-4","4-5","5-10"),] #ensure correct order
  
  # start creating 
  comparison <- merge(data_prop,hosp_age_cat,by= "age_cat")
  row.names(comparison) <- comparison[,1]
  comparison <- comparison[,-1]
  main_comparison <- comparison[,1:2]
  main_comparison <- as.matrix(main_comparison)
  main_comparison<- apply(t(main_comparison),2,rev)
  main_comparison <- main_comparison[c("data","mean"),]
  row.names(main_comparison) <- c("Data","Predicted (Mean)") 
  
  # create bar plot
  pdf("figures/Errorbars_barplot.pdf")
  base_bar <- barplot(main_comparison,beside = TRUE,legend.text = TRUE,
                      ylim = c(0,ceiling(max(hosp_age_cat[,"upper_bd"])*10)/10),
                      main = paste("Proportion of Reported RSV by Age Group in", State_name),
                      col = c(palette_4[1],palette_4[2]),
                      xlab = "Age Category",
                      ylab = "Proportion of Cases")
  arrows(x0 = base_bar[2, ], 
         y0 = hosp_age_cat[,"lower_bd"] , 
         x1 = base_bar[2, ], 
         y1 = hosp_age_cat[,"upper_bd"], 
         angle = 90, 
         code = 3, 
         length = 0.1, 
         col = 'black',
         lwd = 3)
  dev.off()
}