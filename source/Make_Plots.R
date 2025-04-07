# This function makes the basic plots without error bars
Make_Plots <- function(param_set,pop_size,data_sum,data_age_cat, State_name){
  # Clean param data
  param_set <- param_set[order(param_set$Loss, decreasing = FALSE),]
  parms <- as.list(param_set[1,]) # take only the best parameter set
  print(parms) 
  
  # Get solution using best params
  solution <- ODE_solve_SIR(parms,States_IC_vec,20,522)
  
  # Get hospitalization data from solution
  hospitalizations <- Track_Hospitalizations(solution = solution,parms = parms)
  hospitalizations[,2:ncol(hospitalizations)] <- hospitalizations[,2:ncol(hospitalizations)]*pop_size
  
  # Calculate total population
  N_tot <- rowSums(solution[,2:ncol(solution)])
  
  # Rescale timescales to start weeks at 0
  solution[,1] <- solution[,1]- min(solution[,1])
  hospitalizations[,1]<- hospitalizations[,1]-min(hospitalizations[,1])
  
  # Make sure the data is correctly ordered for plotting
  data_sum <- data_sum[order(data_sum$AWEEK1),]
  
  # Plot total population
  plot(solution[,1],N_tot,type = 'l', xlab = "Weeks",ylab = "Total Population")
  
  
  # Plot hospitalizations against data
  plot(data_sum$AWEEK1,data_sum$cases,type = 'l', 
       main = paste("Total Hospitalizations in", State_name),xlab = "Weeks",ylab = "Count", col = 'black',lwd =2,
       bty = 'l', ylim = c(0,max(c(rowSums(hospitalizations[,2:ncol(hospitalizations)]),data_sum$cases)))
  )
  lines(hospitalizations[,1],rowSums(hospitalizations[,2:ncol(hospitalizations)]),col=palette_4[1],lwd =2)
  legend("bottomright",legend = c("data","predicted"),col = c('black',palette_4[1]),lty=c(1,1),cex=0.8,bg = 'white')
  
  
  # print the minimum cases
  print(paste("Minimum achieved:",min(rowSums(hospitalizations[,2:ncol(hospitalizations)]))))
  
  # plot for first six age categories
  for (i in 2:7) {
    sub <- data_age_cat[data_age_cat$age_cat==age_cats[i-1],]
    plot(sub$AWEEK1,sub$cases,col = 'black',type = 'l',
         main =paste("Total Hospitalizations for",colnames(hospitalizations)[i], 
                     "year olds in", State_name), xlab = "Weeks", ylab = "Count",lwd=2)
    lines(hospitalizations[,1],hospitalizations[,i],,col= palette_4[1],lwd=2)
    legend("bottomright",legend = c("data","predicted"),col = c("black",palette_4[1]),lty=c(1,1),cex=0.8)
  }
  
  # Modify simulated data for plotting
  age_infected <- colSums(hospitalizations[,2:ncol(hospitalizations)]) # get total cases for each age group
  age_infected <- age_infected/sum(age_infected) # calculate as proportion
  age_infected <- as.data.frame(age_infected) # correct format
  colnames(age_infected) <- "predicted" # set name for plotting
  age_infected$age_cat <- age_cats # add in age cat labels for plotting
  
  # Modify actual data for plotting
  data_prop <- data_age_cat %>%
    group_by(age_cat)%>%
    summarise(data = sum(cases))
  data_prop$data <- data_prop$data/sum(data_prop$data) # calculate as proportion
  # Merge data with simulated data
  comparison <- merge(data_prop,age_infected,by= "age_cat")
  # fix format so that we can plot
  comparison <- comparison[match(age_cats,comparison$age_cat),]
  row.names(comparison) <- comparison[,1]
  comparison <- comparison[,-1]
  comparison <- as.matrix(comparison)
  comparison<- apply(t(comparison),2,rev)
  comparison <- comparison[c("data","predicted"),]
  row.names(comparison) <- c("Data","Predicted")
  
  # create plot
  barplot(comparison,beside = TRUE,legend.text = TRUE,
          ylim = c(0,0.8),
          main = paste("Proportion of Reported RSV by Age Group in", State_name,"Best Case"),
          col = c(palette_4[1],palette_4[2]),
          xlab = "Age Category",
          ylab = "Proportion of Cases",
          las =2)
  
  # # Create plot of only under 10
  comparison_young <- comparison[,colnames(comparison) %in% c("0-1","1-2","2-3","3-4",
                                                              "4-5","5-10")]
  
  barplot(comparison_young,beside = TRUE,legend.text = TRUE,
          ylim = c(0,0.8),
          main = paste("Proportion of Reported RSV by Age Group in", State_name,"Best Case"),
          col = c(palette_4[1],palette_4[2]),
          xlab = "Age Category",
          ylab = "Proportion of Cases")
  
}