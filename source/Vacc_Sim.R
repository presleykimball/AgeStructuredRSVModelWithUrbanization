# This function simulates vaccination and returns lists of data
Vacc_Sim <- function(parms,tau_1,tau_2,time_of_vacc,time_post_vacc, # put in simulation parameters
                     chis, # put in base range of coverage rates
                     mat_vacc = TRUE, baby_vacc = TRUE, # optional: Turn respective vaccination on or off
                     vacc_ratio = 1 # optional: leverage ratio of maternal:baby vaccine to change coverage rates
                     ){
  # Set burnout period of 20 years+ 5 years want to compare to without vaccine + time_post_vacc
  times <- seq(0,(ceiling(365.25/7*(20+time_of_vacc+time_post_vacc))),by = 0.5)
  # set tau and time_of_vacc as given
  parms$tau_1 <- tau_1
  parms$tau_2 <- tau_2
  parms$time_of_vacc <- ceiling(365.25/7*(20+time_of_vacc))
  
  # Create time vectors for the ode solver
  post_burn_week <- ceiling(365.25/7*(20))
  post_burn_times <- times[times>post_burn_week]
  
  # Create matrices for later storage
  # under1 stores under 1 cases at each time against each chi
  under1 <- matrix(nrow = length(post_burn_times), ncol = (length(chis)+1))
  colnames(under1) <- c("time",chis)
  under1[,1] <- post_burn_times
  # total_stores total hospitalizations for each time against each chi
  total_hosp <- matrix(nrow = length(post_burn_times), ncol = (length(chis)+1))
  colnames(total_hosp) <- c("time",chis)
  total_hosp[,1] <- post_burn_times
  # prop_coverage stores the number of hospitalizations for each age against each chi
  prop_coverage <- matrix(nrow = length(chis),ncol = length(age_cats))
  row.names(prop_coverage) <- chis
  colnames(prop_coverage) <- age_cats
  # max_peak stores the max peak of each year against the coverage rate for 
  max_peak = list()
  # S_list stores susceptible time series
  S_list <- data.frame(post_burn_times)
  I_list <- data.frame(post_burn_times)
    if (mat_vacc == FALSE){
      mat_chis <- chis*0
    } else{
      mat_chis <- chis
    }
    if (baby_vacc == FALSE){
      baby_chis <- chis*0
    } else{
      baby_chis <- chis
    }
  if (vacc_ratio !=1 & mat_vacc == TRUE & baby_vacc == TRUE){
    if (vacc_ratio<1){ # more baby vaccine
      mat_chis <- mat_chis*vacc_ratio
    } else { # more maternal vaccine
      baby_chis <- baby_chis*(1/vacc_ratio)
    }
  }
  # begin for loop to solve over the chi vectors
  for(i in 1:length(chis)){
    # set the chi values
    parms$chi_b <- baby_chis[i]
    parms$chi_m <- mat_chis[i]
    
    
    # solve the ODE over given times
    solution <- ode(y=States_IC_vec, t=times,func=SIR, method = "ode45",
                    parms=parms)
    
    # find the hospitalizations
    hospitalizations <- Track_Hospitalizations(solution,parms)
    
    # make solution and hosp tracker for after burn period
    solution_post_burn <- solution[solution[,1]>ceiling(365.25/7*(20)),]
    hosp_post_burn <- hospitalizations[hospitalizations[,1]>ceiling(365.25/7*(20)),]
    S_vals <- solution_post_burn[,grep("S", colnames(solution))]
    S_vals <- rowSums(S_vals)
    S_list <- cbind(S_list, S_vals)
    I_vals <- solution_post_burn[,grep("I", colnames(solution))]
    I_vals <- rowSums(I_vals)
    I_list <- cbind(I_list, I_vals)
    
    under1[,i+1] <- hosp_post_burn[,2]
    total_hosp[,i+1] <- rowSums(hosp_post_burn[,2:ncol(hosp_post_burn)])
    
    hosp_post_vacc <- hosp_post_burn[hosp_post_burn[,1]>ceiling(365.25/7*(20+time_of_vacc+1)),]
    prop_coverage[i,] <- colSums(hosp_post_vacc[,2:ncol(hosp_post_vacc)])
    
    hosp_post_vacc <- as.data.frame(hosp_post_vacc)
    hosp_post_vacc <- hosp_post_vacc %>%
      mutate(year = as.factor(floor(time*7/365.25)-20-time_of_vacc+1))
    mat <- hosp_post_vacc %>%
      group_by(year) %>%
      summarize(across(matches("-"), max, .names = "{col}"), .groups = "drop")
    mat = mat[1:(nrow(mat)-1),]
    max_peak <- append(max_peak, list(mat))
  }
  under1[,1] <- under1[,1]- min(under1[,1])
  total_hosp[,1] <- total_hosp[,1]-min(total_hosp[,1])
  return(list(under1,total_hosp,prop_coverage,max_peak, S_list, I_list))
}