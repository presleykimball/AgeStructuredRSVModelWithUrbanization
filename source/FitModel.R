# This function calculates the loss for given parameter sets and data
fitmodel <- function(parms,data,pop_size){
  burn_pd = 20 # burn period in years
  
  # Get simulated solution
  solution <- ODE_solve_SIR(parms,States_IC_vec,burn_pd,522)
  
  # Get simulated hospitalizations
  hospitalizations <- Track_Hospitalizations(solution = solution,parms = parms)
  
  # Re-scale hospitalizations to match the data
  hospitalizations[,2:ncol(hospitalizations)] <- hospitalizations[,2:ncol(hospitalizations)]*pop_size
  
  # Format data to look like hospitalizations data
  data <- data %>%
    pivot_wider(names_from = age_cat,values_from = cases)
  data <- data[,c("AWEEK1",age_cats)]
  data <- as.matrix(data)
  # set all NA data values to 0
  data[is.na(data)] <- 0
  # calculate loss as Frobenius norm
  Loss <- norm(hospitalizations[,2:ncol(hospitalizations)] - data[,2:ncol(data)], type = "F")
  return(list(Loss))
}