# This function uses the solution given to track the hospitalizations at each time step.
Track_Hospitalizations <- function(solution,parms){
  # Call in parameters
  alpha_1 <- parms[["alpha_1"]]
  alpha_2 <- parms[["alpha_2"]]
  alpha_3 <- parms[["alpha_3"]]
  q_1 <- parms[["q_1"]]
  q_2 <- parms[["q_2"]]
  q_3 <- parms[["q_3"]]
  q_4 <- parms[["q_4"]]
  q_5 <- parms[["q_5"]]
  q_6 <- parms[["q_6"]]
  
  # call in optional childcare amplification 
  tau_1 <- parms[["tau_1"]]
  tau_2 <- parms[["tau_2"]]
  if (!is.null(parms$Ktau_all) ){
    if( parms$Ktau_all == TRUE) { # if all mult is provided and on
      tau_vec_comp <- rep(tau_2,length(age_cats))} else {
        tau_vec_comp <- tau_vec(tau_1,tau_2)
      }
  } else{ # Ktau_all not provided or is False
    tau_vec_comp <- tau_vec(tau_1,tau_2)
  }
  
  # create q_vec to be the 6 fitted parameters followed by 1's
  q_vec <- c(q_1,q_2,q_3,q_4,q_5,q_6,rep(1,(length(h_1)-6)))
  
  # Initialize matrix
  Hospital <- matrix(NA,nrow = nrow(solution),ncol = length(age_cats)+1)
  colnames(Hospital) <- c("time",age_cats)
  
  # Store times
  Hospital[,"time"] <- solution[,"time"]
  
  # find hospitalizations for each time step
  for (t in solution[,"time"]){
    # collect solution at the given time step
    States <- solution[which(solution[, "time"] == t),]
    # take all data except for the time variable
    States <- States[2:length(States)]
    # convert vector data back to matrix form
    row_names <- sub(" .*", "", names(States))  
    col_names <- sub(".* ", "", names(States)) 
    States <- matrix(States, nrow = length(unique(row_names)), 
                     ncol = length(unique(col_names)), byrow = FALSE,
                     dimnames = list(unique(row_names), unique(col_names)))
    
    # Get susceptible and infected data as these are used to calculate transmission
    # and hospitalization
    S0 <-  States[,'S0']
    S1 <-  States[,'S1']
    S2 <-  States[,'S2']
    S3 <-  States[,'S3']
    I1 <-  States[,'I1']
    I2 <-  States[,'I2']
    I3 <-  States[,'I3']
    I4 <-  States[,'I4']
    
    # calculate transmission
    N_i <- rowSums(States) #number in each age category
    K_mod <- sweep(K,1,N_i,"/") # divide each row of contact matrix by number in that age_cat
    K_mod <- sweep(K_mod,2,tau_vec_comp,"*") # multiply each column by optional amplification tau
    K_mod <- sweep(K_mod,1,tau_vec_comp,"*") # multiply each row by optional amplification tau
    sums <- rho_vec[1]*(K_mod %*% I1)+rho_vec[2]*(K_mod %*% I2)+
      rho_vec[3]*(K_mod %*% I3)+rho_vec[4]*(K_mod %*% I4)
    b_vec <- alpha_1*(alpha_2*cos((2*pi*t-alpha_3)/365.25*7)+1)*sums
    
    # calculate hospitalizations
    Hospital[which(solution[, "time"] == t),2:ncol(Hospital)] <-
      b_vec*(sigma_vec[1]*S0*h_1*q_vec+sigma_vec[2]*S1*h_2*q_vec+
               sigma_vec[3]*S2*h_3+sigma_vec[4]*S3*h_3) 
  }
  
  # return hospitalizations for each time step.
  return(Hospital)
}