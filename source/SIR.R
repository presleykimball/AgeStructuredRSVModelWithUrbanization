# this function produces the set of derivatives at a given time point (t) and 
# previous step (States_vec)
SIR <- function(t,States_vec,parms){
  # convert States_vec back into matrix
  row_names <- sub(" .*", "", names(States_vec))  
  col_names <- sub(".* ", "", names(States_vec)) 
  States <- matrix(States_vec, nrow = length(unique(row_names)), 
                   ncol = length(unique(col_names)), byrow = FALSE,
                   dimnames = list(unique(row_names), unique(col_names)))
  
  # Call in fitting parameters
  alpha_1 <- parms[["alpha_1"]]
  alpha_2 <- parms[["alpha_2"]]
  alpha_3 <- parms[["alpha_3"]]
  
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

  if (is.null(parms$step_function) == FALSE){
    if (parms$step_function == TRUE){
    if (((t+alpha_3) %% (365.25/7))>(141.25/7) &((t+alpha_3) %% (365.25/7))<(171.25/7)){ # half transmission between [21,25) weeks
      alpha_1 <-alpha_1*0.5
  }}}
  # turn vaccination on or off under the following conditions:
  if(is.null(parms[["time_of_vacc"]])){ # if no time_of_vacc is given, set vaccination to 0
    chi_m <- 0
    chi_b <- 0
  } else{
    if (year_round_vacc == FALSE){
      time_of_vacc <- parms[["time_of_vacc"]]+90.25/7+3 # adding this so that vaccines are only implemented starting in the later half of the year i.e. at the beginning of an RSV season.
      # switch to turn vaccination on/off: If we are past the time_of_vacc, turn vaccination on; else, leave it off.
      if(t > time_of_vacc){
        if (is.null(parms[['vacc_shift']])){
          shift <- 0
        } else{
          shift <- parms[['vacc_shift']]
        }
        # call in vaccination coverage parameters. Set vaccination only to be on during RSV season.
        # From CDC, maternal vaccination is only given September through January;
        # baby vaccination is given October through March.
        if(t %% (365.25/7) <= 31/7+shift | t %% (365.25/7) >=242.25/7-shift){
          chi_m <- parms[["chi_m"]]
        } else{
          chi_m <- 0
        }
        if(t %% (365.25/7) <= 90.25/7+shift| t %% (365.25/7) >=272.25/7-shift){
          chi_b <- parms[["chi_b"]]
        } else{
          chi_b <- 0
        }
      } else {
        chi_m <- 0
        chi_b <- 0
      }
    } else{# vaccination on year round
      if(is.null(parms[["time_of_vacc"]])){ # if no time_of_vacc is given, set vaccination to 0
        chi_m <- 0
        chi_b <- 0
      }else{ # if time_of_vacc is given, then do the following
        time_of_vacc <- parms[["time_of_vacc"]]+272.25/7
        # switch to turn vaccination on/off: If we are past the time_of_vacc, turn vaccination on; else, leave it off.
        if(t > time_of_vacc){
          # turn vaccination on all the time
          chi_m <- parms[["chi_m"]]
          chi_b <- parms[["chi_b"]]
        }else{
          chi_m <- 0
          chi_b <- 0
        }
      }
    }
  }
  
  
  
  #Pull out the states  for the model as vectors
  S0 <-  States[,'S0']
  S1 <-  States[,'S1']
  S2 <-  States[,'S2']
  S3 <-  States[,'S3']
  I1 <-  States[,'I1']
  I2 <-  States[,'I2']
  I3 <-  States[,'I3']
  I4 <-  States[,'I4']
  R <-  States[,'R']
  Vm <-  States[,'Vm']
  Vb <-  States[,'Vb']
  
  # calculate infection rate (each entry of the vector is the infection rate for that age group)
  N_i <- rowSums(States) #number in each age category
  K_mod <- sweep(K,1,N_i,"/") # divide each row of contact matrix by number in that age_cat
  K_mod <- sweep(K_mod,2,tau_vec_comp,"*") # multiply each column by optional amplification tau
  K_mod <- sweep(K_mod,1,tau_vec_comp,"*") # multiply each row by optional amplification tau
  sums <- rho_vec[1]*(K_mod %*% I1)+rho_vec[2]*(K_mod %*% I2)+rho_vec[3]*(K_mod %*% I3)+rho_vec[4]*(K_mod %*% I4)
  b_vec <- alpha_1*(alpha_2*cos((2*pi*t-alpha_3)/365.25*7)+1)*sums
  
  ### Calculate derivatives
  
  # Create empty matrix for storage
  f <- States
  f[,]<- 0
  
  f[,'S0'] <- (1-chi_m)*Lambda+c(0,r_vec[1:(length(S0)-1)])*c(0,S0[1:(length(S0)-1)])+nu_b*Vb+nu_m*Vm-
    (sigma_vec[1]*b_vec+0.75*c(chi_b, rep(0,length(S0)-1))+r_vec+mu)*S0
  f[,'S1'] <- c(0,r_vec[1:(length(S1)-1)])*c(0,S1[1:(length(S1)-1)])+gamma_vec[1]*I1 -
    (sigma_vec[2]*b_vec+r_vec+mu)*S1
  f[,'S2'] <- c(0,r_vec[1:(length(S2)-1)])*c(0,S2[1:(length(S2)-1)])+gamma_vec[2]*I2 -
    (sigma_vec[3]*b_vec+r_vec+mu)*S2
  f[,'S3'] <- c(0,r_vec[1:(length(S3)-1)])*c(0,S3[1:(length(S3)-1)])+
    gamma_vec[3]*I3 +(1-p_R)*gamma_vec[4]*I4-
    (sigma_vec[4]*b_vec+r_vec+mu)*S3
  f[,'I1'] <- c(0,r_vec[1:(length(I1)-1)])*c(0,I1[1:(length(I1)-1)])+
    sigma_vec[1]*b_vec*S0+iota_b*b_vec*Vb+iota_m*b_vec*Vm -(gamma_vec[1]+r_vec+mu)*I1
  f[,'I2'] <- c(0,r_vec[1:(length(I2)-1)])*c(0,I2[1:(length(I2)-1)])+sigma_vec[2]*b_vec*S1 -
    (gamma_vec[2]+r_vec+mu)*I2
  f[,'I3'] <- c(0,r_vec[1:(length(I3)-1)])*c(0,I3[1:(length(I3)-1)])+sigma_vec[3]*b_vec*S2 - 
    (gamma_vec[3]+r_vec+mu)*I3
  f[,'I4'] <- c(0,r_vec[1:(length(I4)-1)])*c(0,I4[1:(length(I4)-1)])+sigma_vec[4]*b_vec*S3 - 
    (gamma_vec[4]+r_vec+mu)*I4
  f[,'R'] <- c(0,r_vec[1:(length(R)-1)])*c(0,R[1:(length(R)-1)])+p_R*gamma_vec[4]*I4-
    (r_vec+mu)*R
  f[,'Vm'] <- chi_m* Lambda + c(0,r_vec[1:(length(Vm)-1)])*c(0,Vm[1:(length(Vm)-1)]) - 
    (nu_m + iota_m*b_vec+r_vec+mu)*Vm
  f[,'Vb'] <- c(0,r_vec[1:(length(Vb)-1)])*c(0,Vb[1:(length(Vb)-1)]) + 
    0.75*c(chi_b, rep(0,length(S0)-1))*S0 - (nu_b +iota_b*b_vec +r_vec+mu)*Vb
  
  # return matrix to vector form so that ODE45 can read
  dStates <- as.vector(f)
  return(list(dStates))
}