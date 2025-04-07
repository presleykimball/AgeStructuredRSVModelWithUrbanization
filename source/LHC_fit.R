# This applies the LHC fitting algorithm as described in the main manuscript
LHC_fit <- function(data_age_cat,pop_size){
  ####### FIRST FITTING ALPHA_{1,2,3} ###########
  # Set seed and sample size for first LHC fitting
  sample_size <- 1000
  
  # Draw 3 RS on unif(0,1)
  parm_set_1 <- randomLHS(sample_size,3)
  parm_set_1 <- as.data.frame(parm_set_1)
  colnames(parm_set_1) <- c("alpha_1","alpha_2","alpha_3")
  
  # rescale parameter values to be within an expected range.
  parm_set_1$alpha_1<-1.5*parm_set_1$alpha_1
  parm_set_1$alpha_2<-0.5*parm_set_1$alpha_2
  parm_set_1$alpha_3<-pi/3*(365.25/7)*parm_set_1$alpha_3+pi/3*(365.25/7)
  
  # fix the rescaling parameters
  parm_set_1$q_1 <- rep(1,nrow(parm_set_1))
  parm_set_1$q_2 <- rep(1,nrow(parm_set_1))
  parm_set_1$q_3 <- rep(1,nrow(parm_set_1))
  parm_set_1$q_4 <- rep(1,nrow(parm_set_1))
  parm_set_1$q_5 <- rep(1,nrow(parm_set_1))
  parm_set_1$q_6 <- rep(1,nrow(parm_set_1))
  
  # fix daycare amplification parm
  parm_set_1$tau_1 <- 0
  parm_set_1$tau_2<- 1
  
  # set vaccination to 0
  parm_set_1$chi_b <- 0
  parm_set_1$chi_m <- 0
  # create place for storage
  parm_set_1$Loss <- numeric(nrow(parm_set_1))
  
  # find the loss of each parameter combination
  for (i in 1:nrow(parm_set_1)){
    trial_parms <- as.list(parm_set_1[i,])
    parm_set_1$Loss[i] <- fitmodel(trial_parms,data_age_cat,pop_size)
  }
  
  parm_set_1$Loss <- as.numeric(parm_set_1$Loss)
  parm_set_1 <- parm_set_1[order(parm_set_1$Loss,decreasing = FALSE),]
  
  ####### SECOND REFINED FITTING ALPHA_{1,2,3} ###########
  top_parms <- parm_set_1[1:(0.05*sample_size),]
  
  # Draw 3 RS on unif(0,1) 1000 times
  parm_set_2 <- randomLHS(sample_size,3)
  parm_set_2 <- as.data.frame(parm_set_2)
  colnames(parm_set_2) <- c("alpha_1","alpha_2","alpha_3")
  
  # rescale parameter values to be within an expected range.
  parm_set_2$alpha_1<-(max(top_parms$alpha_1)-min(top_parms$alpha_1))*parm_set_2$alpha_1+
    min(top_parms$alpha_1)
  parm_set_2$alpha_2<-(max(top_parms$alpha_2)-min(top_parms$alpha_2))*parm_set_2$alpha_2+
    min(top_parms$alpha_2)
  parm_set_2$alpha_3<-(max(top_parms$alpha_3)-min(top_parms$alpha_3))*parm_set_2$alpha_3+
    min(top_parms$alpha_3)
  
  # fix the rescaling parameters
  parm_set_2$q_1 <- rep(1,nrow(parm_set_2))
  parm_set_2$q_2 <- rep(1,nrow(parm_set_2))
  parm_set_2$q_3 <- rep(1,nrow(parm_set_2))
  parm_set_2$q_4 <- rep(1,nrow(parm_set_2))
  parm_set_2$q_5 <- rep(1,nrow(parm_set_2))
  parm_set_2$q_6 <- rep(1,nrow(parm_set_2))
  
  # fix daycare amplification parm
  parm_set_2$tau_1 <- 0
  parm_set_2$tau_2 <- 1
  
  # set vaccination to 0
  parm_set_2$chi_b <- 0
  parm_set_2$chi_m <- 0
  
  # create place for storage
  parm_set_2$Loss <- numeric(nrow(parm_set_2))
  
  # find the loss of each parameter combination
  for (i in 1:nrow(parm_set_2)){
    trial_parms <- as.list(parm_set_2[i,])
    parm_set_2$Loss[i] <- fitmodel(trial_parms,data_age_cat,pop_size)
  }
  
  parm_set_2$Loss <- as.numeric(parm_set_2$Loss)
  parm_set_2 <- parm_set_2[order(parm_set_2$Loss,decreasing = FALSE),]
  ####### FITTING OF q VECTOR ###########
  # Draw 3 RS on unif(0,1)
  parm_set_3 <- randomLHS(sample_size*2,6)
  parm_set_3 <- as.data.frame(parm_set_3)
  colnames(parm_set_3) <- c("q_1","q_2","q_3","q_4","q_5","q_6")
  
  # rescale parameter values to be within an expected range.
  parm_set_3$q_1 <- parm_set_3$q_1+0.5
  parm_set_3$q_2 <- parm_set_3$q_2*2+1 # search [1,3]
  parm_set_3$q_3 <- parm_set_3$q_3*2.5+1.5 # search [1.5,4]
  parm_set_3$q_4 <- parm_set_3$q_4*2.5+1.5 # search [1.5,4]
  parm_set_3$q_5 <- parm_set_3$q_5*4+3 # search [3,7]
  parm_set_3$q_6 <- parm_set_3$q_6*5+10 # search[10,15]
  
  # fix the the infectious parameters
  parm_set_3$alpha_1<-rep(parm_set_2$alpha_1[1],nrow(parm_set_3))
  parm_set_3$alpha_2<-rep(parm_set_2$alpha_2[1],nrow(parm_set_3))
  parm_set_3$alpha_3<-rep(parm_set_2$alpha_3[1],nrow(parm_set_3))
  
  # fix daycare amplification parm
  parm_set_3$tau_1 <- 0
  parm_set_3$tau_2 <- 1
  
  # set vaccination to 0
  parm_set_3$chi_b <- 0
  parm_set_3$chi_m <- 0
  
  # create place for storage
  parm_set_3$Loss <- numeric(nrow(parm_set_3))
  
  # find the loss of each parameter combination
  for (i in 1:nrow(parm_set_3)){
    trial_parms <- as.list(parm_set_3[i,])
    parm_set_3$Loss[i] <- fitmodel(trial_parms,data_age_cat,pop_size)
  }
  
  parm_set_3$Loss <- as.numeric(parm_set_3$Loss)
  parm_set_3 <- parm_set_3[order(parm_set_3$Loss,decreasing = FALSE),]
  
  ####### LAST FITTING ALPHA_{1,2,3} ###########
  top_parms <- parm_set_2[1:(0.025*sample_size),]
  top_scaling <- parm_set_3[1:(0.05*sample_size),]
  # Draw 3 RS on unif(0,1) 1000 times
  parm_set_4 <- randomLHS(sample_size*10,9)
  parm_set_4 <- as.data.frame(parm_set_4)
  colnames(parm_set_4) <- c("alpha_1","alpha_2","alpha_3","q_1","q_2","q_3","q_4","q_5","q_6")
  
  # rescale parameter values to be within a range
  parm_set_4$alpha_1<-(max(top_parms$alpha_1)-min(top_parms$alpha_1))*parm_set_4$alpha_1+
    min(top_parms$alpha_1)
  parm_set_4$alpha_2<-(max(top_parms$alpha_2)-min(top_parms$alpha_2))*parm_set_4$alpha_2+
    min(top_parms$alpha_2)
  parm_set_4$alpha_3<-(max(top_parms$alpha_3)-min(top_parms$alpha_3))*parm_set_4$alpha_3+
    min(top_parms$alpha_3)
  parm_set_4$q_1 <- (max(top_scaling$q_1)-min(top_scaling$q_1))*parm_set_4$q_1+
    min(top_scaling$q_1)
  parm_set_4$q_2 <- (max(top_scaling$q_2)-min(top_scaling$q_2))*parm_set_4$q_2+
    min(top_scaling$q_2)
  parm_set_4$q_3 <- (max(top_scaling$q_3)-min(top_scaling$q_3))*parm_set_4$q_3+
    min(top_scaling$q_3)
  parm_set_4$q_4 <- (max(top_scaling$q_4)-min(top_scaling$q_4))*parm_set_4$q_4+
    min(top_scaling$q_4)
  parm_set_4$q_5 <- (max(top_scaling$q_5)-min(top_scaling$q_5))*parm_set_4$q_5+
    min(top_scaling$q_5)
  parm_set_4$q_6 <- (max(top_scaling$q_6)-min(top_scaling$q_6))*parm_set_4$q_6+
    min(top_scaling$q_6)
  
  # fix daycare amplification parm
  parm_set_4$tau_1 <- 0
  parm_set_4$tau_2 <- 1
  
  # set vaccination to 0
  parm_set_4$chi_b <- 0
  parm_set_4$chi_m <- 0
  
  # create place for storage
  parm_set_4$Loss <- numeric(nrow(parm_set_4))
  
  # find the loss of each parameter combination
  for (i in 1:nrow(parm_set_4)){
    trial_parms <- as.list(parm_set_4[i,])
    parm_set_4$Loss[i] <- fitmodel(trial_parms,data_age_cat,pop_size)
  }
  
  parm_set_4$Loss <- as.numeric(parm_set_4$Loss)
  parm_set_4 <- parm_set_4[order(parm_set_4$Loss,decreasing = FALSE),]
  
  #Return fitted param sets
  return(list(parm_set_1,parm_set_2,parm_set_3,parm_set_4))
}