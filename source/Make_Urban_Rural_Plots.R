Make_Urban_Rural_Plots <- function(urban_tau,rural_tau,param_set,urban_pop_size,
                                   rural_pop_size,State_name,urban_age_cat,rural_age_cat){
  
  # Clean params data
  param_set <- param_set[order(param_set$Loss, decreasing = FALSE),]
  parms <- as.list(param_set[1,])
  print(parms)
  
  # Get solution and hospitalizations for urban
  parms$tau_1 <- 0
  parms$tau_2 <- urban_tau
  urban_solution <- ODE_solve_SIR(parms,States_IC_vec,20,522)
  urban_hospitalizations <- Track_Hospitalizations(solution = urban_solution,parms = parms)
  urban_hospitalizations[,2:ncol(urban_hospitalizations)] <-
    urban_hospitalizations[,2:ncol(urban_hospitalizations)]*urban_pop_size
  
  # Get solution and hospitalizations for urban
  parms$tau_2 <- rural_tau
  rural_solution <- ODE_solve_SIR(parms,States_IC_vec,20,522)
  rural_hospitalizations <- Track_Hospitalizations(solution = rural_solution,parms = parms)
  rural_hospitalizations[,2:ncol(rural_hospitalizations)] <-
    rural_hospitalizations[,2:ncol(rural_hospitalizations)]*rural_pop_size
  
  # Modify data for plotting
  urban_age_hosps <- colSums(urban_hospitalizations[,2:ncol(urban_hospitalizations)])
  rural_age_hosps <- colSums(rural_hospitalizations[,2:ncol(rural_hospitalizations)])
  urban_age_hosps <- urban_age_hosps/sum(urban_age_hosps)
  rural_age_hosps <- rural_age_hosps/sum(rural_age_hosps)
  urban_age_hosps <- as.data.frame(urban_age_hosps)
  rural_age_hosps <- as.data.frame(rural_age_hosps)
  colnames(urban_age_hosps) <- "urban"
  colnames(rural_age_hosps) <- "rural"
  urban_age_hosps$age_cat <- age_cats
  rural_age_hosps$age_cat <- age_cats
  
  urban_prop <- urban_age_cat %>%
    group_by(age_cat)%>%
    summarise(data = sum(cases))
  urban_prop$data <- urban_prop$data/sum(urban_prop$data)
  urban_prop <- as.data.frame(urban_prop)
  colnames(urban_prop) <- c("age_cat","urban_data")
  
  rural_prop <- rural_age_cat %>%
    group_by(age_cat)%>%
    summarise(data = sum(cases))
  rural_prop$data <- rural_prop$data/sum(rural_prop$data)
  rural_prop <- as.data.frame(rural_prop)
  colnames(rural_prop) <- c("age_cat","rural_data")
  
  comparison <- merge(urban_age_hosps,rural_age_hosps,by= "age_cat")
  data_comparison <- merge(urban_prop,rural_prop, by= "age_cat")
  data_comparison <- merge(comparison,data_comparison, by = "age_cat")
  comparison <- comparison[match(age_cats,comparison$age_cat),]
  data_comparison <- data_comparison[match(age_cats,data_comparison$age_cat),]
  row.names(comparison) <- comparison[,1]
  row.names(data_comparison) <- data_comparison[,1]
  comparison <- comparison[,-1]
  data_comparison <- data_comparison[,-1]
  comparison <- as.matrix(comparison)
  data_comparison <- as.matrix(data_comparison)
  comparison <- apply(t(comparison),2,rev)
  data_comparison<- apply(t(data_comparison),2,rev)
  comparison <- comparison[c("urban","rural"),]
  data_comparison <- data_comparison[c("urban","urban_data","rural","rural_data"),]
  row.names(comparison) <- c("Predicted Urban","Predicted Rural")
  row.names(data_comparison) <- c("Predicted Urban","Data Urban","Predicted Rural", "Data Rural")
  
  barplot(data_comparison,beside = TRUE,legend.text = TRUE,
          ylim = c(0,0.8),
          main = paste("Proportion of Reported RSV by Age Group in", State_name,"Best Case"),
          col = c(palette_4[1],palette_4[2],palette_4[3],palette_4[4]),
          xlab = "Age Category",
          ylab = "Proportion of Cases",
          las =2)
  barplot(comparison,beside = TRUE,legend.text = TRUE,
          ylim = c(0,0.8),
          main = paste("Proportion of Reported RSV by Age Group in", State_name,
                       "Best Case"),
          col = c(palette_4[1],palette_4[2]),
          xlab = "Age Category",
          ylab = "Proportion of Cases",
          las =2)
  
  data_comparison_young <- data_comparison[,colnames(data_comparison) %in% c("0-1","1-2","2-3","3-4",
                                                                             "4-5","5-10")]
  comparison_young <- comparison[,colnames(comparison) %in% c("0-1","1-2","2-3","3-4",
                                                              "4-5","5-10")]
  barplot(data_comparison_young,beside = TRUE,legend.text = TRUE,
          ylim = c(0,0.8),
          main = paste("Proportion of Reported RSV by Age Group in", State_name,"Best Case"),
          col = c(palette_4[1],palette_4[2],palette_4[3],palette_4[4]),
          xlab = "Age Category",
          ylab = "Proportion of Cases")
  barplot(comparison_young,beside = TRUE,legend.text = TRUE,
          ylim = c(0,0.8),
          main = paste("Proportion of Reported RSV by Age Group in", State_name,"Best Case"),
          col = c(palette_4[1],palette_4[2]),
          xlab = "Age Category",
          ylab = "Proportion of Cases")
}