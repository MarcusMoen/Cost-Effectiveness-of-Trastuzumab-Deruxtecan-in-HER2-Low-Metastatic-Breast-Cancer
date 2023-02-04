




#NOTE:
#The transitions from rate to probability is found in https://journals.sagepub.com/doi/pdf/10.1177/0272989X05282637 Figure 3


#Function for calculating the standard error of a survival curve given 
#number at risk, n_risk, number of deaths at each timestep, n_deaths, and the survival probabilities, surv_prob


library(expm)

#TODO: Function currently not in use. Need to find n_deaths first.
calc_serror <- function(n_risk, n_deaths, surv_prob){
  serror <- c()
  sum_temp <- 0
  for(i in 1:length(n_risk)){
    sum_temp <- sum_temp + n_deaths[i]/(n_risk[i]*(n_risk[i]-n_deaths[i]))
    serror <- append(serror, surv_prob*sqrt(sum_temp))
  }
  return(serror)
}


km_pf_chemo <- c(1, 0.983, 0.756, 0.62, 0.603, 0.501, 0.439, 0.398, 0.34, 0.277, 0.268, 0.226, 0.217, 0.185) 
n_risk_pf_chemo <- c(184, 166, 119, 93, 90, 73, 60, 51, 45, 34, 32, 29, 26, 22)
n_deaths_pf_chemo <- c()#?????










#Function for finding the transition matrices such that they can be imported to the CEA file
#Takes as its input the optimal values for the eval_variables function and outputs a list of transition matrices
transition_matrices <- function(input_var, calibrate=T){
  #print("RUNNING")
  
  r_prog <- input_var[1] #Rate 1
  r_prog_death <- input_var[2] #Rate 2
  hr_pf2prog_chemo2tdxd <- input_var[3]
  
  gamma_chemo <- input_var[4]
  alpha_chemo <- input_var[5]
  
  gamma_tdxd <- input_var[6]
  alpha_tdxd <- input_var[7]
  
  gamma_ILD <- input_var[8]
  alpha_ILD <- input_var[9]
  
  #Rate of dying from other causes (background mortality)
  r_OC_death <- rep(c(0.000455117,0.000667712,0.000939403,0.001487732,0.002481213, 0.004368355, 0.008008381,0.014505686, 0.024617632)
                    ,each=(12*5))

  if(calibrate){
    n <- 40
  }else{
    n <- length(r_OC_death)
  }
  #The rate of getting an AE and being discontinued. (Needs updating)
  Ft_chemo <- c()
  Ft_tdxd <- c()
  Ft_tdxd_ILD <- c()
  
  Ht_chemo <- c()
  Ht_tdxd <- c()
  Ht_tdxd_ILD <- c()
  
  r_AE <- c()
  r_AE_tdxd <- c()
  r_ILD <- c()
  r_ILD_tdxd <- c()
  for(t in 1:(n+1)){
    Ft_chemo <- append(Ft_chemo, alpha_chemo*(1-exp(-gamma_chemo*t)))
    Ft_tdxd <- append(Ft_tdxd, alpha_tdxd*(1-exp(-gamma_tdxd*t)))
    Ft_tdxd_ILD <- append(Ft_tdxd_ILD, alpha_ILD*(1-exp(-gamma_ILD*t)))
    
    Ht_chemo <- append(Ht_chemo, -log(1-Ft_chemo[t]))
    Ht_tdxd <- append(Ht_tdxd, -log(1-Ft_tdxd[t]))
    Ht_tdxd_ILD <- append(Ht_tdxd_ILD, -log(1-Ft_tdxd_ILD[t]))
    
    if(t > 1){
      r_AE <- append(r_AE, Ht_chemo[t]-Ht_chemo[t-1])
      r_AE_tdxd <- append(r_AE_tdxd, Ht_tdxd[t]-Ht_tdxd[t-1])
      r_ILD <- append(r_ILD, 0)
      r_ILD_tdxd <- append(r_ILD_tdxd, Ht_tdxd_ILD[t]-Ht_tdxd_ILD[t-1])
    }
  }
  
  
  #Defining some hazard ratios (Dummy variables for now, maybe changed later)
  hr_progAE <- 1
  hr_progAE_death <- 1
  hr_pfae2prog_chemo2tdxd <- 1
  hr_prog2death_chemo2tdxd <- 1
  hr_progAE2death_chemo2tdxd <- 1
  hr_pf2progILD_chemo2tdxd <- 1
  hr_progILD2death_chemo2tdxd <- 1
  
  
  #Finding the relationship between the rates
  r_progAE <- r_prog*hr_progAE
  r_progILD <- r_prog*0
  r_progILD_death <- r_prog_death*0
  r_progAE_death <- r_prog_death*hr_progAE_death
  
  
  r_prog_tdxd <- r_prog*hr_pf2prog_chemo2tdxd
  r_progAE_tdxd <- r_prog*hr_pfae2prog_chemo2tdxd
  r_prog_death_tdxd <- r_prog_death*hr_prog2death_chemo2tdxd
  r_progAE_death_tdxd <- r_prog_death*hr_progAE2death_chemo2tdxd
  r_progILD_tdxd <- r_prog*hr_pf2progILD_chemo2tdxd
  r_progILD_death_tdxd <- r_prog_death*hr_progILD2death_chemo2tdxd
  
  P_chemo <- list()
  P_tdxd <- list()
  
  for(t in 1:n){
    #Creating the transition matrices
    G_chemo <- matrix(c(-(r_prog + r_OC_death[t] + r_AE[t] + r_ILD[t]), r_AE[t], r_ILD[t], r_prog, 0, 0, r_OC_death[t],
                        0, -(r_progAE + r_OC_death[t]), 0, 0, r_progAE, 0, r_OC_death[t],
                        0, 0, -(r_progILD + r_OC_death[t]), 0, 0, r_progILD, r_OC_death[t],
                        0, 0, 0, -(r_OC_death[t] + r_prog_death), 0, 0, (r_OC_death[t] + r_prog_death),
                        0, 0, 0, 0, -(r_OC_death[t] + r_progAE_death), 0, (r_OC_death[t] + r_progAE_death),
                        0, 0, 0, 0, 0, -(r_OC_death[t] + r_progILD_death), (r_OC_death[t] + r_progILD_death),
                        0, 0, 0, 0, 0, 0, 0), 
                      ncol = 7, nrow = 7,
                      dimnames = list(c("PF", "PF_AE", "PF_ILD", "Prog", "Prog_AE", "Prog_ILD", "Death"),
                                      c("PF", "PF_AE", "PF_ILD", "Prog", "Prog_AE", "Prog_ILD", "Death")),
                      byrow = T
    )
    #print("This is G_chemo:")
    #print(G_chemo)
    
    G_tdxd <- matrix(c(-(r_prog_tdxd + r_OC_death[t] + r_AE_tdxd[t] + r_ILD_tdxd[t]), r_AE_tdxd[t], r_ILD_tdxd[t], r_prog_tdxd, 0, 0, r_OC_death[t],
                        0, -(r_progAE_tdxd + r_OC_death[t]), 0, 0, r_progAE_tdxd, 0, r_OC_death[t],
                        0, 0, -(r_progILD_tdxd + r_OC_death[t]), 0, 0, r_progILD_tdxd, r_OC_death[t],
                        0, 0, 0, -(r_OC_death[t] + r_prog_death_tdxd), 0, 0, (r_OC_death[t] + r_prog_death_tdxd),
                        0, 0, 0, 0, -(r_OC_death[t] + r_progAE_death_tdxd), 0, (r_OC_death[t] + r_progAE_death_tdxd),
                        0, 0, 0, 0, 0, -(r_OC_death[t] + r_progILD_death_tdxd), (r_OC_death[t] + r_progILD_death_tdxd),
                        0, 0, 0, 0, 0, 0, 0), 
                      ncol = 7, nrow = 7,
                      dimnames = list(c("PF", "PF_AE", "PF_ILD", "Prog", "Prog_AE", "Prog_ILD", "Death"),
                                      c("PF", "PF_AE", "PF_ILD", "Prog", "Prog_AE", "Prog_ILD", "Death")),
                      byrow = T
    )
    
    A_chemo <- expm(G_chemo)
    #print("This is A_chemo:")
    #print(A_chemo)
    A_tdxd <- expm(G_tdxd)
    
    P_chemo[[t]] <- A_chemo
    P_tdxd[[t]] <- A_tdxd
    
  }
  
  matrices <- list(P_chemo, P_tdxd)
  return(matrices)
}



#Input:
#input_var = a list of three values, the rate from PF to P, the rate from P to D, and the HR.
#Output:
#The square error when using the given inputs.
eval_variables <- function(input_var){
  
  
  r_prog <- input_var[1] #Rate 1
  r_prog_death <- input_var[2] #Rate 2
  hr_pf2prog_chemo2tdxd <- input_var[3]
  
  gamma_chemo <- input_var[4]
  alpha_chemo <- input_var[5]
  
  gamma_tdxd <- input_var[6]
  alpha_tdxd <- input_var[7]
  
  gamma_ILD <- input_var[8]
  alpha_ILD <- input_var[9]
  
  
  matrix_list <- transition_matrices(input_var)
  A_chemo_list <- matrix_list[[1]]
  A_tdxd_list <- matrix_list[[2]]
  
  
  
  #The extracted Kaplan-Meier values
  km_os_chemo <- c(1, 0.986, 0.981, 0.957, 0.94, 0.927, 0.884, 0.841, 0.793, 0.738, 0.712, 0.683, 0.669, 0.634)#, 0.610, 0.566, 0.520, 0.484, 0.459, 0.460, 0.457, 0.419, 0.375)
  km_pf_chemo <- c(1, 0.983, 0.756, 0.62, 0.603, 0.501, 0.439, 0.398, 0.34, 0.277, 0.268, 0.226, 0.217, 0.185) 
  km_os_tdxd <- c(1, 0.992, 0.985, 0.967, 0.958, 0.939, 0.928, 0.899, 0.87, 0.857, 0.829, 0.815, 0.793, 0.774)#, 0.742, 0.731, 0.688, 0.661, 0.623, 0.586, 0.566, 0.539, 0.507, 0.500, 0.484)
  km_pf_tdxd <- c(1, 0.988, 0.896, 0.82, 0.809, 0.765, 0.681, 0.633, 0.592, 0.554, 0.491, 0.457, 0.422, 0.386)#, 0.371, 0.363, 0.339, 0.31, 0.295, 0.274, 0.259) 
  
  cum_AE_chemo <- c(0, 0.0145, 0.0266, 0.0358, 0.044, 0.05, 0.056, 0.0603, 0.064, 0.067, 0.069, 0.071, 0.073, 0.074, 0.075, 0.076, 0.077, 0.0774, 0.0778)
  cum_AE_tdxd <- c(0, 0.0145, 0.0266, 0.0358, 0.044, 0.05, 0.056, 0.0603, 0.064, 0.067, 0.069, 0.071, 0.073, 0.074, 0.075, 0.076, 0.077, 0.0774, 0.0778)
  cum_ILD <- c(0, 0.0145, 0.0266, 0.0358, 0.044, 0.05, 0.056, 0.0603, 0.064, 0.067, 0.069, 0.071, 0.073, 0.074, 0.075, 0.076, 0.077, 0.0774, 0.0778)
  
  #Vectors for storing the estimates of the KM curves
  km_os_chemo_model <- c()
  km_pf_chemo_model <- c()
  km_os_tdxd_model <- c()
  km_pf_tdxd_model <- c() 
  
  cum_AE_chemo_model <- c()
  cum_AE_tdxd_model <- c()
  cum_ILD_model <- c()
  
  
  n = max(length(km_os_chemo), length(km_pf_chemo), length(km_os_tdxd), length(km_pf_tdxd), length(cum_ILD)) #Find the maximum length
  
  #Start states
  s0_chemo <- c(1,0,0,0,0,0,0)
  s0_tdxd <- c(1,0,0,0,0,0,0)
  
  survival_chemo <- s0_chemo[1]+s0_chemo[4]
  km_os_chemo_model <- append(km_os_chemo_model, survival_chemo)
  km_pf_chemo_model <- append(km_pf_chemo_model, s0_chemo[1])
  
  survival_tdxd <- s0_tdxd[1]+s0_tdxd[4]
  km_os_tdxd_model <- append(km_os_tdxd_model, survival_tdxd)
  km_pf_tdxd_model <- append(km_pf_tdxd_model, s0_tdxd[1])
  
  cum_AE_chemo_model <- append(cum_AE_chemo_model, 0)
  cum_AE_tdxd_model <- append(cum_AE_tdxd_model, 0)
  cum_ILD_model <- append(cum_ILD_model, 0)
  
  
  #Do the simulations
  for(t in 1:n){
    
    A_chemo <- A_chemo_list[[t]]
    A_tdxd <- A_tdxd_list[[t]]
    
    
    
    s1_chemo <- s0_chemo %*% A_chemo
    s1_tdxd <- s0_tdxd %*% A_tdxd
    
    #survival_chemo <- s1_chemo[1]+s1_chemo[4]
    #km_os_chemo_model <- append(km_os_chemo_model, survival_chemo)
    #km_pf_chemo_model <- append(km_pf_chemo_model, s1_chemo[1])
    
    #survival_tdxd <- s1_tdxd[1]+s1_tdxd[4]
    #km_os_tdxd_model <- append(km_os_tdxd_model, survival_tdxd)
    #km_pf_tdxd_model <- append(km_pf_tdxd_model, s1_tdxd[1])
    
    survival_chemo <- km_os_chemo_model[t]*(s1_chemo[1]+s1_chemo[4])/(s1_chemo[1]+s1_chemo[4]+s0_chemo[1]*A_chemo[1,7]+s0_chemo[4]*A_chemo[4,7])
    km_os_chemo_model <- append(km_os_chemo_model, survival_chemo)
    km_pf_chemo_model <- append(km_pf_chemo_model, km_pf_chemo_model[t]*s1_chemo[1]/(s1_chemo[1]+s0_chemo[1]*A_chemo[1,4]+s0_chemo[1]*A_chemo[1,7]))
    
    survival_tdxd <- km_os_tdxd_model[t]*(s1_tdxd[1]+s1_tdxd[4])/(s1_tdxd[1]+s1_tdxd[4]+s0_tdxd[1]*A_tdxd[1,7]+s0_tdxd[4]*A_tdxd[4,7])
    km_os_tdxd_model <- append(km_os_tdxd_model, survival_tdxd)
    km_pf_tdxd_model <- append(km_pf_tdxd_model, km_pf_tdxd_model[t]*(s1_tdxd[1])/(s1_tdxd[1]+s0_tdxd[1]*A_tdxd[1,4]+s0_tdxd[1]*A_tdxd[1,7]))
    
    cum_AE_chemo_model <- append(cum_AE_chemo_model, cum_AE_chemo_model[t]+s0_chemo[1]*A_chemo[1,2])
    cum_AE_tdxd_model <- append(cum_AE_tdxd_model, cum_AE_tdxd_model[t]+s0_tdxd[1]*A_tdxd[1,2])
    cum_ILD_model <- append(cum_ILD_model, cum_ILD_model[t]+s0_tdxd[1]*A_tdxd[1,3])
    
    
    s0_chemo <- s1_chemo
    s0_tdxd <- s1_tdxd
    
  }
  
  
  
  #Calculate the error
  #serror <- sum((km_os_chemo[2:(length(km_os_chemo))] - km_os_chemo_model[2:(length(km_os_chemo))])^2) + 
   # sum((km_pf_chemo[2:(length(km_pf_chemo))] - km_pf_chemo_model[2:(length(km_pf_chemo))])^2) + 
    #sum((km_os_tdxd[2:(length(km_os_tdxd))] - km_os_tdxd_model[2:(length(km_os_tdxd))])^2) + 
    #sum((km_pf_tdxd[2:(length(km_pf_tdxd))] - km_pf_tdxd_model[2:(length(km_pf_tdxd))])^2) +
    #sum((cum_AE_chemo[2:(length(cum_AE_chemo_model))] - cum_AE_chemo_model[2:(length(cum_AE_chemo_model))])^2) +
    #sum((cum_AE_tdxd[2:(length(cum_AE_tdxd_model))] - cum_AE_tdxd_model[2:(length(cum_AE_tdxd_model))])^2) +
    #sum((cum_ILD[2:(length(cum_ILD_model))] - cum_ILD_model[2:(length(cum_ILD_model))])^2)
  
  
  
  serror <- sum(abs(km_os_chemo[2:(length(km_os_chemo))] - km_os_chemo_model[2:(length(km_os_chemo))])) + 
    sum(abs(km_pf_chemo[2:(length(km_pf_chemo))] - km_pf_chemo_model[2:(length(km_pf_chemo))])) + 
    sum(abs(km_os_tdxd[2:(length(km_os_tdxd))] - km_os_tdxd_model[2:(length(km_os_tdxd))])) + 
    sum(abs(km_pf_tdxd[2:(length(km_pf_tdxd))] - km_pf_tdxd_model[2:(length(km_pf_tdxd))])) +
    sum(abs(cum_AE_chemo[2:(length(cum_AE_chemo))] - cum_AE_chemo_model[2:(length(cum_AE_chemo))])) +
    sum(abs(cum_AE_tdxd[2:(length(cum_AE_tdxd))] - cum_AE_tdxd_model[2:(length(cum_AE_tdxd))])) +
    sum(abs(cum_ILD[2:(length(cum_ILD))] - cum_ILD_model[2:(length(cum_ILD))]))
  
  #print(sum(abs(km_os_chemo[2:(length(km_os_chemo))] - km_os_chemo_model[2:(length(km_os_chemo))])))
  #print(sum(abs(km_pf_chemo[2:(length(km_pf_chemo))] - km_pf_chemo_model[2:(length(km_pf_chemo))])))
  #print(sum(abs(km_os_tdxd[2:(length(km_os_tdxd))] - km_os_tdxd_model[2:(length(km_os_tdxd))])))
  #print(sum(abs(km_pf_tdxd[2:(length(km_pf_tdxd))] - km_pf_tdxd_model[2:(length(km_pf_tdxd))])))
  #print(sum(abs(cum_AE_chemo[2:(length(cum_AE_chemo))] - cum_AE_chemo_model[2:(length(cum_AE_chemo))])))
  #print(sum(abs(cum_AE_tdxd[2:(length(cum_AE_tdxd))] - cum_AE_tdxd_model[2:(length(cum_AE_tdxd))])))
  #print(sum(abs(cum_ILD[2:(length(cum_ILD))] - cum_ILD_model[2:(length(cum_ILD))])))
  
  return(serror)
  
  
  
}

#Find the optimal values (Here we could do a wider search for starting values since we might have missed the global optimum.)
fit_out <- optim(c(0.15, 0.08, 0.5, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2), 
                 eval_variables,
                 hessian = T)


#Extract the optimal values for our three variables
opt_var <- fit_out$par

opt_var



#Find the optimal transition matrices
tm <- transition_matrices(opt_var, calibrate=F)







#Plot the correct KM-curves vs the cureves from our model. This part of the file is used for validation
library(ggplot2)
library(ggpubr)



#The transitions from rate to probability is found in https://journals.sagepub.com/doi/pdf/10.1177/0272989X05282637 Figure 3
plot_function <- function(input_var){
  
  r_pf2p_chemo <- input_var[1] 
  r_p2d_chemo <- input_var[2] 
  hr <- input_var[3]
  
  gamma_chemo <- input_var[4]
  alpha_chemo <- input_var[5]
  
  gamma_tdxd <- input_var[6]
  alpha_tdxd <- input_var[7]
  
  gamma_ILD <- input_var[8]
  alpha_ILD <- input_var[9]
  
  matrix_list <- transition_matrices(input_var)
  A_chemo_list <- matrix_list[[1]]
  A_tdxd_list <- matrix_list[[2]]
  
  
  #The extracted Kaplan-Meier values
  km_os_chemo <- c(1, 0.986, 0.981, 0.957, 0.94, 0.927, 0.884, 0.841, 0.793, 0.738, 0.712, 0.683, 0.669, 0.634)#, 0.610, 0.566, 0.520, 0.484, 0.459, 0.460, 0.457, 0.419, 0.375)
  km_pf_chemo <- c(1, 0.983, 0.756, 0.62, 0.603, 0.501, 0.439, 0.398, 0.34, 0.277, 0.268, 0.226, 0.217, 0.185) 
  km_os_tdxd <- c(1, 0.992, 0.985, 0.967, 0.958, 0.939, 0.928, 0.899, 0.87, 0.857, 0.829, 0.815, 0.793, 0.774)#, 0.742, 0.731, 0.688, 0.661, 0.623, 0.586, 0.566, 0.539, 0.507, 0.500, 0.484)
  km_pf_tdxd <- c(1, 0.988, 0.896, 0.82, 0.809, 0.765, 0.681, 0.633, 0.592, 0.554, 0.491, 0.457, 0.422, 0.386)#, 0.371, 0.363, 0.339, 0.31, 0.295, 0.274, 0.259) 
  
  #Vectors for storing the estimates of the KM curves
  km_os_chemo_model <- c()
  km_pf_chemo_model <- c()
  km_os_tdxd_model <- c()
  km_pf_tdxd_model <- c() 
  
  
  n = max(length(km_os_chemo), length(km_pf_chemo), length(km_os_tdxd), length(km_pf_tdxd))
  
  #Start states
  s0_chemo <- c(1,0,0,0,0,0,0)
  s0_tdxd <- c(1,0,0,0,0,0,0)
  
  survival_chemo <- s0_chemo[1]+s0_chemo[4]
  km_os_chemo_model <- append(km_os_chemo_model, survival_chemo)
  km_pf_chemo_model <- append(km_pf_chemo_model, s0_chemo[1])
  
  survival_tdxd <- s0_tdxd[1]+s0_tdxd[4]
  km_os_tdxd_model <- append(km_os_tdxd_model, survival_tdxd)
  km_pf_tdxd_model <- append(km_pf_tdxd_model, s0_tdxd[1])
  
  ae_test_chemo <- c(0)
  ae_test_tdxd <- c(0)
  
  
  for(t in 1:n){
    A_chemo <- A_chemo_list[[t]]
    A_tdxd <- A_tdxd_list[[t]]
    
    s1_chemo <- s0_chemo %*% A_chemo
    s1_tdxd <- s0_tdxd %*% A_tdxd
    
    #survival_chemo <- s1_chemo[1]+s1_chemo[4]
    #km_os_chemo_model <- append(km_os_chemo_model, survival_chemo)
    #km_pf_chemo_model <- append(km_pf_chemo_model, s1_chemo[1])
    
    #survival_tdxd <- s1_tdxd[1]+s1_tdxd[4]
    #km_os_tdxd_model <- append(km_os_tdxd_model, survival_tdxd)
    #km_pf_tdxd_model <- append(km_pf_tdxd_model, s1_tdxd[1])
    
    survival_chemo <- km_os_chemo_model[t]*(s1_chemo[1]+s1_chemo[4])/(s1_chemo[1]+s1_chemo[4]+s0_chemo[1]*A_chemo[1,7]+s0_chemo[4]*A_chemo[4,7])
    km_os_chemo_model <- append(km_os_chemo_model, survival_chemo)
    km_pf_chemo_model <- append(km_pf_chemo_model, km_pf_chemo_model[t]*s1_chemo[1]/(s1_chemo[1]+s0_chemo[1]*A_chemo[1,4]+s0_chemo[1]*A_chemo[1,7]))
    
    survival_tdxd <- km_os_tdxd_model[t]*(s1_tdxd[1]+s1_tdxd[4])/(s1_tdxd[1]+s1_tdxd[4]+s0_tdxd[1]*A_tdxd[1,7]+s0_tdxd[4]*A_tdxd[4,7])
    km_os_tdxd_model <- append(km_os_tdxd_model, survival_tdxd)
    km_pf_tdxd_model <- append(km_pf_tdxd_model, km_pf_tdxd_model[t]*(s1_tdxd[1])/(s1_tdxd[1]+s0_tdxd[1]*A_tdxd[1,4]+s0_tdxd[1]*A_tdxd[1,7]))
    
    ae_test_chemo <- append(ae_test_chemo, ae_test_chemo[t]+s0_chemo[1]*A_chemo[1,2])
    ae_test_tdxd <- c(ae_test_tdxd, ae_test_tdxd[t]+s0_tdxd[1]*A_tdxd[1,2])
    
    s0_chemo <- s1_chemo
    s0_tdxd <- s1_tdxd
    
  }
  
  
  
  m = min(length(km_os_chemo), length(km_pf_chemo), length(km_os_tdxd), length(km_pf_tdxd))
  idx <- 1:m
  df <- data.frame(idx, km_os_chemo_model[1:m], km_os_tdxd_model[1:m], km_os_chemo[1:m])
  
  
  
  
  plot1 <- ggplot(df, aes(x=idx)) + 
    geom_line(aes(y = km_os_chemo_model[1:m]), color = "red") + 
    geom_line(aes(y = km_os_chemo[1:m]), color="blue", linetype="twodash") +
    ylim(0, 1) +
    ylab("Probability") +
    xlab("Month") +
    ggtitle("OS physician's choice")+
    scale_x_continuous(breaks = round(seq(0, 18, by = 3),1))+theme(axis.text=element_text(size=12),
                                                                   axis.title=element_text(size=14))
  
  plot2 <- ggplot(df, aes(x=idx)) + 
    geom_line(aes(y = km_os_tdxd_model[1:m]), color = "red") + 
    geom_line(aes(y = km_os_tdxd[1:m]), color="blue", linetype="twodash") +
    ylim(0, 1) +
    ylab("Probability")+
    xlab("Month")+
    ggtitle("OS T-DxD")+
    scale_x_continuous(breaks = round(seq(0, 18, by = 3),1))+theme(axis.text=element_text(size=12),
                                                                   axis.title=element_text(size=14))
  
  plot3 <- ggplot(df, aes(x=idx)) + 
    geom_line(aes(y = km_pf_chemo_model[1:m]), color = "red") + 
    geom_line(aes(y = km_pf_chemo[1:m]), color="blue", linetype="twodash") +
    ylim(0, 1) +
    ylab("Probability")+
    xlab("Month")+
    ggtitle("PFS physician's choice")+
    scale_x_continuous(breaks = round(seq(0, 18, by = 3),1))+theme(axis.text=element_text(size=12),
                                                                   axis.title=element_text(size=14))
  
  plot4 <- ggplot(df, aes(x=idx)) + 
    geom_line(aes(y = km_pf_tdxd_model[1:m]), color = "red") + 
    geom_line(aes(y = km_pf_tdxd[1:m]), color="blue", linetype="twodash") +
    ylim(0, 1) +
    ylab("Probability")+
    xlab("Month")+
    ggtitle("PFS T-DxD")+
    scale_x_continuous(breaks = round(seq(0, 18, by = 3),1))+theme(axis.text=element_text(size=12),
                                                                   axis.title=element_text(size=14))
  
  print(ae_test_chemo)
  print(ae_test_tdxd)
  plot(0:(length(ae_test_chemo)-1), ae_test_chemo, type = "l")
  plot(0:(length(ae_test_chemo)-1), ae_test_tdxd, type = "l")
  
  ggarrange(plot1, plot2, plot3, plot4,
            ncol = 2, nrow = 2, common.legend = TRUE,legend="bottom") 
  
}


plot_function(opt_var)






