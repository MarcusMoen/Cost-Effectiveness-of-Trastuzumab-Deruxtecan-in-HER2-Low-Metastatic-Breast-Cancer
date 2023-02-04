
#Question regarding the code:
#Should we extract the same length for each KM-curve? Or should we add all the data we can?


#NOTE:
#The transitions from rate to probability is found in https://journals.sagepub.com/doi/pdf/10.1177/0272989X05282637 Figure 3


#Function for calculating the standard error of a survival curve given 
#number at risk, n_risk, number of deaths at each timestep, n_deaths, and the survival probabilities, surv_prob

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





#Input:
#input_var = a list of three values, the rate from PF to P, the rate from P to D, and the HR.
#Output:
#The square error when using the given inputs.
eval_variables <- function(input_var){
  
  r_pf2p_chemo <- input_var[1] 
  r_p2d_chemo <- input_var[2] 
  hr <- input_var[3]
  
  #Rate of dying from other causes (background mortality)
  r_d_oc <- 0.0004
  
  #Turning chemo rates to probabilities
  lambda1 <- r_pf2p_chemo + r_d_oc
  p_pf2p_chemo <- r_pf2p_chemo*exp(-r_p2d_chemo)*(1-exp(-(lambda1-r_p2d_chemo)))/(r_pf2p_chemo+r_d_oc-r_p2d_chemo)
  p_pf2d_chemo <- (1-exp(-lambda1)) -p_pf2p_chemo
  p_p2d_chemo <- 1-exp(-(r_p2d_chemo + r_d_oc))
  
  #Turning rates and HR into probabilities for T-DxD
  r_pf2p_tdxd <- hr*r_pf2p_chemo
  lambda2 <- r_pf2p_tdxd + r_d_oc
  p_pf2p_tdxd <- r_pf2p_tdxd*exp(-r_p2d_chemo)*(1-exp(-(lambda2-r_p2d_chemo)))/(r_pf2p_tdxd+r_d_oc-r_p2d_chemo)
  p_pf2d_tdxd <- 1-exp(-lambda2)-p_pf2p_tdxd
  p_p2d_tdxd <- p_p2d_chemo
  
  
  #Creating the transition matrices
  A_chemo <- matrix(c(exp(-lambda1), 0, 0, 
                      p_pf2p_chemo,  1-p_p2d_chemo, 0, 
                      p_pf2d_chemo, p_p2d_chemo, 1), 
                    nrow = 3)
  A_tdxd <- matrix(c(exp(-lambda2), 0, 0, 
                     p_pf2p_tdxd,  1-p_p2d_tdxd, 0, 
                     p_pf2d_tdxd, p_p2d_tdxd, 1), 
                   nrow = 3)
  
  
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
  
  
  n = max(length(km_os_chemo), length(km_pf_chemo), length(km_os_tdxd), length(km_pf_tdxd)) #Find the maximum length
  
  #Start states
  s0_chemo <- c(1,0,0)
  s0_tdxd <- c(1,0,0)
  
  survival_chemo <- s0_chemo[1]+s0_chemo[2]
  km_os_chemo_model <- append(km_os_chemo_model, survival_chemo)
  km_pf_chemo_model <- append(km_pf_chemo_model, s0_chemo[1])
  
  survival_tdxd <- s0_tdxd[1]+s0_tdxd[2]
  km_os_tdxd_model <- append(km_os_tdxd_model, survival_tdxd)
  km_pf_tdxd_model <- append(km_pf_tdxd_model, s0_tdxd[1])
  
  
  #Do the simulations
  for(t in 1:n){
    s1_chemo <- s0_chemo %*% A_chemo
    survival_chemo <- s1_chemo[1]+s1_chemo[2]
    km_os_chemo_model <- append(km_os_chemo_model, survival_chemo)
    km_pf_chemo_model <- append(km_pf_chemo_model, s1_chemo[1])
    s0_chemo <- s1_chemo
    
    
    s1_tdxd <- s0_tdxd %*% A_tdxd
    survival_tdxd <- s1_tdxd[1]+s1_tdxd[2]
    km_os_tdxd_model <- append(km_os_tdxd_model, survival_tdxd)
    km_pf_tdxd_model <- append(km_pf_tdxd_model, s1_tdxd[1])
    s0_tdxd <- s1_tdxd
    
  }
  
  
  
  #Calculate the error
  serror <- sum((km_os_chemo[2:(length(km_os_chemo))] - km_os_chemo_model[2:(length(km_os_chemo))])^2) + 
    sum((km_pf_chemo[2:(length(km_pf_chemo))] - km_pf_chemo_model[2:(length(km_pf_chemo))])^2) + 
    sum((km_os_tdxd[2:(length(km_os_tdxd))] - km_os_tdxd_model[2:(length(km_os_tdxd))])^2) + 
    sum((km_pf_tdxd[2:(length(km_pf_tdxd))] - km_pf_tdxd_model[2:(length(km_pf_tdxd))])^2)
  
  
  return(serror)
  
  
  
}

#Find the optimal values (Here we could do a wider search for starting values since we might have missed the global optimum.)
fit_out <- optim(c(0.15, 0.08, 0.5), 
                 eval_variables,
                 hessian = T)

#Extract the optimal values for our three varaibles
opt_var <- fit_out$par

opt_var




#Function for finding the transition matrices such that they can be imported to the CEA file
#Takes as its input the optimal values for the eval_variables function and outputs a list of transition matrices
transition_matrices <- function(input_var){
  
  r_pf2p_chemo <- input_var[1] 
  r_p2d_chemo <- input_var[2] 
  hr <- input_var[3]
  
  #Rate of dying from other causes 
  r_d_oc <- 0.0004 #This value is not correct yet.
  
  #Turning chemo rates to probabilities
  lambda1 <- r_pf2p_chemo + r_d_oc
  p_pf2p_chemo <- r_pf2p_chemo*exp(-r_p2d_chemo)*(1-exp(-(lambda1-r_p2d_chemo)))/(r_pf2p_chemo+r_d_oc-r_p2d_chemo)
  p_pf2d_chemo <- (1-exp(-lambda1)) -p_pf2p_chemo
  p_p2d_chemo <- 1-exp(-(r_p2d_chemo + r_d_oc))
  
  #Turning rates and HR into probabilities for T-DxD
  r_pf2p_tdxd <- hr*r_pf2p_chemo
  lambda2 <- r_pf2p_tdxd + r_d_oc
  p_pf2p_tdxd <- r_pf2p_tdxd*exp(-r_p2d_chemo)*(1-exp(-(lambda2-r_p2d_chemo)))/(r_pf2p_tdxd+r_d_oc-r_p2d_chemo)
  p_pf2d_tdxd <- 1-exp(-lambda2)-p_pf2p_tdxd
  p_p2d_tdxd <- p_p2d_chemo
  
  
  #Creating the transition matrices
  A_chemo <- matrix(c(exp(-lambda1), 0, 0, 
                      p_pf2p_chemo,  1-p_p2d_chemo, 0, 
                      p_pf2d_chemo, p_p2d_chemo, 1), 
                    nrow = 3)
  A_tdxd <- matrix(c(exp(-lambda2), 0, 0, 
                     p_pf2p_tdxd,  1-p_p2d_tdxd, 0, 
                     p_pf2d_tdxd, p_p2d_tdxd, 1), 
                   nrow = 3)
  
  matrices <- list(A_chemo, A_tdxd)
  return(matrices)
}


#Find the optimal transition matrices
tm <- transition_matrices(opt_var)







#Plot the correct KM-curves vs the cureves from our model. This part of the file is used for validation
library(ggplot2)
library(ggpubr)



#The transitions from rate to probability is found in https://journals.sagepub.com/doi/pdf/10.1177/0272989X05282637 Figure 3
plot_function <- function(input_var){
  
  r_pf2p_chemo <- input_var[1] 
  r_p2d_chemo <- input_var[2] 
  hr <- input_var[3]
  
  #Rate of dying from other causes 
  r_d_oc <- 0.0004 #This value is not correct yet.
  
  #Turning chemo rates to probabilities
  lambda1 <- r_pf2p_chemo + r_d_oc
  p_pf2p_chemo <- r_pf2p_chemo*exp(-r_p2d_chemo)*(1-exp(-(lambda1-r_p2d_chemo)))/(r_pf2p_chemo+r_d_oc-r_p2d_chemo)
  p_pf2d_chemo <- (1-exp(-lambda1)) -p_pf2p_chemo
  p_p2d_chemo <- 1-exp(-(r_p2d_chemo + r_d_oc))
  
  #Turning rates and HR into probabilities for T-DxD
  r_pf2p_tdxd <- hr*r_pf2p_chemo
  lambda2 <- r_pf2p_tdxd + r_d_oc
  p_pf2p_tdxd <- r_pf2p_tdxd*exp(-r_p2d_chemo)*(1-exp(-(lambda2-r_p2d_chemo)))/(r_pf2p_tdxd+r_d_oc-r_p2d_chemo)
  p_pf2d_tdxd <- 1-exp(-lambda2)-p_pf2p_tdxd
  p_p2d_tdxd <- p_p2d_chemo
  
  
  #Creating the transition matrices
  A_chemo <- matrix(c(exp(-lambda1), 0, 0, 
                      p_pf2p_chemo,  1-p_p2d_chemo, 0, 
                      p_pf2d_chemo, p_p2d_chemo, 1), 
                    nrow = 3)
  A_tdxd <- matrix(c(exp(-lambda2), 0, 0, 
                     p_pf2p_tdxd,  1-p_p2d_tdxd, 0, 
                     p_pf2d_tdxd, p_p2d_tdxd, 1), 
                   nrow = 3)
  
  
  #The extracted Kaplan-Meier values
  km_os_chemo <- c(1, 0.986, 0.981, 0.957, 0.94, 0.927, 0.884, 0.841, 0.793, 0.738, 0.712, 0.683, 0.669, 0.634, 0.610, 0.566, 0.520, 0.484, 0.459, 0.460, 0.457, 0.419, 0.375)
  km_pf_chemo <- c(1, 0.983, 0.756, 0.62, 0.603, 0.501, 0.439, 0.398, 0.34, 0.277, 0.268, 0.226, 0.217, 0.185) 
  km_os_tdxd <- c(1, 0.992, 0.985, 0.967, 0.958, 0.939, 0.928, 0.899, 0.87, 0.857, 0.829, 0.815, 0.793, 0.774, 0.742, 0.731, 0.688, 0.661, 0.623, 0.586, 0.566, 0.539, 0.507, 0.500, 0.484)
  km_pf_tdxd <- c(1, 0.988, 0.896, 0.82, 0.809, 0.765, 0.681, 0.633, 0.592, 0.554, 0.491, 0.457, 0.422, 0.386, 0.371, 0.363, 0.339, 0.31, 0.295, 0.274, 0.259) 
  
  #Vectors for storing the estimates of the KM curves
  km_os_chemo_model <- c()
  km_pf_chemo_model <- c()
  km_os_tdxd_model <- c()
  km_pf_tdxd_model <- c() 
  
  
  n = max(length(km_os_chemo), length(km_pf_chemo), length(km_os_tdxd), length(km_pf_tdxd))
  
  #Start states
  s0_chemo <- c(1,0,0)
  s0_tdxd <- c(1,0,0)
  
  survival_chemo <- s0_chemo[1]+s0_chemo[2]
  km_os_chemo_model <- append(km_os_chemo_model, survival_chemo)
  km_pf_chemo_model <- append(km_pf_chemo_model, s0_chemo[1])
  
  survival_tdxd <- s0_tdxd[1]+s0_tdxd[2]
  km_os_tdxd_model <- append(km_os_tdxd_model, survival_tdxd)
  km_pf_tdxd_model <- append(km_pf_tdxd_model, s0_tdxd[1])
  
  
  for(t in 1:n){
    s1_chemo <- s0_chemo %*% A_chemo
    survival_chemo <- s1_chemo[1]+s1_chemo[2]
    km_os_chemo_model <- append(km_os_chemo_model, survival_chemo)
    km_pf_chemo_model <- append(km_pf_chemo_model, s1_chemo[1])
    s0_chemo <- s1_chemo
    
    
    s1_tdxd <- s0_tdxd %*% A_tdxd
    survival_tdxd <- s1_tdxd[1]+s1_tdxd[2]
    km_os_tdxd_model <- append(km_os_tdxd_model, survival_tdxd)
    km_pf_tdxd_model <- append(km_pf_tdxd_model, s1_tdxd[1])
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
  
  ggarrange(plot1, plot2, plot3, plot4,
            ncol = 2, nrow = 2, common.legend = TRUE,legend="bottom") 
}


plot_function(opt_var)






