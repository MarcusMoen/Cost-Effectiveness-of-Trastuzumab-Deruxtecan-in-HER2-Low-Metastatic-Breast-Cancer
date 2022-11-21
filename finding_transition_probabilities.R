
#Question regarding the code:
#Should we extract the same length for each KM-curve? Or should we add all the data we can?






#The transitions from rate to probability is found in https://journals.sagepub.com/doi/pdf/10.1177/0272989X05282637 Figure 3
finding_transition_prob <- function(input_var){
  
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
  km_os_chemo <- c(0.986, 0.981, 0.957, 0.94, 0.927, 0.884, 0.841, 0.793, 0.738, 0.712, 0.683, 0.669, 0.634, 0.610, 0.566, 0.520, 0.484, 0.459, 0.460, 0.457, 0.419, 0.375)
  km_pf_chemo <- c(0.983, 0.756, 0.62, 0.603, 0.501, 0.439, 0.398, 0.34, 0.277, 0.268, 0.226, 0.217, 0.185) 
  km_os_tdxd <- c(0.992, 0.985, 0.967, 0.958, 0.939, 0.928, 0.899, 0.87, 0.857, 0.829, 0.815, 0.793, 0.774, 0.742, 0.731, 0.688, 0.661, 0.623, 0.586, 0.566, 0.539, 0.507, 0.500, 0.484)
  km_pf_tdxd <- c(0.988, 0.896, 0.82, 0.809, 0.765, 0.681, 0.633, 0.592, 0.554, 0.491, 0.457, 0.422, 0.386, 0.371, 0.363, 0.339, 0.31, 0.295, 0.274, 0.259) 
  
  #Vectors for storing the estimates of the KM curves
  km_os_chemo_model <- c()
  km_pf_chemo_model <- c()
  km_os_tdxd_model <- c()
  km_pf_tdxd_model <- c() 
  
  
  n = max(length(km_os_chemo), length(km_pf_chemo), length(km_os_tdxd), length(km_pf_tdxd)) #Find the maximum length
  
  #Start states
  s0_chemo <- c(1,0,0)
  s0_tdxd <- c(1,0,0)
  
  
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
  error2 <- sum((km_os_chemo-km_os_chemo_model[1:length(km_os_chemo)])^2) + 
    sum((km_pf_chemo-km_pf_chemo_model[1:length(km_pf_chemo)])^2) + 
    sum((km_os_tdxd-km_os_tdxd_model[1:length(km_os_tdxd)])^2) + 
    sum((km_pf_tdxd-km_pf_tdxd_model[1:length(km_pf_tdxd)])^2)
  
  
  return(error2)
  
  
}


#Find the optimal values
fit_out <- optim(c(0.03, 0.04, 0.5), 
                 finding_transition_prob,
                 hessian = T)
opt_var <- fit_out$par

opt_var




#Function for finding the transition matrices such that they can be imported to the CEA file
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
  r_d_oc <- 0.000004 #This value is not correct yet.
  
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
  km_os_chemo <- c(0.986, 0.981, 0.957, 0.94, 0.927, 0.884, 0.841, 0.793, 0.738, 0.712, 0.683, 0.669, 0.634, 0.610, 0.566, 0.520, 0.484, 0.459, 0.460, 0.457, 0.419, 0.375)
  km_pf_chemo <- c(0.983, 0.756, 0.62, 0.603, 0.501, 0.439, 0.398, 0.34, 0.277, 0.268, 0.226, 0.217, 0.185) 
  km_os_tdxd <- c(0.992, 0.985, 0.967, 0.958, 0.939, 0.928, 0.899, 0.87, 0.857, 0.829, 0.815, 0.793, 0.774, 0.742, 0.731, 0.688, 0.661, 0.623, 0.586, 0.566, 0.539, 0.507, 0.500, 0.484)
  km_pf_tdxd <- c(0.988, 0.896, 0.82, 0.809, 0.765, 0.681, 0.633, 0.592, 0.554, 0.491, 0.457, 0.422, 0.386, 0.371, 0.363, 0.339, 0.31, 0.295, 0.274, 0.259) 
  
  #Vectors for storing the estimates of the KM curves
  km_os_chemo_model <- c()
  km_pf_chemo_model <- c()
  km_os_tdxd_model <- c()
  km_pf_tdxd_model <- c() 
  
  
  n = max(length(km_os_chemo), length(km_pf_chemo), length(km_os_tdxd), length(km_pf_tdxd))
  
  #Start states
  s0_chemo <- c(1,0,0)
  s0_tdxd <- c(1,0,0)
  
  
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
  
  #ggplot(data=df, aes(x=idx, y=km_os_chemo_model)) +
   # geom_line(linetype = "dashed")+
    #geom_point()
  
  
  
  plot1 <- ggplot(df, aes(x=idx)) + 
    geom_line(aes(y = km_os_chemo_model[1:m]), color = "darkred") + 
    geom_line(aes(y = km_os_chemo[1:m]), color="steelblue", linetype="twodash") +
    ylim(0, 1) +
    ylab("OS chemo")
  
  plot2 <- ggplot(df, aes(x=idx)) + 
    geom_line(aes(y = km_os_tdxd_model[1:m]), color = "darkred") + 
    geom_line(aes(y = km_os_tdxd[1:m]), color="steelblue", linetype="twodash") +
    ylim(0, 1) +
    ylab("OS T-DxD")
  
  plot3 <- ggplot(df, aes(x=idx)) + 
    geom_line(aes(y = km_pf_chemo_model[1:m]), color = "darkred") + 
    geom_line(aes(y = km_pf_chemo[1:m]), color="steelblue", linetype="twodash") +
    ylim(0, 1) +
    ylab("PFS chemo")
  
  plot4 <- ggplot(df, aes(x=idx)) + 
    geom_line(aes(y = km_pf_tdxd_model[1:m]), color = "darkred") + 
    geom_line(aes(y = km_pf_tdxd[1:m]), color="steelblue", linetype="twodash") +
    ylim(0, 1) +
    ylab("PFS T-DxD")
  
  ggarrange(plot1, plot2, plot3, plot4,
            #labels = c("OS chemo", "OS T-DxD", "PFS chemo", "PFS T-DxD"),
            ncol = 2, nrow = 2)
}


plot_function(opt_var)







