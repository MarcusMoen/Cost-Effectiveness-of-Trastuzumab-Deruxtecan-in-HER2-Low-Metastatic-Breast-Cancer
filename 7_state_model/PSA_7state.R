# Please reference the manuscript for mathematical notations and variable descriptions.
CorrelateUtils <- function(U, Q, epsilon, delta){
  n <- nrow(U) #number of PSA samples
  s <- ncol(U) #number of states
  R <- matrix(rnorm(n*s,0,1), n, s) #the reference matrix.
  
  C <- matrix(0,s,s) #a place holder for the correlation matrix
  Viol <- matrix(0,s,s) #violations matrix.
  for (j in 2:s){ # {j,k} is the selected pair of state comparisons
    for (k in 1:(j-1)){
      rho <- 1 # #bivariate correlations
      X = U[,c(j,k)] #selected columns of U
      Y = R[,c(j,k)] #selected columns of R
      viol <- 0
      while(viol<epsilon & rho>=0){ #if these conditions are met, continue
        rho <- rho - delta #reduce correltion
        Xstar = induceRankCorrelation(X, Y, rho) #correlated utilities.
        viol = mean((Q[j,k] * Xstar[,1]) < (Q[j,k] * Xstar[,2])) #compute %violations between the col. vectors.
      }
      #Viol[j,k] <- viol
      C[j,k] <- rho + delta #record the desired correlation.
    }
    #print(j) #just to show the column indices.
    #print(k)
  }
  #Fill in the other elements of C.
  C = C + t(C)
  for (j in 1:s){
    C[j,j] <- 1 # % the diagonal of ones.
  }
  ## Eigenvectors and Eigenvalues correction of C
  eigenResults <- eigen(C)
  B <- eigenResults$values
  V <- eigenResults$vectors
  B[B<=0] <- 0.0001 #to make sure C is positive definite, set eigenvalues<=0 to a very small positive number
  Cstar <- V %*% diag(B) %*% solve(V) #reconstruct C
  Ustar <- induceRankCorrelation(U, R, Cstar) #similar to above, induce the correlation.
  return(Ustar)
}
## To induce Rank correlation: inputs X: QoL vectors, Y is the reference vectors, and Sigma is the correlation matrix.
induceRankCorrelation <- function(X, Y, Sigma){
  if (length(Sigma)==1){ #if Sigma is a single value, convert it to a 2x2 matrix.
    Sigma <- matrix(c(1, Sigma,
                      Sigma, 1), 2, 2)
  }
  n <- nrow(X)
  s <- ncol(X)
  #Initialize matrices.
  Xsorted <- matrix(0, n, s)
  Yrank <- matrix(0, n, s)
  Xstar <- matrix(0, n, s)
  
  P <- chol(Sigma) #compute the upper triangular matrix
  Ystar <- Y %*% P #Sort the values in the reference vectors by multiplying by P
  cor(Ystar)
  for (j in 1:s){
    Xsorted[,j] <- sort(X[,j]) #Sort each variable
    Yrank[order(Ystar[,j]),j] <- seq(1:n) #Reverse sort
    Xstar[,j]=Xsorted[Yrank[,j],j] #sort Xsorted to have the same ranks as Ystar.
  }
  return(Xstar) #return the sorted vectors.
} 


calc_alpha_beta_gamma <- function(mean, p_of_mean){
  return(c(1/p_of_mean^2, mean*p_of_mean^2))
}




source("finding_transition_probabilities_7state.R")
source("CEA_7state.R")



generate_psa_qaly_params <- function(n_sim, seed = 0){
  #Draw values for QALYs
  qaly_pf_chemo <- rbeta(n_sim, shape1 = 2, shape2 = 2)*(0.81-0.52) + rep(0.52, n_sim)
  qaly_p_chemo <- rbeta(n_sim, shape1 = 2, shape2 = 2)*(0.65-0.4) + rep(0.4, n_sim)
  qaly_pfAE_chemo <- rbeta(n_sim, shape1 = 2, shape2 = 2)*(0.707-0.417) + rep(0.417, n_sim)
  qaly_pAE_chemo <- rbeta(n_sim, shape1 = 2, shape2 = 2)*(0.65-0.297) + rep(0.297, n_sim)
  qaly_pfILD_tdxd <- rbeta(n_sim, shape1 = 2, shape2 = 2)*(0.707-0.417) + rep(0.417, n_sim)
  qaly_pILD_tdxd <- rbeta(n_sim, shape1 = 2, shape2 = 2)*(0.65-0.297) + rep(0.297, n_sim)
  
  U <- matrix(0, n_sim, 6,
              dimnames = list(seq(1,n_sim,1), c("pf", "pfAE", "pfILD", "p", "pAE", "pILD")))
  
  U[,1] <- qaly_pf_chemo
  U[,2] <- qaly_pfAE_chemo
  U[,3] <- qaly_pfILD_tdxd
  U[,4] <- qaly_p_chemo
  U[,5] <- qaly_pAE_chemo
  U[,6] <- qaly_pILD_tdxd
  
  Q <- matrix(c( 0, 0, 0, 0, 0, 0, #preference order matrix.
                 -1, 0, 0, 0, 0, 0,
                 -1, 0, 0, 0, 0, 0,
                 -1, 0, 0, 0, 0, 0,
                 -1, -1, 0, -1, 0, 0, 
                 -1, 0, -1, -1, 0, 0
  ), 6, 6, byrow = TRUE) 
  
  Ustar <- CorrelateUtils(U, Q, 0.05, 0.1)
  
  qaly_pf_chemo <- c()
  qaly_pfAE_chemo <- c()
  qaly_pfILD_tdxd <- c()
  qaly_p_chemo <- c()
  qaly_pAE_chemo <- c()
  qaly_pILD_tdxd <- c()
  
  for(i in 1:n_sim){
    error <- 0
    for(j in 1:6){
      a <- Ustar[i,]*Q[j,]
      if(max(a[a<0]) > Ustar[i,j]*-1){
        error <- error +1
      }
    }
    
    if(error == 0){
      qaly_pf_chemo <- append(qaly_pf_chemo, Ustar[i,1])
      qaly_pfAE_chemo <- append(qaly_pfAE_chemo, Ustar[i,2])
      qaly_pfILD_tdxd <- append(qaly_pfILD_tdxd, Ustar[i,3])
      qaly_p_chemo <- append(qaly_p_chemo, Ustar[i,4])
      qaly_pAE_chemo <- append(qaly_pAE_chemo, Ustar[i,5])
      qaly_pILD_tdxd <- append(qaly_pILD_tdxd, Ustar[i,6])
    }
  }
  
  
  df_qaly_params <- data.frame(
    
    qaly_pf_chemo = qaly_pf_chemo,
    qaly_p_chemo = qaly_p_chemo,
    qaly_pfAE_chemo = qaly_pfAE_chemo,
    qaly_pAE_chemo = qaly_pAE_chemo,
    qaly_pf_tdxd = qaly_pf_chemo,
    qaly_p_tdxd = qaly_p_chemo,
    qaly_pfAE_tdxd = qaly_pfAE_chemo,
    qaly_pAE_tdxd = qaly_pAE_chemo,
    qaly_pfILD_tdxd = qaly_pfILD_tdxd,
    qaly_pILD_tdxd = qaly_pILD_tdxd
  )
  
  return(df_qaly_params)
  
}



#Params: all costs, all QALYs, 
generate_psa_cost_params <- function(n_sim, seed = 0){
  set.seed(seed)

  #Calculate the cost of PF in chemo:
  price_cape <- ((rbeta(n_sim, shape1 = 2, shape2 = 2)*((801.2+200)-(801.2-200)) + rep((801.2-200), n_sim))*28*2325/9000)*0.201
  
  price_eribulin <- ((rbeta(n_sim, shape1 = 2, shape2 = 2)*((1224.17+300)-(1224.17-300)) + rep((1224.17-300), n_sim))*2*2.6/1)
  
  price_abraxane <- ((rbeta(n_sim, shape1 = 2, shape2 = 8)*((935.81+900)-935.81) + rep(935.81, n_sim))*483.6/100)*0.103
  
  price_paclitaxel <- ((rbeta(n_sim, shape1 = 8, shape2 = 2)*((155+20)-(155-95)) + rep((155-95), n_sim))*325.5/300)*0.082
  
  price_gemcitabine <- ((rbeta(n_sim, shape1 = 2, shape2 = 8)*((33.81+100)-(33.81-15)) + rep((33.81-15), n_sim))*2*2325/2000)*0.103
  
  admin_price <- 303.86*0.201+831.91*0.511+454.14*0.103+454.14*0.082+831.91*0.103
  cost_pf_chemo <- (price_cape + price_eribulin*0.511 + price_abraxane + price_paclitaxel + price_gemcitabine)*4/3 + rep(admin_price, n_sim)
  
  cost_pf_eribulin <- price_eribulin*4/3 + rep(831.91, n_sim)
  
  
  #Calculate the cost of PF in tdxd: 
  drug_price_tdxd <- (rbeta(n_sim, shape1 = 2, shape2 = 2)*((2435.71+600)-(2435.71-600)) + rep((2435.71-600), n_sim))*418.5/100
  cost_pf_tdxd <- (drug_price_tdxd)*4/3 + rep(522.43, n_sim)

  
  #Calculate the cost of P in chemo and tdxd:
  cost_p_chemo <- rgamma(n_sim, shape = (10882.601)^2/(2114.74)^2, scale = ((2114.74)^2/10882.601))
  cost_p_tdxd <- cost_p_chemo
  
  #Calculate the cost of the AE and ILD states for both chemo and tdxd:
  cost_pfAE_chemo <- rbeta(n_sim, shape1 = 2, shape2 = 2)*((4518.18*1.25)-(518.18*0.75)) + rep(518.18, n_sim) + rep(575.19, n_sim)
  cost_pAE_chemo <- cost_p_chemo
  cost_pfAE_tdxd <- cost_pfAE_chemo
  cost_pAE_tdxd <- cost_p_chemo
  cost_pfILD_tdxd <-cost_pfAE_chemo
  cost_pILD_tdxd <-cost_p_chemo
  
  
  df_cost_params <- data.frame(
    
    cost_pf_tdxd = cost_pf_tdxd,
    cost_pf_chemo = cost_pf_chemo,
    cost_p_tdxd = cost_p_tdxd,
    cost_p_chemo = cost_p_chemo,
    cost_pfAE_tdxd = cost_pfAE_tdxd,
    cost_pfAE_chemo = cost_pfAE_chemo,
    cost_pAE_tdxd = cost_pAE_tdxd,
    cost_pAE_chemo = cost_pAE_chemo,
    cost_pfILD_tdxd = cost_pfILD_tdxd,
    cost_pILD_tdxd = cost_pILD_tdxd,
    cost_pf_eribulin = cost_pf_eribulin
  )
  
  return(df_cost_params)
}







PSA_params <- function(n_sim, seed = 0){
  df_qaly_params <- generate_psa_qaly_params(n_sim)
  n_sim2 <- length(df_qaly_params$qaly_pf_chemo)
  df_cost_params <- generate_psa_cost_params(n_sim2)
  
  df_params <- cbind(df_qaly_params, df_cost_params)
  
  return(df_params)
}


draw_log_normal <- function(mu, sigma){
  location <- log(mu^2 / sqrt(sigma^2 + mu^2))
  shape <- sqrt(log(1 + (sigma^2 / mu^2)))
  draw <- rlnorm(n=1, location, shape)
  return(draw)
}



run_PSA <- function(df_params){
  
  median_os_chemo_v <- c()
  median_pf_chemo_v <- c()
  median_os_tdxd_v <- c()
  median_pf_tdxd_v <- c()
  target_ae_chemo_v <- c()
  target_ae_tdxd_v <- c()
  target_ild_tdxd_v <- c()
  
  n_sims <- length(df_params$qaly_pf_chemo)
  
  df_psa_chemo <- data.frame(
    DiscountedCost = 1:n_sims,
    DiscountedQALY = 1:n_sims
  )
  
  df_psa_tdxd <- data.frame(
    DiscountedCost = 1:n_sims,
    DiscountedQALY = 1:n_sims
  )
  
  df_psa_res <- data.frame(
    sim = 1:(n_sims*2),
    DiscountedCost = 1:(n_sims*2),
    DiscountedQALY = 1:(n_sims*2),
    group = 1:(n_sims*2),
    ICER = 1:(n_sims*2)
  )
  
  df_psa_eribulin <- data.frame(
    sim = 1:(n_sims*2),
    DiscountedCost = 1:(n_sims*2),
    DiscountedQALY = 1:(n_sims*2),
    group = 1:(n_sims*2),
    ICER = 1:(n_sims*2)
  )
  
  df_psa_cheap_tdxd <- data.frame(
    sim = 1:(n_sims*2),
    DiscountedCost = 1:(n_sims*2),
    DiscountedQALY = 1:(n_sims*2),
    group = 1:(n_sims*2),
    ICER = 1:(n_sims*2)
  )
  
  df_psa_cheap_tdxd_eribulin <- data.frame(
    sim = 1:(n_sims*2),
    DiscountedCost = 1:(n_sims*2),
    DiscountedQALY = 1:(n_sims*2),
    group = 1:(n_sims*2),
    ICER = 1:(n_sims*2)
    
  )
  
  n_cycles <- 120
  
  for(i in 1:n_sims){
    print(i)
    
    
    median_os_chemo <- draw_log_normal(mu = 16.8, sigma = 1.4)
    median_pf_chemo <- draw_log_normal(mu = 5.1, sigma = 0.66)
    median_os_tdxd <- draw_log_normal(mu = 23.4, sigma = 1.22)
    median_pf_tdxd <- draw_log_normal(mu = 9.9, sigma = 0.59)
    
    target_ae_chemo <- draw_log_normal(mu = 0.081, sigma = 0.0203)
    target_ae_tdxd <- draw_log_normal(mu = 0.0757, sigma = 0.0189)
    target_ild_tdxd <- draw_log_normal(mu = 0.0863, sigma = 0.0216)
    
    median_os_chemo_v <- append(median_os_chemo_v, median_os_chemo)
    median_pf_chemo_v <- append(median_pf_chemo_v, median_pf_chemo)
    median_os_tdxd_v <- append(median_os_tdxd_v, median_os_tdxd)
    median_pf_tdxd_v <- append(median_pf_tdxd_v, median_pf_tdxd)
    target_ae_chemo_v <- append(target_ae_chemo_v, target_ae_chemo)
    target_ae_tdxd_v <- append(target_ae_tdxd_v, target_ae_tdxd)
    target_ild_tdxd_v <- append(target_ild_tdxd_v, target_ild_tdxd)
    
    opt_var <- find_opt_var(target_var2 = c(median_os_chemo, median_pf_chemo, median_os_tdxd, median_pf_tdxd, target_ae_chemo, target_ae_tdxd, target_ild_tdxd))
    
    tm <- transition_matrices(opt_var, calibrate = F)
    A_chemo_list <- tm[[1]]
    A_tdxd_list <- tm[[2]]
    
    #print(length(A_chemo_list))
    
    #Finding the values for chemo
    df_list_chemo <- cea7state(A_chemo_list, n_cycles)
    df_chemo <- df_list_chemo[[1]]
    
    #Finding the values for tdxd
    df_list_tdxd <- cea7state(A_tdxd_list, n_cycles)
    df_tdxd <- df_list_tdxd[[1]]
    
    dr_v <- df_list_chemo[[3]]
    
    
    
    
    #The cost values
    cost_pf_chemo <- df_params[i,]$cost_pf_chemo
    cost_p_chemo <- df_params[i,]$cost_p_chemo
    cost_pfAE_chemo <- df_params[i,]$cost_pfAE_chemo
    cost_pAE_chemo <- df_params[i,]$cost_pAE_chemo
    
    cost_pf_eribulin <- df_params[i,]$cost_pf_eribulin
    
    
    cost_pf_tdxd <- df_params[i,]$cost_pf_tdxd
    cost_p_tdxd <- df_params[i,]$cost_p_tdxd
    cost_pfAE_tdxd <- df_params[i,]$cost_pfAE_tdxd
    cost_pAE_tdxd <- df_params[i,]$cost_pAE_tdxd
    
    cost_pfILD <- df_params[i,]$cost_pfILD_tdxd
    cost_pILD <- df_params[i,]$cost_pILD_tdxd
    
    
    #Calculating the cost
    
    cost_chemo_d <- (sum(df_chemo$ProgressionFree*dr_v)*cost_pf_chemo + sum(df_chemo$Progressed*dr_v)*cost_p_chemo +
                       sum(df_chemo$ProgressionFreeAE*dr_v)*cost_pfAE_chemo + sum(df_chemo$ProgressedAE*dr_v)*cost_pAE_chemo
                     - (df_chemo[1,]$ProgressionFree/2)*cost_pf_chemo) + 10375.8444
    
    cost_tdxd_d <- (sum(df_tdxd$ProgressionFree*dr_v)*cost_pf_tdxd + sum(df_chemo$Progressed*dr_v)*cost_p_tdxd +
                      sum(df_tdxd$ProgressionFreeAE*dr_v)*cost_pfAE_tdxd + sum(df_chemo$ProgressedAE*dr_v)*cost_pAE_tdxd +
                      sum(df_tdxd$ProgressionFreeILD*dr_v)*cost_pfILD + sum(df_chemo$ProgressedILD*dr_v)*cost_pILD 
                    - (df_tdxd[1,]$ProgressionFree/2)*cost_pf_tdxd)+7270.62446
    
    cost_eribulin_d <- (sum(df_chemo$ProgressionFree*dr_v)*cost_pf_eribulin + sum(df_chemo$Progressed*dr_v)*cost_p_chemo +
                       sum(df_chemo$ProgressionFreeAE*dr_v)*cost_pf_eribulin + sum(df_chemo$ProgressedAE*dr_v)*cost_pAE_chemo
                     - (df_chemo[1,]$ProgressionFree/2)*cost_pf_eribulin) + 10375.8444
    
    cost_tdxd_eribulin_d <- (sum(df_tdxd$ProgressionFree*dr_v)*cost_pf_eribulin + sum(df_chemo$Progressed*dr_v)*cost_p_tdxd +
                      sum(df_tdxd$ProgressionFreeAE*dr_v)*cost_pf_eribulin + sum(df_chemo$ProgressedAE*dr_v)*cost_pAE_tdxd +
                      sum(df_tdxd$ProgressionFreeILD*dr_v)*cost_pf_eribulin + sum(df_chemo$ProgressedILD*dr_v)*cost_pILD 
                    - (df_tdxd[1,]$ProgressionFree/2)*cost_pf_eribulin)+7270.62446
    
    
    #The QALY values
    qaly_pf_chemo <- df_params[i,]$qaly_pf_chemo
    qaly_p_chemo <- df_params[i,]$qaly_p_chemo
    qaly_pfAE_chemo <- df_params[i,]$qaly_pfAE_chemo
    qaly_pAE_chemo <- df_params[i,]$qaly_pAE_chemo
    
    
    qaly_pf_tdxd <- df_params[i,]$qaly_pf_tdxd
    qaly_p_tdxd <- df_params[i,]$qaly_p_tdxd
    qaly_pfAE_tdxd <- df_params[i,]$qaly_pfAE_tdxd
    qaly_pAE_tdxd <- df_params[i,]$qaly_pAE_tdxd
    
    qaly_pfILD <- df_params[i,]$qaly_pfILD_tdxd
    qaly_pILD <- df_params[i,]$qaly_pILD_tdxd
    
    
    #Calculate the qalys
    
    qaly_chemo_d <- (sum(df_chemo$ProgressionFree*dr_v)*qaly_pf_chemo + sum(df_chemo$Progressed*dr_v)*qaly_p_chemo +
                       sum(df_chemo$ProgressionFreeAE*dr_v)*qaly_pfAE_chemo + sum(df_chemo$ProgressedAE*dr_v)*qaly_pAE_chemo
                     - (df_chemo[1,]$ProgressionFree/2)*qaly_pf_chemo)/12
    
    qaly_tdxd_d <- (sum(df_tdxd$ProgressionFree*dr_v)*qaly_pf_tdxd + sum(df_chemo$Progressed*dr_v)*qaly_p_tdxd +
                      sum(df_tdxd$ProgressionFreeAE*dr_v)*qaly_pfAE_tdxd + sum(df_chemo$ProgressedAE*dr_v)*qaly_pAE_tdxd +
                      sum(df_tdxd$ProgressionFreeILD*dr_v)*qaly_pfILD + sum(df_chemo$ProgressedILD*dr_v)*qaly_pILD 
                    - (df_tdxd[1,]$ProgressionFree/2)*qaly_pf_tdxd)/12
    
    
    
    
    
    
    
    #df_psa_chemo[i,]$DiscountedCost <- cost_chemo_d
    #df_psa_chemo[i,]$DiscountedQALY <- qaly_chemo_d
    
    #df_psa_tdxd[i,]$DiscountedCost <- cost_tdxd_d
    #df_psa_tdxd[i,]$DiscountedQALY <- qaly_tdxd_d
    
    #T-DxD vs Mixed chemo
    df_psa_res[(2*i-1),]$DiscountedCost <- cost_chemo_d
    df_psa_res[(2*i-1),]$DiscountedQALY <- qaly_chemo_d
    df_psa_res[(2*i-1),]$group <- "Mixed Chemo"
    df_psa_res[(2*i-1),]$sim <- i
    
    df_psa_res[(2*i),]$DiscountedCost <- cost_tdxd_d
    df_psa_res[(2*i),]$DiscountedQALY <- qaly_tdxd_d
    df_psa_res[(2*i),]$group <- "T-DxD"
    df_psa_res[(2*i),]$sim <- i
    
    df_psa_res[(2*i-1),]$ICER <- (cost_tdxd_d-cost_chemo_d)/(qaly_tdxd_d-qaly_chemo_d)
    df_psa_res[(2*i),]$ICER <- (cost_tdxd_d-cost_chemo_d)/(qaly_tdxd_d-qaly_chemo_d)
    
    #T-DxD vs Eribulin
    df_psa_eribulin[(2*i-1),]$DiscountedCost <- cost_eribulin_d
    df_psa_eribulin[(2*i-1),]$DiscountedQALY <- qaly_chemo_d
    df_psa_eribulin[(2*i-1),]$group <- "Eribulin"
    df_psa_eribulin[(2*i-1),]$sim <- i
    
    df_psa_eribulin[(2*i),]$DiscountedCost <- cost_tdxd_d
    df_psa_eribulin[(2*i),]$DiscountedQALY <- qaly_tdxd_d
    df_psa_eribulin[(2*i),]$group <- "T-DxD"
    df_psa_eribulin[(2*i),]$sim <- i
    
    df_psa_eribulin[(2*i-1),]$ICER <- (cost_tdxd_d-cost_eribulin_d)/(qaly_tdxd_d-qaly_chemo_d)
    df_psa_eribulin[(2*i),]$ICER <- (cost_tdxd_d-cost_eribulin_d)/(qaly_tdxd_d-qaly_chemo_d)
    
    
    #Cheap T-DxD vs Mixed Chemo
    df_psa_cheap_tdxd[(2*i-1),]$DiscountedCost <- cost_chemo_d
    df_psa_cheap_tdxd[(2*i-1),]$DiscountedQALY <- qaly_chemo_d
    df_psa_cheap_tdxd[(2*i-1),]$group <- "Mixed Chemo"
    df_psa_cheap_tdxd[(2*i-1),]$sim <- i
    
    df_psa_cheap_tdxd[(2*i),]$DiscountedCost <- cost_tdxd_eribulin_d
    df_psa_cheap_tdxd[(2*i),]$DiscountedQALY <- qaly_tdxd_d
    df_psa_cheap_tdxd[(2*i),]$group <- "Cheap T-DxD"
    df_psa_cheap_tdxd[(2*i),]$sim <- i
    
    df_psa_cheap_tdxd[(2*i-1),]$ICER <- (cost_tdxd_eribulin_d-cost_chemo_d)/(qaly_tdxd_d-qaly_chemo_d)
    df_psa_cheap_tdxd[(2*i),]$ICER <- (cost_tdxd_eribulin_d-cost_chemo_d)/(qaly_tdxd_d-qaly_chemo_d)
    
    
    #Cheap T-DxD vs Eribulin
    df_psa_cheap_tdxd_eribulin[(2*i-1),]$DiscountedCost <- cost_eribulin_d
    df_psa_cheap_tdxd_eribulin[(2*i-1),]$DiscountedQALY <- qaly_chemo_d
    df_psa_cheap_tdxd_eribulin[(2*i-1),]$group <- "Eribulin"
    df_psa_cheap_tdxd_eribulin[(2*i-1),]$sim <- i
    
    df_psa_cheap_tdxd_eribulin[(2*i),]$DiscountedCost <- cost_tdxd_eribulin_d
    df_psa_cheap_tdxd_eribulin[(2*i),]$DiscountedQALY <- qaly_tdxd_d
    df_psa_cheap_tdxd_eribulin[(2*i),]$group <- "Cheap T-DxD"
    df_psa_cheap_tdxd_eribulin[(2*i),]$sim <- i
    
    df_psa_cheap_tdxd_eribulin[(2*i-1),]$ICER <- (cost_tdxd_eribulin_d-cost_eribulin_d)/(qaly_tdxd_d-qaly_chemo_d)
    df_psa_cheap_tdxd_eribulin[(2*i),]$ICER <- (cost_tdxd_eribulin_d-cost_eribulin_d)/(qaly_tdxd_d-qaly_chemo_d)
    
    
  }
  
  df_params_full <- df_params
  df_params_full["median_os_chemo"] <- median_os_chemo_v
  df_params_full["median_pf_chemo"] <- median_pf_chemo_v
  df_params_full["median_os_tdxd"] <- median_os_tdxd_v
  df_params_full["median_pf_tdxd"] <- median_pf_tdxd_v
  df_params_full["target_ae_chemo"] <- target_ae_chemo_v
  df_params_full["target_ae_tdxd"] <- target_ae_tdxd_v
  df_params_full["target_ild_tdxd"] <- target_ild_tdxd_v
  


  return(list(df_psa_res, df_psa_eribulin, df_psa_cheap_tdxd, df_psa_cheap_tdxd_eribulin, df_params_full))
}








n = 50
df_param <- PSA_params(n)


print("Number of simulations:")
print(nrow(df_param))

res <- run_PSA(df_param)





plot_psa_scatter <- function(df_psa_res, group_name1, group_name2){
  X<-split(df_psa_res, df_psa_res$group)
  means <- data.frame(mean_cost = c(mean(X[[group_name1]]["DiscountedCost"][[1]]), mean(X[[group_name2]]["DiscountedCost"][[1]])),
                      mean_qaly = c(mean(X[[group_name1]]["DiscountedQALY"][[1]]), mean(X[[group_name2]]["DiscountedQALY"][[1]])),
                      group = c(group_name1, group_name2))
  
  
  ggplot(df_psa_res, aes(x=DiscountedCost, y=DiscountedQALY, col = group)) + 
    geom_point() + geom_point(data=means,  mapping=aes(x = mean_cost, y = mean_qaly, size = 1)) +
    labs(x = "Discounted Cost", y = "Discounted QALY") 
}

 

plot_psa_ceac <- function(df_psa_res, group_name1, group_name2){ 
  
  df_ICER = df_psa_res[seq(1, nrow(df_psa_res), 2), ][["ICER"]]
  n <- length(df_ICER)
  
  start <- floor(min(df_ICER)/10000)*10000 - 20000
  end <- ceiling(max(df_ICER)/10000)*10000 + 20000
  m <- ceiling((end-start)/10000)
  
  
  df_ceac <- data.frame(
    x = 1:(m*2),
    y = 1:(m*2),
    group = 1:(m*2)
  )
  
  for(i in 1:m){
    wtp <- start + (i-1)*10000
    df_ceac[(2*i-1),]$x <- wtp
    df_ceac[(2*i),]$x <- wtp
    
    df_ceac[(2*i-1),]$y <- length(df_ICER[df_ICER>wtp])/n
    df_ceac[(2*i),]$y <- length(df_ICER[df_ICER<=wtp])/n
    
    df_ceac[(2*i-1),]$group <- group_name1
    df_ceac[(2*i),]$group <- group_name2
  }
  ggplot(df_ceac, aes(x=x, y=y, col = group)) + 
    geom_line()  +
    labs(x = "WTP Threshold", y = "Probability of cost-effective") 
  
}




plot_psa_scatter(res[[1]], "Mixed Chemo","T-DxD")
plot_psa_ceac(res[[1]], "Mixed Chemo","T-DxD")

plot_psa_scatter(res[[2]], "Eribulin","T-DxD")
plot_psa_ceac(res[[2]], "Eribulin","T-DxD")

plot_psa_scatter(res[[3]], "Mixed Chemo","Cheap T-DxD")
plot_psa_ceac(res[[3]] , "Mixed Chemo","Cheap T-DxD")

plot_psa_scatter(res[[4]], "Eribulin","Cheap T-DxD")
plot_psa_ceac(res[[4]], "Eribulin","Cheap T-DxD")







