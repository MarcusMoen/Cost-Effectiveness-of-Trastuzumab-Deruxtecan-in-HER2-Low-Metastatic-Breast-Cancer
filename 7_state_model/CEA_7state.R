
#Input:
#A = a 7x7 transition matrix
#n_cycles = number of cycles
#Function for doing n_cycles given a 7x7 transition matrix
#Output:
#A list of two DataFrames. One with the probability values for each of the three states at each cycle. And one which is used for plotting.
cea7state <- function(A_list, n_cycles){
  
  #Make dataframes to store values
  df <- data.frame(
    cycle = 0:n_cycles,
    ProgressionFree = 0:n_cycles,
    ProgressionFreeAE = 0:n_cycles,
    ProgressionFreeILD = 0:n_cycles,
    Progressed = 0:n_cycles,
    ProgressedAE = 0:n_cycles,
    ProgressedILD = 0:n_cycles,
    Dead = 0:n_cycles
  )
  
  df_plot <- data.frame(
    cycle = rep(0:n_cycles, 7),
    state = 0:(n_cycles*7+6),
    value = 0:(n_cycles*7+6)
  )
  
  #The discount rate per cycle (assuming an annual discount rate of 3%)
  dr_c <- (1+0.03)^(1/12) - 1
  #A discount rate vector (For storing the discount at each timestep)
  dr_v <- 0:n_cycles
  dr_v[1] <- 1
  
  ly <- 1:n_cycles
  
  #The initial state
  s0 <- matrix(c(1,0,0,0,0,0,0), ncol = 7)
  
  cl <- 1/12
  
  #The Markov simulation
  for (t in 1:(n_cycles)) {
    
    #Update the discount rate vector
    dr_v[t+1] <- (1/(1+0.03))^((t)*cl)
    
    #Do one cycle
    #print(A_list[t])
    s1 <- s0 %*% A_list[[t]]
    #Update the dataframes
    df["ProgressionFree"][df["ProgressionFree"] == t] <- s1[1,1]
    df["ProgressionFreeAE"][df["ProgressionFreeAE"] == t] <- s1[1,2]
    df["ProgressionFreeILD"][df["ProgressionFreeILD"] == t] <- s1[1,3]
    df["Progressed"][df["Progressed"] == t] <- s1[1,4]
    df["ProgressedAE"][df["ProgressedAE"] == t] <- s1[1,5]
    df["ProgressedILD"][df["ProgressedILD"] == t] <- s1[1,6]
    df["Dead"][df["Dead"] == t] <- s1[1,7] 
    
    df_plot["state"][df_plot["value"] == t] <- "ProgressionFree"
    df_plot["state"][df_plot["value"] == t + n_cycles+1] <- "ProgressionFreeAE"
    df_plot["state"][df_plot["value"] == t + n_cycles*2+2] <- "ProgressionFreeILD"
    df_plot["state"][df_plot["value"] == t + n_cycles*3+3] <- "Progressed"
    df_plot["state"][df_plot["value"] == t + n_cycles*4+4] <- "ProgressedAE"
    df_plot["state"][df_plot["value"] == t + n_cycles*5+5] <- "ProgressedILD"
    df_plot["state"][df_plot["value"] == t + n_cycles*6+6] <- "Dead"
    df_plot["value"][df_plot["value"] == t] <- s1[1,1]
    df_plot["value"][df_plot["value"] == t + n_cycles+1] <- s1[1,2]
    df_plot["value"][df_plot["value"] == t + n_cycles*2+2] <- s1[1,3]
    df_plot["value"][df_plot["value"] == t + n_cycles*3+3] <- s1[1,4]
    df_plot["value"][df_plot["value"] == t + n_cycles*4+4] <- s1[1,5]
    df_plot["value"][df_plot["value"] == t + n_cycles*5+5] <- s1[1,6]
    df_plot["value"][df_plot["value"] == t + n_cycles*6+6] <- s1[1,7]
    
    s0 <- s1
    
  }
  
  df["ProgressionFree"][df["ProgressionFree"] == 0] <- 1
  df_plot["value"][df_plot["state"] == 0] <- 1
  df_plot["value"][df_plot["state"] == n_cycles+1] <- 0
  df_plot["value"][df_plot["state"] == n_cycles*2+2] <- 0
  df_plot["value"][df_plot["state"] == n_cycles*3+3] <- 0
  df_plot["value"][df_plot["state"] == n_cycles*4+4] <- 0
  df_plot["value"][df_plot["state"] == n_cycles*5+5] <- 0
  df_plot["value"][df_plot["state"] == n_cycles*6+6] <- 0
  
  
  df_plot["state"][df_plot["state"] == 0] <- "ProgressionFree"
  df_plot["state"][df_plot["state"] == n_cycles+1] <- "ProgressionFreeAE"
  df_plot["state"][df_plot["state"] == n_cycles*2+2] <- "ProgressionFreeILD"
  df_plot["state"][df_plot["state"] == n_cycles*3+3] <- "Progressed"
  df_plot["state"][df_plot["state"] == n_cycles*4+4] <- "ProgressedAE"
  df_plot["state"][df_plot["state"] == n_cycles*5+5] <- "ProgressedILD"
  df_plot["state"][df_plot["state"] == n_cycles*6+6] <- "Dead"
  
  
  return(list(df, df_plot, dr_v))
}




#Import the transition matrices from finding_transition_probabilities.R
source("finding_transition_probabilities_7state.R")
A_chemo_list <- tm[[1]]
A_tdxd_list <- tm[[2]]

n_cycles <- 120

#Finding the values for chemo
df_list_chemo <- cea7state(A_chemo_list, n_cycles)
df_chemo <- df_list_chemo[[1]]
df_plot_chemo <- df_list_chemo[[2]]

#Finding the values for tdxd
df_list_tdxd <- cea7state(A_tdxd_list, n_cycles)
df_tdxd <- df_list_tdxd[[1]]
df_plot_tdxd <- df_list_tdxd[[2]]



#Plot results
library(ggplot2)
p1 <- ggplot(df_plot_chemo, aes(x = cycle,y = value, group = state, color = state)) +
  geom_line() + 
  ggtitle("Evolution of patients with physician's choice treatment") +
  xlab("Months from start") + ylab("Probability") + 
  theme(
    plot.title = element_text(color="dodgerblue4", size=14, face="bold.italic"),
    axis.title.x = element_text(color="darkolivegreen4", size=14, face="bold"),
    axis.title.y = element_text(color="#993333", size=14, face="bold")
    )

p2 <- ggplot(df_plot_tdxd, aes(x = cycle,y = value, group = state, color = state)) +
  geom_line() + 
  ggtitle("Evolution of patients with T-DxD treatment") +
  xlab("Months from start") + ylab("Probability") + 
  theme(
  plot.title = element_text(color="dodgerblue4", size=14, face="bold.italic"),
  axis.title.x = element_text(color="darkolivegreen4", size=14, face="bold"),
  axis.title.y = element_text(color="#993333", size=14, face="bold")
  )


#p1
#p2

ggarrange(p1, p2,
          #labels = c("OS chemo", "OS T-DxD", "PFS chemo", "PFS T-DxD"),
          nrow = 2)

p_pf <- ggplot(NULL, aes(x = cycle, y = ProgressionFree)) + 
  geom_point(data = df_chemo, color = "blue") +
  geom_point(data = df_tdxd, color = "red") +
  ggtitle("% Progression Free") +
  xlab("Months from start") + ylab("Probability") +
  scale_x_continuous(breaks = round(seq(0, 120, by = 12),1))+theme(axis.text=element_text(size=12),
                                                                   axis.title=element_text(size=14),
                                                                   plot.title=element_text(size=20))

p_p <- ggplot(NULL, aes(x = cycle, y = Progressed)) + 
    geom_point(data = df_chemo, color = "blue") +
    geom_point(data = df_tdxd, color = "red") +
    ggtitle("% Progressed") +
    xlab("Months from start") + ylab("Probability") +
  scale_x_continuous(breaks = round(seq(0, 120, by = 12),1))+theme(axis.text=element_text(size=12),
                                                                   axis.title=element_text(size=14),
                                                                   plot.title=element_text(size=20))

p_a <- ggplot(NULL, aes(x = cycle, y = 1-Dead)) + 
  geom_point(data = df_chemo, color = "blue") +
  geom_point(data = df_tdxd, color = "red") +
  ggtitle("% Alive") +
  xlab("Months from start") + ylab("Probability") +
  scale_x_continuous(breaks = round(seq(0, 120, by = 12),1))+theme(axis.text=element_text(size=12),
                                                                 axis.title=element_text(size=14),
                                                                 plot.title=element_text(size=20))


p_pf

p_p

p_a

ggarrange(p_pf, p_p, p_a)




#Everything under her needs to be updated!!!


#The cost values
cost_pf_chemo <- 7203.56
#cost_pf_chemo <- 11000
cost_p_chemo <- 10882.6
cost_pfAE_chemo <- 5093.37
cost_pAE_chemo <- 10882.6


cost_pf_tdxd <- 14113.69
cost_p_tdxd <- 10882.6
cost_pfAE_tdxd <- 5093.37
cost_pAE_tdxd <- 10882.6

cost_pfILD <- 5093.37
cost_pILD <- 10882.6

dr_v <- df_list_chemo[[3]]

#Calculating the cost
cost_chemo <- (sum(df_chemo$ProgressionFree)*cost_pf_chemo + sum(df_chemo$Progressed)*cost_p_chemo +
              sum(df_chemo$ProgressionFreeAE)*cost_pfAE_chemo + sum(df_chemo$ProgressedAE)*cost_pAE_chemo
               - (df_chemo[1,]$ProgressionFree/2)*cost_pf_chemo) + 10375.8444

cost_tdxd <- (sum(df_tdxd$ProgressionFree)*cost_pf_tdxd + sum(df_chemo$Progressed)*cost_p_tdxd +
                sum(df_tdxd$ProgressionFreeAE)*cost_pfAE_tdxd + sum(df_chemo$ProgressedAE)*cost_pAE_tdxd +
                sum(df_tdxd$ProgressionFreeILD)*cost_pfILD + sum(df_chemo$ProgressedILD)*cost_pILD 
              - (df_tdxd[1,]$ProgressionFree/2)*cost_pf_tdxd)+7270.62446

cost_chemo_d <- (sum(df_chemo$ProgressionFree*dr_v)*cost_pf_chemo + sum(df_chemo$Progressed*dr_v)*cost_p_chemo +
                   sum(df_chemo$ProgressionFreeAE*dr_v)*cost_pfAE_chemo + sum(df_chemo$ProgressedAE*dr_v)*cost_pAE_chemo
                 - (df_chemo[1,]$ProgressionFree/2)*cost_pf_chemo) + 10375.8444

cost_tdxd_d <- (sum(df_tdxd$ProgressionFree*dr_v)*cost_pf_tdxd + sum(df_chemo$Progressed*dr_v)*cost_p_tdxd +
                     sum(df_tdxd$ProgressionFreeAE*dr_v)*cost_pfAE_tdxd + sum(df_chemo$ProgressedAE*dr_v)*cost_pAE_tdxd +
                     sum(df_tdxd$ProgressionFreeILD*dr_v)*cost_pfILD + sum(df_chemo$ProgressedILD*dr_v)*cost_pILD 
                   - (df_tdxd[1,]$ProgressionFree/2)*cost_pf_tdxd)+7270.62446


#The QALY values
qaly_pf_chemo <- 0.65-0.0715527
qaly_p_chemo <- 0.54
qaly_pfAE_chemo <- 0.547
qaly_pAE_chemo <- 0.437


qaly_pf_tdxd <- 0.65-0.05670907
qaly_p_tdxd <- 0.54
qaly_pfAE_tdxd <- 0.547
qaly_pAE_tdxd <- 0.437
  
qaly_pfILD <- 0.547
qaly_pILD <- 0.437


#Calculate the qalys

qaly_chemo <- (sum(df_chemo$ProgressionFree)*qaly_pf_chemo + sum(df_chemo$Progressed)*qaly_p_chemo +
                 sum(df_chemo$ProgressionFreeAE)*qaly_pfAE_chemo + sum(df_chemo$ProgressedAE)*qaly_pAE_chemo
               - (df_chemo[1,]$ProgressionFree/2)*qaly_pf_chemo)/12

qaly_tdxd <- (sum(df_tdxd$ProgressionFree)*qaly_pf_tdxd + sum(df_chemo$Progressed)*qaly_p_tdxd +
                sum(df_tdxd$ProgressionFreeAE)*qaly_pfAE_tdxd + sum(df_chemo$ProgressedAE)*qaly_pAE_tdxd +
                sum(df_tdxd$ProgressionFreeILD)*qaly_pfILD + sum(df_chemo$ProgressedILD)*qaly_pILD 
              - (df_tdxd[1,]$ProgressionFree/2)*qaly_pf_tdxd)/12

qaly_chemo_d <- (sum(df_chemo$ProgressionFree*dr_v)*qaly_pf_chemo + sum(df_chemo$Progressed*dr_v)*qaly_p_chemo +
                   sum(df_chemo$ProgressionFreeAE*dr_v)*qaly_pfAE_chemo + sum(df_chemo$ProgressedAE*dr_v)*qaly_pAE_chemo
                 - (df_chemo[1,]$ProgressionFree/2)*qaly_pf_chemo)/12

qaly_tdxd_d <- (sum(df_tdxd$ProgressionFree*dr_v)*qaly_pf_tdxd + sum(df_chemo$Progressed*dr_v)*qaly_p_tdxd +
                  sum(df_tdxd$ProgressionFreeAE*dr_v)*qaly_pfAE_tdxd + sum(df_chemo$ProgressedAE*dr_v)*qaly_pAE_tdxd +
                  sum(df_tdxd$ProgressionFreeILD*dr_v)*qaly_pfILD + sum(df_chemo$ProgressedILD*dr_v)*qaly_pILD 
                - (df_tdxd[1,]$ProgressionFree/2)*qaly_pf_tdxd)/12


#Calculate total life years
ly_chemo <- (sum(df_chemo$ProgressionFree) + sum(df_chemo$Progressed) +
               sum(df_chemo$ProgressionFreeAE) + sum(df_chemo$ProgressedAE)
             - df_chemo[1,]$ProgressionFree/2 )/12
ly_tdxd <- (sum(df_tdxd$ProgressionFree) + sum(df_tdxd$Progressed) +
              sum(df_tdxd$ProgressionFreeAE) + sum(df_tdxd$ProgressedAE)+
              sum(df_tdxd$ProgressionFreeILD) + sum(df_tdxd$ProgressedILD)
            - df_tdxd[1,]$ProgressionFree/2 )/12

ly_chemo_notAE <- (sum(df_chemo$ProgressionFree) + sum(df_chemo$Progressed)
             - df_chemo[1,]$ProgressionFree/2 )/12

ly_tdxd_notAE <- (sum(df_tdxd$ProgressionFree) + sum(df_tdxd$Progressed))/12

ly_chemo_AE <- (
               sum(df_chemo$ProgressionFreeAE) + sum(df_chemo$ProgressedAE))/12
ly_tdxd_AE <- (
              sum(df_tdxd$ProgressionFreeAE) + sum(df_tdxd$ProgressedAE))/12

ly_tdxd_ILD <- (
              sum(df_tdxd$ProgressionFreeILD) + sum(df_tdxd$ProgressedILD))/12

#Summary of the results
df_res <- data.frame(
  Strategy = c("Chemo:", "T-Dxd:"),
  Costs = c(cost_chemo, cost_tdxd),
  DiscountedCost = c(cost_chemo_d, cost_tdxd_d),
  LifeYears = c(ly_chemo, ly_tdxd),
  QALYs = c(qaly_chemo, qaly_tdxd),
  DiscountedQALY = c(qaly_chemo_d, qaly_tdxd_d),
  Incremental_Costs = c(0, cost_tdxd_d-cost_chemo_d),
  Incremental_QALYs = c(0, qaly_tdxd_d-qaly_chemo_d),
  ICER = c(0, (cost_tdxd_d-cost_chemo_d)/(qaly_tdxd_d-qaly_chemo_d))
)

df_res




















#ggplot(data=df_res, aes(x=Costs, y=QALYs)) +
 # geom_line(linetype = "dashed")+
  #geom_point()+ 
  #ggtitle("Cost-effectiveness plane") +
  #xlab("Cost") + ylab("QALY")
  
