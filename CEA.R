
#Input:
#A = a 3x3 transition matrix
#n_cycles = number of cycles
#Function for doing n_cycles given a 3x3 transition matrix
#Output:
#A list of two DataFrames. One with the probability values for each of the three states at each cycle. And one which is used for plotting.
cea3state <- function(A, n_cycles){
  
  #Make dataframes to store values
  df <- data.frame(
    cycle = 0:n_cycles,
    ProgressionFree = 0:n_cycles,
    Progressed = 0:n_cycles,
    Dead = 0:n_cycles
  )
  
  df_plot <- data.frame(
    cycle = rep(0:n_cycles, 3),
    state = 0:(n_cycles*3+2),
    value = 0:(n_cycles*3+2)
  )
  
  #The discount rate per cycle (assuming an annual discount rate of 3%)
  dr_c <- (1+0.03)^(1/12) - 1
  #A discount rate vector (For storing the discount at each timestep)
  dr_v <- 0:n_cycles
  dr_v[1] <- 1
  
  ly <- 1:n_cycles
  
  #The initial state
  s0 <- matrix(c(1,0,0), ncol = 3)
  
  cl <- 1/12
  
  #The Markov simulation
  for (t in 0:(n_cycles - 1)) {
    
    #Update the discount rate vector
    dr_v[t+2] <- (1/(1+0.03))^((t+1)*cl)
    
    #Do one cycle
    s1 <- s0 %*% A
    #Update the dataframes
    df["ProgressionFree"][df["ProgressionFree"] == t+1] <- s1[1,1]
    df["Progressed"][df["Progressed"] == t+1] <- s1[1,2]
    df["Dead"][df["Dead"] == t+1] <- s1[1,3] 
    
    df_plot["state"][df_plot["value"] == t+1] <- "ProgressionFree"
    df_plot["state"][df_plot["value"] == t+1 + n_cycles+1] <- "Progressed"
    df_plot["state"][df_plot["value"] == t+1 + n_cycles*2+2] <- "Dead"
    df_plot["value"][df_plot["value"] == t+1] <- s1[1,1]
    df_plot["value"][df_plot["value"] == t+1 + n_cycles+1] <- s1[1,2]
    df_plot["value"][df_plot["value"] == t+1 + n_cycles*2+2] <- s1[1,3]
    
    s0 <- s1
    
  }
  
  df["ProgressionFree"][df["ProgressionFree"] == 0] <- 1
  df_plot["value"][df_plot["state"] == 0] <- 1
  df_plot["value"][df_plot["state"] == n_cycles+1] <- 0
  df_plot["value"][df_plot["state"] == n_cycles*2+2] <- 0
  df_plot["state"][df_plot["state"] == 0] <- "ProgressionFree"
  df_plot["state"][df_plot["state"] == n_cycles+1] <- "Progressed"
  df_plot["state"][df_plot["state"] == n_cycles*2+2] <- "Dead"
  
  
  return(list(df, df_plot, dr_v))
}



#Import the transition matrices from finding_transition_probabilities.R
source("finding_transition_probabilities.R")
A_chemo <- tm[[1]]
A_tdxd <- tm[[2]]

n_cycles <- 120

#Finding the values for chemo
df_list_chemo <- cea3state(A_chemo, n_cycles)
df_chemo <- df_list_chemo[[1]]
df_plot_chemo <- df_list_chemo[[2]]

#Finding the values for tdxd
df_list_tdxd <- cea3state(A_tdxd, n_cycles)
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




#These values need fixing
cost_pf_chemo <- 7776
cost_p_chemo <- 11337.97
cost_pf_tdxd <- 14602
cost_p_tdxd <- 11337.97

dr_v <- df_list_chemo[[3]]

#Calculating the cost
cost_chemo <- (sum(df_chemo$ProgressionFree)*cost_pf_chemo + sum(df_chemo$Progressed)*cost_p_chemo - (df_chemo[1,]$ProgressionFree/2)*cost_pf_chemo)
cost_tdxd <- (sum(df_tdxd$ProgressionFree)*cost_pf_tdxd + sum(df_chemo$Progressed)*cost_p_tdxd - (df_tdxd[1,]$ProgressionFree/2)*cost_pf_tdxd)

cost_chemo_d <- (sum(df_chemo$ProgressionFree*dr_v)*cost_pf_chemo + sum(df_chemo$Progressed*dr_v)*cost_p_chemo - (df_chemo[1,]$ProgressionFree/2)*cost_pf_chemo)
cost_tdxd_d <- (sum(df_tdxd$ProgressionFree*dr_v)*cost_pf_tdxd + sum(df_chemo$Progressed*dr_v)*cost_p_tdxd - (df_tdxd[1,]$ProgressionFree/2)*cost_pf_tdxd)


#These values need fixing
qaly_pf_chemo <- 0.601
qaly_p_chemo <- 0.54
qaly_pf_tdxd <- 0.602
qaly_p_tdxd <- 0.54

#These values need fixing
#qaly_pf_chemo <- 1
#qaly_p_chemo <- 0.54
#qaly_pf_tdxd <- 1
#qaly_p_tdxd <- 0.54

#Calculate the qalys
qaly_chemo <- (sum(df_chemo$ProgressionFree)*qaly_pf_chemo + sum(df_chemo$Progressed)*qaly_p_chemo - df_chemo[1,]$ProgressionFree*qaly_pf_chemo/2 )/12
qaly_tdxd <- (sum(df_tdxd$ProgressionFree)*qaly_pf_tdxd + sum(df_tdxd$Progressed)*qaly_p_tdxd - df_tdxd[1,]$ProgressionFree*qaly_pf_tdxd/2 )/12

qaly_chemo_d <- (sum(df_chemo$ProgressionFree*dr_v)*qaly_pf_chemo + sum(df_chemo$Progressed*dr_v)*qaly_p_chemo - df_chemo[1,]$ProgressionFree*qaly_pf_chemo/2 )/12
qaly_tdxd_d <- (sum(df_tdxd$ProgressionFree*dr_v)*qaly_pf_tdxd + sum(df_tdxd$Progressed*dr_v)*qaly_p_tdxd - df_tdxd[1,]$ProgressionFree*qaly_pf_tdxd/2 )/12

#Calculate total life years
ly_chemo <- (sum(df_chemo$ProgressionFree) + sum(df_chemo$Progressed) - df_chemo[1,]$ProgressionFree/2 )/12
ly_tdxd <- (sum(df_tdxd$ProgressionFree) + sum(df_tdxd$Progressed) - df_tdxd[1,]$ProgressionFree/2 )/12


#Summary of the results
df_res <- data.frame(
  Strategy = c("Chemo:", "T-Dxd:"),
  Costs = c(cost_chemo, cost_tdxd),
  DiscountedCost = c(cost_chemo_d, cost_tdxd_d),
  LifeYears = c(ly_chemo, ly_tdxd),
  QALYs = c(qaly_chemo, qaly_tdxd),
  DiscountedQALY = c(qaly_chemo_d, qaly_tdxd_d),
  Incremental_Costs = c(0, cost_tdxd-cost_chemo),
  Incremental_QALYs = c(0, qaly_tdxd-qaly_chemo),
  ICER = c(0, (cost_tdxd_d-cost_chemo_d)/(qaly_tdxd_d-qaly_chemo_d))
)

df_res




















#ggplot(data=df_res, aes(x=Costs, y=QALYs)) +
 # geom_line(linetype = "dashed")+
  #geom_point()+ 
  #ggtitle("Cost-effectiveness plane") +
  #xlab("Cost") + ylab("QALY")
  
