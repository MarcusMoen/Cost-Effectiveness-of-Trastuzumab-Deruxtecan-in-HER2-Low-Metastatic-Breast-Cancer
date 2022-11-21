


cea <- function(A, n_cycles){
  
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
  
  #The Markov simulation
  for (t in 0:(n_cycles - 1)) {
    
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
  
  
  return(list(df, df_plot))
}



#Import the transition matrices from finding_transition_probabilities.R
source("finding_transition_probabilities.R")
A_chemo <- tm[[1]]
A_tdxd <- tm[[2]]



#Finding the values for chemo
df_list_chemo <- cea(A_chemo, 60)
df_chemo <- df_list_chemo[[1]]
df_plot_chemo <- df_list_chemo[[2]]

#Finding the values for tdxd
df_list_tdxd <- cea(A_tdxd, 60)
df_tdxd <- df_list_tdxd[[1]]
df_plot_tdxd <- df_list_tdxd[[2]]



#Plot results
library(ggplot2)
ggplot(df_plot_chemo, aes(x = cycle,y = value, group = state, color = state)) +
  geom_line()

ggplot(df_plot_tdxd, aes(x = cycle,y = value, group = state, color = state)) +
  geom_line()
























print("The LY's:")
print((sum(df[2:601,]$healthy) + sum(df[2:601,]$acuteEXD) + df[1,]$healthy/2 )/12)

print("The QALY's:")
print((sum(df$healthy)*1 + sum(df$acuteEXD)*0.7 - df[1,]$healthy/2)/12)

print("The cost:")
print((sum(df$healthy)*1500 + sum(df$acuteEXD)*3000 - (df[1,]$healthy/2)*1500)/12)

print("The discounted LY's:")
print((sum(df$healthy*dr_v) + sum(df$acuteEXD*dr_v) - df[1,]$healthy/2 )/12)

print("The discounted QALY's:")
print((sum(df$healthy*dr_v) + sum(df$acuteEXD*dr_v)*0.7 - df[1,]$healthy/2)/12)

print("The discounted costs:")
print((sum(df$healthy*dr_v)*1500 + sum(df$acuteEXD*dr_v)*3000 - (df[1,]$healthy/2)*1500)/12)

print("1b:")
print("Prevelence at age 55:")
print(df[60,]$acuteEXD/(df[60,]$acuteEXD+df[60,]$healthy))