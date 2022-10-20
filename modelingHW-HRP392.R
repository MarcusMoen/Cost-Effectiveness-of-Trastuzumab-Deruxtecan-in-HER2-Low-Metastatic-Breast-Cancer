#The age-specific mortality rates
m <- c(0.005, 0.01, 0.02, 0.05, 0.07, 0.1, 0.13, 0.17, 0.2, 0.3, 1)

#The age-specific mortality rates
m_exd <- c(0.005, 0.01, 0.02, 0.05, 0.07, 0.1, 0.13, 0.17, 0.2, 0.3, 1)
m_exd2 <- c(0.005, 0.01, 0.02, 0.05, 0.07, 0.1, 0.13, 0.17, 0.2, 0.3, 1)

#The cycle length
cl <- 1/12
#Annual costs
ac_h = 1500
ac_exd = 3000
#Utility
u_h = 1
u_exd = 0.7

#Probability
p <- c(0.15, 0.25)

#Correct the rates to cycle length such that it is for each cycle
for( i in 1:11){
  r_a <- -log(1-m[i])
  r_c <- r_a * cl
  p_c <- 1-exp(-r_c)
  m[i] <- p_c
  
  
  m_exd[i] <- 1-(1-p_c)^4
  
  rr <- (1-exp(-4*r_c))/(1-exp(-r_c))
  m_exd2[i] <- rr*p_c
}

for( i in 1:2){
  r_a <- -log(1-p[i])
  r_c <- r_a * cl
  p_c = 1-exp(-r_c)
  p[i] = p_c
}


#The discount rate per cycle (assuming an annual discount rate of 3%)
dr_c <- (1+0.03)^(1/12) - 1

dr_v <- 0:600
dr_v[1] <- 1




df <- data.frame(
  cycle = 0:600,
  healthy = 0:600,
  acuteEXD = 0:600,
  dead = 0:600
)

df_plot <- data.frame(
  cycle = rep(0:600, 3),
  state = 0:1802,
  value = 0:1802
)


i <- 1 #Variable for counting where we are in the m array

ly <- 1:600

#The initial state
s0 <- matrix(c(1,0,0), ncol = 3)

#change <- c(0, 59, 119, 179, 239, 299, 359, 419, 479, 539, 599)
change <- c(0, 60, 120, 180, 240, 300, 360, 420, 480, 540)

#The Markov simulation
for (t in 0:(12*50 - 1)) {
  
  dr_v[t+2] <- (1/(1+0.03))^((t+1)*cl)
  
  #Update the transition matrix every 5 year
  if(t %in% change){
    l <- 1-m[i]
    A <- matrix(c(l*(1-p[1]), (1-m_exd[i])*p[2], 0, l*p[1],  1-m_exd[i]-(1-m_exd[i])*p[2], 0, m[i], m_exd[i], 1), nrow = 3)
    i <- i + 1
    #print("This is A:")
    #print(A)
  }
  
  #Do one cycle
  s1 <- s0 %*% A
  #Update the data-frame
  df["healthy"][df["healthy"] == t+1] <- s1[1,1]
  df["acuteEXD"][df["acuteEXD"] == t+1] <- s1[1,2]
  df["dead"][df["dead"] == t+1] <- s1[1,3] 
  
  df_plot["state"][df_plot["value"] == t+1] <- "Healthy"
  df_plot["state"][df_plot["value"] == t+1 + 601] <- "Acute EXD"
  df_plot["state"][df_plot["value"] == t+1 + 1202] <- "Dead"
  df_plot["value"][df_plot["value"] == t+1] <- s1[1,1]
  df_plot["value"][df_plot["value"] == t+1 + 601] <- s1[1,2]
  df_plot["value"][df_plot["value"] == t+1 + 1202] <- s1[1,3]
  
  s0 <- s1
  
}

df["healthy"][df["healthy"] == 0] <- 1
df_plot["value"][df_plot["state"] == 0] <- 1
df_plot["value"][df_plot["state"] == 601] <- 0
df_plot["value"][df_plot["state"] == 1202] <- 0
df_plot["state"][df_plot["state"] == 0] <- "Healthy"
df_plot["state"][df_plot["state"] == 601] <- "Acute EXD"
df_plot["state"][df_plot["state"] == 1202] <- "Dead"


#plot(df$cycle, df$healthy, type = "l", col = 1)
#lines(df$cycle, df$acuteEXD, type = "l", col = 2)
#lines(df$cycle, df$dead, type = "l", col = 3)


library(ggplot2)


ggplot(df_plot, aes(x = cycle,y = value, group = state, color = state)) +
  geom_line()

print((sum(df[2:601,]$healthy) + sum(df[2:601,]$acuteEXD) + df[1,]$healthy/2 )/12)

print((sum(df$healthy)*1 + sum(df$acuteEXD)*0.7 - df[1,]$healthy/2)/12)

print((sum(df$healthy)*1500 + sum(df$acuteEXD)*3000 - (df[1,]$healthy/2)*1500)/12)


print((sum(df$healthy*dr_v) + sum(df$acuteEXD*dr_v) - df[1,]$healthy/2 )/12)

print((sum(df$healthy*dr_v) + sum(df$acuteEXD*dr_v)*0.7 - df[1,]$healthy/2)/12)

print((sum(df$healthy*dr_v)*1500 + sum(df$acuteEXD*dr_v)*3000 - (df[1,]$healthy/2)*1500)/12)

