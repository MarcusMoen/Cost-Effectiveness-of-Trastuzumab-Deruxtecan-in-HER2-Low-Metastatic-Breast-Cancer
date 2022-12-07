source("finding_transition_probabilities.R")

opt_val <- opt_var

tm_opt <- transition_matrices(opt_val)

hr <- seq(0.4, 0.9, by=0.01)
cost_pf_tdxd_l <- seq(8000, 20000, by=250)

n_cycles <- 120

df <- data.frame(
  Trial = 1:(length(hr)*length(cost_pf_tdxd_l)),
  HR = 1:(length(hr)*length(cost_pf_tdxd_l)),
  CostPF_tdxd_l = 1:(length(hr)*length(cost_pf_tdxd_l)),
  ICER = 1:(length(hr)*length(cost_pf_tdxd_l)),
  ICER_Range = 1:(length(hr)*length(cost_pf_tdxd_l))
)

df_eribulin <- data.frame(
  Trial = 1:(length(hr)*length(cost_pf_tdxd_l)),
  HR = 1:(length(hr)*length(cost_pf_tdxd_l)),
  CostPF_tdxd_l = 1:(length(hr)*length(cost_pf_tdxd_l)),
  ICER = 1:(length(hr)*length(cost_pf_tdxd_l)),
  ICER_Range = 1:(length(hr)*length(cost_pf_tdxd_l))
)

df_paclitaxel <- data.frame(
  Trial = 1:(length(hr)*length(cost_pf_tdxd_l)),
  HR = 1:(length(hr)*length(cost_pf_tdxd_l)),
  CostPF_tdxd_l = 1:(length(hr)*length(cost_pf_tdxd_l)),
  ICER = 1:(length(hr)*length(cost_pf_tdxd_l)),
  ICER_Range = 1:(length(hr)*length(cost_pf_tdxd_l))
)


for(i in 1:length(hr)){
  if(i%%50 == 0){
    print(i)
  }
  vals <- opt_var
  vals[3] <- hr[i]
  tm <- transition_matrices(vals)
  A_chemo <- tm[[1]]
  A_tdxd <- tm[[2]]
  df_list_chemo <- cea3state(A_chemo, n_cycles)
  df_chemo <- df_list_chemo[[1]]
  df_list_tdxd <- cea3state(A_tdxd, n_cycles)
  df_tdxd <- df_list_tdxd[[1]]
  
  for(j in 1:length(cost_pf_tdxd_l)){
    cost_pf_chemo <- 7776
    cost_pf_eribulin <- 9854.82
    cost_pf_paclitaxel <- 980.31
    cost_p_chemo <- 16058.83
    cost_pf_tdxd <- cost_pf_tdxd_l[j]
    cost_p_tdxd <- 16058.83
    dr_v <- df_list_chemo[[3]]
    
    cost_chemo <- (sum(df_chemo$ProgressionFree)*cost_pf_chemo + sum(df_chemo$Progressed)*cost_p_chemo - (df_chemo[1,]$ProgressionFree/2)*cost_pf_chemo)/12
    cost_eribulin <- (sum(df_chemo$ProgressionFree)*cost_pf_eribulin + sum(df_chemo$Progressed)*cost_p_chemo - (df_chemo[1,]$ProgressionFree/2)*cost_pf_eribulin)/12
    cost_paclitaxel <- (sum(df_chemo$ProgressionFree)*cost_pf_paclitaxel + sum(df_chemo$Progressed)*cost_p_chemo - (df_chemo[1,]$ProgressionFree/2)*cost_pf_paclitaxel)/12
    cost_tdxd <- (sum(df_tdxd$ProgressionFree)*cost_pf_tdxd + sum(df_chemo$Progressed)*cost_p_tdxd - (df_tdxd[1,]$ProgressionFree/2)*cost_pf_tdxd)/12
    
    cost_chemo_d <- (sum(df_chemo$ProgressionFree*dr_v)*cost_pf_chemo + sum(df_chemo$Progressed*dr_v)*cost_p_chemo - (df_chemo[1,]$ProgressionFree/2)*cost_pf_chemo)/12
    cost_eribulin_d <- (sum(df_chemo$ProgressionFree*dr_v)*cost_pf_eribulin + sum(df_chemo$Progressed*dr_v)*cost_p_chemo - (df_chemo[1,]$ProgressionFree/2)*cost_pf_eribulin)/12
    cost_paclitaxel_d <- (sum(df_chemo$ProgressionFree*dr_v)*cost_pf_paclitaxel + sum(df_chemo$Progressed*dr_v)*cost_p_chemo - (df_chemo[1,]$ProgressionFree/2)*cost_pf_paclitaxel)/12
    cost_tdxd_d <- (sum(df_tdxd$ProgressionFree*dr_v)*cost_pf_tdxd + sum(df_chemo$Progressed*dr_v)*cost_p_tdxd - (df_tdxd[1,]$ProgressionFree/2)*cost_pf_tdxd)/12
    
    qaly_pf_chemo <- 0.601
    qaly_p_chemo <- 0.54
    qaly_pf_tdxd <- 0.602
    qaly_p_tdxd <- 0.54
    
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
    
    df["HR"][df["HR"] == (i-1)*length(cost_pf_tdxd_l)+j] <- hr[i]
    df["CostPF_tdxd_l"][df["CostPF_tdxd_l"] == (i-1)*length(cost_pf_tdxd_l)+j] <- cost_pf_tdxd_l[j]
    df["ICER"][df["ICER"] == (i-1)*length(cost_pf_tdxd_l)+j] <- df_res$ICER[2]
    
    df_res_eribulin <- data.frame(
      Strategy = c("Chemo:", "T-Dxd:"),
      Costs = c(cost_chemo, cost_tdxd),
      DiscountedCost = c(cost_chemo_d, cost_tdxd_d),
      LifeYears = c(ly_chemo, ly_tdxd),
      QALYs = c(qaly_chemo, qaly_tdxd),
      DiscountedQALY = c(qaly_chemo_d, qaly_tdxd_d),
      Incremental_Costs = c(0, cost_tdxd-cost_chemo),
      Incremental_QALYs = c(0, qaly_tdxd-qaly_chemo),
      ICER = c(0, (cost_tdxd_d-cost_eribulin_d)/(qaly_tdxd_d-qaly_chemo_d))
    )
    
    df_eribulin["HR"][df_eribulin["HR"] == (i-1)*length(cost_pf_tdxd_l)+j] <- hr[i]
    df_eribulin["CostPF_tdxd_l"][df_eribulin["CostPF_tdxd_l"] == (i-1)*length(cost_pf_tdxd_l)+j] <- cost_pf_tdxd_l[j]
    df_eribulin["ICER"][df_eribulin["ICER"] == (i-1)*length(cost_pf_tdxd_l)+j] <- df_res_eribulin$ICER[2]
    
    df_res_paclitaxel <- data.frame(
      Strategy = c("Chemo:", "T-Dxd:"),
      Costs = c(cost_chemo, cost_tdxd),
      DiscountedCost = c(cost_chemo_d, cost_tdxd_d),
      LifeYears = c(ly_chemo, ly_tdxd),
      QALYs = c(qaly_chemo, qaly_tdxd),
      DiscountedQALY = c(qaly_chemo_d, qaly_tdxd_d),
      Incremental_Costs = c(0, cost_tdxd-cost_chemo),
      Incremental_QALYs = c(0, qaly_tdxd-qaly_chemo),
      ICER = c(0, (cost_tdxd_d-cost_paclitaxel_d)/(qaly_tdxd_d-qaly_chemo_d))
    )
    
    df_paclitaxel["HR"][df_paclitaxel["HR"] == (i-1)*length(cost_pf_tdxd_l)+j] <- hr[i]
    df_paclitaxel["CostPF_tdxd_l"][df_paclitaxel["CostPF_tdxd_l"] == (i-1)*length(cost_pf_tdxd_l)+j] <- cost_pf_tdxd_l[j]
    df_paclitaxel["ICER"][df_paclitaxel["ICER"] == (i-1)*length(cost_pf_tdxd_l)+j] <- df_res_paclitaxel$ICER[2]

  }
}


for(i in 1:nrow(df)){
  if(df["ICER"][df["ICER_Range"] == i] < 50000){
    df["ICER_Range"][df["ICER_Range"] == i] = "Under 50k"
  }else if(df["ICER"][df["ICER_Range"] == i] < 75000){
    df["ICER_Range"][df["ICER_Range"] == i] = "Between 50k and 75k"
  }else if(df["ICER"][df["ICER_Range"] == i] < 100000){
    df["ICER_Range"][df["ICER_Range"] == i] = "Between 75k and 100k"
  }else{
    df["ICER_Range"][df["ICER_Range"] == i] = "Over 100k"
  }
}

df["ICER_Range"][df["HR"] == 0.5 & df["CostPF_tdxd_l"] == 14500] = "Base Case"

for(i in 1:nrow(df_eribulin)){
  if(df_eribulin["ICER"][df_eribulin["ICER_Range"] == i] < 50000){
    df_eribulin["ICER_Range"][df_eribulin["ICER_Range"] == i] = "Under 50k"
  }else if(df_eribulin["ICER"][df_eribulin == i] < 75000){
    df_eribulin["ICER_Range"][df_eribulin["ICER_Range"] == i] = "Between 50k and 75k"
  }else if(df_eribulin["ICER"][df_eribulin["ICER_Range"] == i] < 100000){
    df_eribulin["ICER_Range"][df_eribulin["ICER_Range"] == i] = "Between 75k and 100k"
  }else{
    df_eribulin["ICER_Range"][df_eribulin["ICER_Range"] == i] = "Over 100k"
  }
}

df_eribulin["ICER_Range"][df_eribulin["HR"] == 0.5 & df_eribulin["CostPF_tdxd_l"] == 14500] = "Base Case"

for(i in 1:nrow(df_paclitaxel)){
  if(df_paclitaxel["ICER"][df_paclitaxel["ICER_Range"] == i] < 50000){
    df_paclitaxel["ICER_Range"][df_paclitaxel["ICER_Range"] == i] = "Under 50k"
  }else if(df_paclitaxel["ICER"][df_paclitaxel["ICER_Range"] == i] < 75000){
    df_paclitaxel["ICER_Range"][df_paclitaxel["ICER_Range"] == i] = "Between 50k and 75k"
  }else if(df_paclitaxel["ICER"][df_paclitaxel["ICER_Range"] == i] < 100000){
    df_paclitaxel["ICER_Range"][df_paclitaxel["ICER_Range"] == i] = "Between 75k and 100k"
  }else{
    df_paclitaxel["ICER_Range"][df_paclitaxel["ICER_Range"] == i] = "Over 100k"
  }
}

df_paclitaxel["ICER_Range"][df_paclitaxel["HR"] == 0.5 & df_paclitaxel["CostPF_tdxd_l"] == 14500] = "Base Case"


library(ggplot2)
p1 <- ggplot(df, aes(x = HR, y = CostPF_tdxd_l, group = ICER_Range, color = ICER_Range)) +
  geom_point(shape = 15, size = 4.5) + 
  #geom_point(aes(x=0.5, y=14000), colour="black")+
  ggtitle("Two-way sensitivity analysis (T-DxD vs Physican's Choice)") +
  xlab("HR") + ylab("Cost of T-DxD") + 
  theme(
    plot.title = element_text(color="dodgerblue4", size=16, face="bold.italic"),
    axis.title.x = element_text(color="darkolivegreen4", size=14, face="bold"),
    axis.title.y = element_text(color="#993333", size=14, face="bold")
  )

p1 + scale_color_brewer(palette = "Spectral")


p2 <- ggplot(df_eribulin, aes(x = HR, y = CostPF_tdxd_l, group = ICER_Range, color = ICER_Range)) +
  geom_point(shape = 15, size = 4.5) + 
  #geom_point(aes(x=0.5, y=14000), colour="black")+
  ggtitle("Two-way sensitivity analysis (T-DxD vs Expensive Chemo)") +
  xlab("HR") + ylab("Cost of T-DxD") + 
  theme(
    plot.title = element_text(color="dodgerblue4", size=16, face="bold.italic"),
    axis.title.x = element_text(color="darkolivegreen4", size=14, face="bold"),
    axis.title.y = element_text(color="#993333", size=14, face="bold")
  )

p2 + scale_color_brewer(palette = "RdBu")


p3 <- ggplot(df_paclitaxel, aes(x = HR, y = CostPF_tdxd_l, group = ICER_Range, color = ICER_Range)) +
  geom_point(shape = 15, size = 4.5) + 
  #geom_point(aes(x=0.5, y=14000), colour="black")+
  ggtitle("Two-way sensitivity analysis (T-DxD vs Cheap Chemo)") +
  xlab("HR") + ylab("Cost of T-DxD") + 
  theme(
    plot.title = element_text(color="dodgerblue4", size=16, face="bold.italic"),
    axis.title.x = element_text(color="darkolivegreen4", size=14, face="bold"),
    axis.title.y = element_text(color="#993333", size=14, face="bold")
  )

p3 + scale_color_brewer(palette = 4)







