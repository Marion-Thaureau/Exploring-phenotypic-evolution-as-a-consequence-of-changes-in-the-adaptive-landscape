######################################
##                                  ##
##  FORAMS - SENSITIVITY ANALYSIS   ##
##         iterations 50 - 75       ##  
##                                  ##
######################################

#R version 4.2.1
.libPaths("/cluster/home/marionth/R")
library(evoTS) #evoTS version 1.0.3

setwd("/cluster/home/marionth/Detrending")


###### GET THE OUTPUT FOR THE FORAMINIFERA BEST MODEL ######
# Load evoTS object and best model output for the foraminifera dataset
load("foraminifera_evoTS_obj.RData")
load("foraminifera_bestmodel.RData")

# Get the parameter estimates for the best model
A_matrix<-foraminifera_best_model$A[1:2, 1:2]
R_matrix<-foraminifera_best_model$R[1:2, 1:2]
R_matrix[2,2] <- 1.20 # Making 'Sigma' positive definite
anc = foraminifera_best_model$anc[1:2]
optima = foraminifera_best_model$opt[1:2]
ts_length = length(ln.bodylength_isotopes_multiobj$xx[,1])

nrep <- 25
winners <- character(nrep)   # store best model name each run
set.seed(250)
i <- 0
attempt <- 0

###### SIMULATE THE TIME SERIES AND FIT MODELS TO THE DATASET ######
while (i < nrep) {   
  attempt <- attempt + 1
  
  # ---- Simulate the time series ----
  data_sim<-sim.multi.OU(ts_length, anc = anc, optima = optima, A = A_matrix, R = R_matrix, vp = 0.000001)
  data_sim$tt[,1] <- data_sim$tt[,1]/(max(data_sim$tt[,1])) # rescale the time vector
  data_sim$tt[,2] <- data_sim$tt[,1]
    
  # ---- Fit models ----
  m1<-fit.multivariate.URW(data_sim, R="diag")
  print(m1$converge)
  if (m1$converge != "Model converged successfully") { message("Attempt ", attempt, ": m1 did not converge, re-simulating."); next }   ###
  m2 <- fit.multivariate.OU(data_sim, A.matrix="diag", R.matrix="diag")
  print(m2$converge)
  if (m2$converge != "Model converged successfully") { message("Attempt ", attempt, ": m2 did not converge, re-simulating."); next }   ###
  m3 <- fit.multivariate.OU(data_sim, A.matrix="upper.tri", R.matrix="diag")
  print(m3$converge)
  if (m3$converge != "Model converged successfully") { message("Attempt ", attempt, ": m3 did not converge, re-simulating."); next }   ###
  
  i <- i + 1   ###
  
  # ---- Extract AICcs ----
  aics <- c(
    m1 = m1$AICc,
    m2 = m2$AICc,
    m3 = m3$AICc
  )
  
  print(i)
  print(aics)
  
  # ---- Store best model ----
  winners[i] <- names(which.min(aics))
}

# ---- Count how often each model wins ----
table(winners)

save(winners, file = "forams_sensitivity_analysis_it50.RData")
