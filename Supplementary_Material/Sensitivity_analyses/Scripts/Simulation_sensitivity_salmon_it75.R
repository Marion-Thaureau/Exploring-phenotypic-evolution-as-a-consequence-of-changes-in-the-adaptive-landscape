######################################
##                                  ##
##  SALMON - SENSITIVITY ANALYSIS   ##
##        iterations 75 - 100       ##  
##                                  ##
######################################

#R version 4.2.1
.libPaths("/cluster/home/marionth/R")
library(evoTS) #evoTS version 1.0.3

setwd("/cluster/home/marionth/Detrending")


###### SIMULATE TIME SERIES WITH DEPENDENT ADAPTATION WITH OUTPUT OF THE SALMON BEST MODEL ######
# Load evoTS object and best model output for the salmon dataset
load("salmon_evoTS_obj.RData")
load("salmon_bestmodel.RData")

# Get the parameter estimates for matrices A and R
A_matrix<-salmon_best_model$A
R_matrix<-salmon_best_model$R

set.seed(123)

# Generate an evoTS object by simulating a multivariate dataset with parameter estimates of the salmon best model
ts_length = length(ln.bodymass_ln.discharge_multiobj$xx[,1])
data_sim<-sim.multi.OU(ts_length, anc = salmon_best_model$anc, optima = salmon_best_model$opt, A = A_matrix, R = R_matrix, vp = 0.000001)

data_sim$tt[,1] <- data_sim$tt[,1]/(max(data_sim$tt[,1])) # rescale the time vector
data_sim$tt[,2] <- data_sim$tt[,1]

nrep <- 25
winners <- character(nrep)   # store best model name each run
set.seed(375)
i <- 0
attempt <- 0

###### FIT MODELS TO THE DATASET WITHOUT DETRENDING ######
while (i < nrep) {   
  attempt <- attempt + 1
    
    # ---- Fit models ----
  m1<-fit.multivariate.URW(data_sim, R="diag")
  print(m1$converge)
  if (m1$converge != "Model converged successfully") { message("Attempt ", attempt, ": m1 did not converge, re-simulating."); next }   ###
  m2 <- fit.multivariate.OU(data_sim, A.matrix="diag", R.matrix="diag")
  print(m2$converge)
  if (m2$converge != "Model converged successfully") { message("Attempt ", attempt, ": m2 did not converge, re-simulating."); next }   ###
  m3 <- fit.multivariate.OU(data_sim, A.matrix="upper.tri", R.matrix="diag", trace = TRUE)
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

save(winners, file = "salmon_sensitivity_analysis_it75.RData")
