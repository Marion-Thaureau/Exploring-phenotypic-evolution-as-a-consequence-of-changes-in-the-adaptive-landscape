###############################################
##                                           ##
##  SALMON RESULT TABLE FOR STANDARD ERRORS  ##
##                                           ## 
###############################################

#R version 4.2.1

###### LOAD BEST RESULT FOR EACH RUN ######
setwd("../../Results/Salmon_results/S_Rdatafiles")

load("salmon_bestrun_model_1a.RData")
load("salmon_bestrun_model_1b.RData")

load("salmon_bestrun_model_2a_it1.RData")
load("salmon_bestrun_model_2a_it50.RData")
load("salmon_bestrun_model_2b_it1.RData")
load("salmon_bestrun_model_2b_it50.RData")

load("salmon_bestrun_model_3a_it1.RData")
load("salmon_bestrun_model_3a_it50.RData")
load("salmon_bestrun_model_3b_it1.RData")
load("salmon_bestrun_model_3b_it50.RData")

load("salmon_bestrun_model_4a_it1.RData")
load("salmon_bestrun_model_4a_it50.RData")

load("salmon_bestrun_model_5a_it1.RData")
load("salmon_bestrun_model_5a_it50.RData")
load("salmon_bestrun_model_5b_it1_sp.RData")
load("salmon_bestrun_model_5b_it50_sp.RData")



###### SELECT THE BEST ITERATION FOR EACH MODEL ######
salmon_model_1a_best <- salmon_model_1a
salmon_model_1a_best$SE.A <- matrix(0, nrow = 2, ncol = 2)
salmon_model_1a_best$SE.optima <- rep(NA, 2)

salmon_model_1b_best <- salmon_model_1b
salmon_model_1b_best$SE.A <- matrix(0, nrow = 2, ncol = 2)
salmon_model_1b_best$SE.optima <- rep(NA, 2)

if (salmon_model_2a_best_it1$result2$AICc < salmon_model_2a_best_it50$result2$AICc) {
  salmon_model_2a_best <- salmon_model_2a_best_it1$result2
} else {
  salmon_model_2a_best <- salmon_model_2a_best_it50$result2
}

if (salmon_model_2b_best_it1$result2$AICc < salmon_model_2b_best_it50$result2$AICc) {
  salmon_model_2b_best <- salmon_model_2b_best_it1$result2
} else {
  salmon_model_2b_best <- salmon_model_2b_best_it50$result2
}

if (salmon_model_3a_best_it1$result2$AICc < salmon_model_3a_best_it50$result2$AICc) {
  salmon_model_3a_best <- salmon_model_3a_best_it1$result2
} else {
  salmon_model_3a_best <- salmon_model_3a_best_it50$result2
}

if (salmon_model_3b_best_it1$result2$AICc < salmon_model_3b_best_it50$result2$AICc) {
  salmon_model_3b_best <- salmon_model_3b_best_it1$result2
} else {
  salmon_model_3b_best <- salmon_model_3b_best_it50$result2
}

if (salmon_model_4a_best_it1$result2$AICc < salmon_model_4a_best_it50$result2$AICc) {
  salmon_model_4a_best <- salmon_model_4a_best_it1$result2
} else {
  salmon_model_4a_best <- salmon_model_4a_best_it50$result2
}

if (salmon_model_5a_best_it1$result2$AICc < salmon_model_5a_best_it50$result2$AICc) {
  salmon_model_5a_best <- salmon_model_5a_best_it1$result2
} else {
  salmon_model_5a_best <- salmon_model_5a_best_it50$result2
}

if (salmon_model_5b_best_it1$result2$AICc < salmon_model_5b_best_it50$result2$AICc) {
  salmon_model_5b_best <- salmon_model_5b_best_it1$result2
} else {
  salmon_model_5b_best <- salmon_model_5b_best_it50$result2
}



###### EXTRACT BEST PARAMETERS FOR EACH MODEL ######
models_list <- c("1a", "1b", "2a", "2b", "3a", "3b", "4a", "5a", "5b")

#Difference in AICc with the best model, DeltaAICc
AICc_list <- c()
for (model in models_list) {
  model_name <- paste0("salmon_model_", model, "_best")
  AICc_list[model] <- get(model_name)$AICc
}

AICc_best = min(AICc_list)

deltaAICc_list <- c()
for (model in models_list) {
  model_name <- paste0("salmon_model_", model, "_best")
  deltaAICc_list[model] <- AICc_list[model] - AICc_best
}

#Save the best model
model_best = names(which.min(AICc_list))
name_best_model = paste0("salmon_model_", model_best, "_best")
salmon_best_model = get(name_best_model)
save(salmon_best_model, file = 'salmon_bestmodel.RData')

#number of parameters k and loglikelihood estimates
k_list <- c()
logL_list <- c()

for (model in models_list) {
  model_name <- paste0("salmon_model_", model, "_best")
  k_list[model] <- get(model_name)$K
  logL_list[model] <- get(model_name)$logL
}

#Ancestral and optimum values for each time series
anc_trait_list <- c()
anc_envi_list <- c()
opt_trait_list <- c()
opt_envi_list <- c()

for (model in models_list) {
  model_name <- paste0("salmon_model_", model, "_best")
  anc_trait_list[model] <- get(model_name)$SE.anc[[1]]
  anc_envi_list[model] <- get(model_name)$SE.anc[[2]]
  opt_trait_list[model] <- get(model_name)$SE.optima[[1]]
  opt_envi_list[model] <- get(model_name)$SE.optima[[2]]
}

#Elements of the A matrix (deterministic values)
α11_list <- c()
α22_list <- c()
α12_list <- c()

for (model in models_list) {
  model_name <- paste0("salmon_model_", model, "_best")
  α11_list[model] <- get(model_name)$SE.A[1,1]
  α22_list[model] <- get(model_name)$SE.A[2,2]
  α12_list[model] <- get(model_name)$SE.A[1,2]
}

#Elements of the R matrix (stochastic values)
σ11_list <- c()
σ22_list <- c()
σ12_list <- c()

for (model in models_list) {
  model_name <- paste0("salmon_model_", model, "_best")
  σ11_list[model] <- get(model_name)$SE.R[1,1]
  σ22_list[model] <- get(model_name)$SE.R[2,2]
  σ12_list[model] <- get(model_name)$SE.R[1,2]
}


#Round values 
deltaAICc_list <- round(deltaAICc_list, 2)
k_list <- round(k_list, 2)
logL_list <- round(logL_list, 2)
anc_trait_list <- round(anc_trait_list, 2)
anc_envi_list <- round(anc_envi_list, 2)
opt_trait_list <- round(opt_trait_list, 2)
opt_envi_list <- round(opt_envi_list, 2)
α11_list <- round(α11_list, 2)
α22_list <- round(α22_list, 2)
α12_list <- round(α12_list, 2)
σ11_list <- round(σ11_list, 2)
σ22_list <- round(σ22_list, 2)
σ12_list <- round(σ12_list, 2)



###### CREATE TABLE OF RESULTS ######
salmon_result_table_SE <- data.frame(
  Model.k = mapply(function(model, k) paste0(model, " (", k, ")"), models_list, k_list),
  DAICc.logL = mapply(function(deltaAICc, logL) paste0(deltaAICc, "/", logL), deltaAICc_list, logL_list),
  anc1.anc2 = mapply(function(anc1, anc2) paste0("(", anc1, ", ", anc2, ")"), anc_trait_list, anc_envi_list),
  opt1.opt2 = mapply(function(opt1, opt2) paste0("(", opt1, ", ", opt2, ")"), opt_trait_list, opt_envi_list),
  a11.a22 = mapply(function(α11, α22) paste0("(", α11, ", ", α22, ")"), α11_list, α22_list),
  a12 = α12_list,
  s11.s22 = mapply(function(σ11, σ22) paste0("(", σ11, ", ", σ22, ")"), σ11_list, σ22_list),
  s12 = σ12_list
)

#Export the table of results
write.csv(salmon_result_table_SE, file = "../salmon_results_table_SE.csv", row.names = FALSE)

