##################################################
##                                              ##
##  COCCOLITH RESULT TABLE FOR STANDARD ERRORS  ##
##                                              ## 
##################################################

#R version 4.2.1

###### LOAD BEST RESULT FOR EACH RUN ######
setwd("../../Results/Coccolith_results/C_Rdatafiles")

load("coccolith_bestrun_model_1a.RData")
load("coccolith_bestrun_model_1b.RData")

load("coccolith_bestrun_model_2a_it1.RData")
load("coccolith_bestrun_model_2a_it50.RData")
load("coccolith_bestrun_model_2b_it1.RData")
load("coccolith_bestrun_model_2b_it10.RData")
load("coccolith_bestrun_model_2b_it20.RData")
load("coccolith_bestrun_model_2b_it30.RData")
load("coccolith_bestrun_model_2b_it40.RData")
load("coccolith_bestrun_model_2b_it50.RData")
load("coccolith_bestrun_model_2b_it60.RData")
load("coccolith_bestrun_model_2b_it70.RData")
load("coccolith_bestrun_model_2b_it80.RData")
load("coccolith_bestrun_model_2b_it90.RData")

load("coccolith_bestrun_model_3a_it1.RData")
load("coccolith_bestrun_model_3a_it50.RData")
load("coccolith_bestrun_model_3b_it1.RData")
load("coccolith_bestrun_model_3b_it50.RData")

load("coccolith_bestrun_model_4a_it1.RData")
load("coccolith_bestrun_model_4a_it50.RData")
load("coccolith_bestrun_model_4b_it1.RData")
load("coccolith_bestrun_model_4b_it50.RData")

load("coccolith_bestrun_model_5a_it1.RData")
load("coccolith_bestrun_model_5a_it50.RData")
load("coccolith_bestrun_model_5b_it1.RData")
load("coccolith_bestrun_model_5b_it25.RData")
load("coccolith_bestrun_model_5b_it50.RData")
load("coccolith_bestrun_model_5b_it75.RData")

load("coccolith_bestrun_model_6a_it1_sp2.RData")
load("coccolith_bestrun_model_6a_it25_sp2.RData")
load("coccolith_bestrun_model_6a_it50_sp2.RData")
load("coccolith_bestrun_model_6a_it75_sp2.RData")
load("coccolith_bestrun_model_6b_it1.RData")
load("coccolith_bestrun_model_6b_it25.RData")
load("coccolith_bestrun_model_6b_it50.RData")
load("coccolith_bestrun_model_6b_it75.RData")



###### SELECT THE BEST ITERATION FOR EACH MODEL ######
coccolith_model_1a_best <- coccolith_model_1a
coccolith_model_1a_best$SE.A <- matrix(0, nrow = 3, ncol = 3)
coccolith_model_1a_best$SE.optima <- rep(NA, 3)

coccolith_model_1b_best <- coccolith_model_1b
coccolith_model_1b_best$SE.A <- matrix(0, nrow = 3, ncol = 3)
coccolith_model_1b_best$SE.optima <- rep(NA, 3)

if (coccolith_model_2a_best_it1$result2$AICc < coccolith_model_2a_best_it50$result2$AICc) {
  coccolith_model_2a_best <- coccolith_model_2a_best_it1$result2
} else {
  coccolith_model_2a_best <- coccolith_model_2a_best_it50$result2
}

#iterate through the 10 subsets of data for the model 2b
model_2b_itsubset <- c("1", "10", "20", "30", "40", "50", "60", "70", "80", "90")
model_2b_AICc_list <- c()
for (itsubset in model_2b_itsubset) {
  it_name <- paste0("coccolith_model_2b_best_it", itsubset)
  model_2b_AICc_list[itsubset] <- get(it_name)$result2$AICc
}
model_2b_AICc_min <- names(which.min(model_2b_AICc_list))
name_coccolith_model_2b_best = paste0("coccolith_model_2b_best_it", model_2b_AICc_min)
coccolith_model_2b_best = get(name_coccolith_model_2b_best)$result2

if (coccolith_model_3a_best_it1$result2$AICc < coccolith_model_3a_best_it50$result2$AICc) {
  coccolith_model_3a_best <- coccolith_model_3a_best_it1$result2
} else {
  coccolith_model_3a_best <- coccolith_model_3a_best_it50$result2
}

if (coccolith_model_3b_best_it1$result2$AICc < coccolith_model_3b_best_it50$result2$AICc) {
  coccolith_model_3b_best <- coccolith_model_3b_best_it1$result2
} else {
  coccolith_model_3b_best <- coccolith_model_3b_best_it50$result2
}

if (coccolith_model_4a_best_it1$result2$AICc < coccolith_model_4a_best_it50$result2$AICc) {
  coccolith_model_4a_best <- coccolith_model_4a_best_it1$result2
} else {
  coccolith_model_4a_best <- coccolith_model_4a_best_it50$result2
}

if (coccolith_model_4b_best_it1$result2$AICc < coccolith_model_4b_best_it50$result2$AICc) {
  coccolith_model_4b_best <- coccolith_model_4b_best_it1$result2
} else {
  coccolith_model_4b_best <- coccolith_model_4b_best_it50$result2
}

if (coccolith_model_5a_best_it1$result2$AICc < coccolith_model_5a_best_it50$result2$AICc) {
  coccolith_model_5a_best <- coccolith_model_5a_best_it1$result2
} else {
  coccolith_model_5a_best <- coccolith_model_5a_best_it50$result2
}

#iterate through the 4 subsets of data for the model 5b
model_5b_itsubset <- c("1", "25", "50", "75")
model_5b_AICc_list <- c()
for (itsubset in model_5b_itsubset) {
  it_name <- paste0("coccolith_model_5b_best_it", itsubset)
  model_5b_AICc_list[itsubset] <- get(it_name)$result2$AICc
}
model_5b_AICc_min <- names(which.min(model_5b_AICc_list))
name_coccolith_model_5b_best = paste0("coccolith_model_5b_best_it", model_5b_AICc_min)
coccolith_model_5b_best = get(name_coccolith_model_5b_best)$result2

#iterate through the 4 subsets of data for the model 6a
model_6a_itsubset <- c("1", "25", "50", "75")
model_6a_AICc_list <- c()
for (itsubset in model_6a_itsubset) {
  it_name <- paste0("coccolith_model_6a_best_it", itsubset)
  model_6a_AICc_list[itsubset] <- get(it_name)$result2$AICc
}
model_6a_AICc_min <- names(which.min(model_6a_AICc_list))
name_coccolith_model_6a_best = paste0("coccolith_model_6a_best_it", model_6a_AICc_min)
coccolith_model_6a_best = get(name_coccolith_model_6a_best)$result2

#iterate through the 4 subsets of data for the model 6b
model_6b_itsubset <- c("1", "25", "50", "75")
model_6b_AICc_list <- c()
for (itsubset in model_6b_itsubset) {
  it_name <- paste0("coccolith_model_6b_best_it", itsubset)
  model_6b_AICc_list[itsubset] <- get(it_name)$result2$AICc
}
model_6b_AICc_min <- names(which.min(model_6b_AICc_list))
name_coccolith_model_6b_best = paste0("coccolith_model_6b_best_it", model_6b_AICc_min)
coccolith_model_6b_best = get(name_coccolith_model_6b_best)$result2



###### EXTRACT BEST PARAMETERS FOR EACH MODEL ######
models_list <- c("1a", "1b", "2a", "2b", "3a", "3b", "4a", "4b", "5a", "5b", "6a", "6b")

#Difference in AICc with the best model, DeltaAICc
AICc_list <- c()
for (model in models_list) {
  model_name <- paste0("coccolith_model_", model, "_best")
  AICc_list[model] <- get(model_name)$AICc
}

AICc_best = min(AICc_list)

deltaAICc_list <- c()
for (model in models_list) {
  model_name <- paste0("coccolith_model_", model, "_best")
  deltaAICc_list[model] <- AICc_list[model] - AICc_best
}

#Save the best model
model_best = names(which.min(AICc_list))
name_best_model = paste0("coccolith_model_", model_best, "_best")
coccolith_best_model = get(name_best_model)
save(coccolith_best_model, file = 'coccolith_bestmodel.RData')

#number of parameters k and loglikelihood estimates
k_list <- c()
logL_list <- c()

for (model in models_list) {
  model_name <- paste0("coccolith_model_", model, "_best")
  k_list[model] <- get(model_name)$K
  logL_list[model] <- get(model_name)$logL
}

#Ancestral and optimum values for each time series
anc_trait_list <- c()
anc_envi1_list <- c()
anc_envi2_list <- c()
opt_trait_list <- c()
opt_envi1_list <- c()
opt_envi2_list <- c()

for (model in models_list) {
  model_name <- paste0("coccolith_model_", model, "_best")
  anc_trait_list[model] <- get(model_name)$SE.anc[[1]]
  anc_envi1_list[model] <- get(model_name)$SE.anc[[2]]
  anc_envi2_list[model] <- get(model_name)$SE.anc[[3]]
  opt_trait_list[model] <- get(model_name)$SE.optima[[1]]
  opt_envi1_list[model] <- get(model_name)$SE.optima[[2]]
  opt_envi2_list[model] <- get(model_name)$SE.optima[[3]]
}

#Elements of the A matrix (deterministic values)
α11_list <- c()
α22_list <- c()
α33_list <- c()
α12_list <- c()
α13_list <- c()

for (model in models_list) {
  model_name <- paste0("coccolith_model_", model, "_best")
  α11_list[model] <- get(model_name)$SE.A[1,1]
  α22_list[model] <- get(model_name)$SE.A[2,2]
  α33_list[model] <- get(model_name)$SE.A[3,3]
  α12_list[model] <- get(model_name)$SE.A[1,2]
  α13_list[model] <- get(model_name)$SE.A[1,3]
}

#Elements of the R matrix (stochastic values)
σ11_list <- c()
σ22_list <- c()
σ33_list <- c()
σ12_list <- c()
σ13_list <- c()
σ23_list <- c()

for (model in models_list) {
  model_name <- paste0("coccolith_model_", model, "_best")
  σ11_list[model] <- get(model_name)$SE.R[1,1]
  σ22_list[model] <- get(model_name)$SE.R[2,2]
  σ33_list[model] <- get(model_name)$SE.R[3,3]
  σ12_list[model] <- get(model_name)$SE.R[1,2]
  σ13_list[model] <- get(model_name)$SE.R[1,3]
  σ23_list[model] <- get(model_name)$SE.R[2,3]
}

#Round values 
deltaAICc_list <- round(deltaAICc_list, 2)
k_list <- round(k_list, 2)
logL_list <- round(logL_list, 2)
anc_trait_list <- round(anc_trait_list, 2)
anc_envi1_list <- round(anc_envi1_list, 2)
anc_envi2_list <- round(anc_envi2_list, 2)
opt_trait_list <- round(opt_trait_list, 2)  
opt_envi1_list <- round(opt_envi1_list, 2)
opt_envi2_list <- round(opt_envi2_list, 2)
α11_list <- round(α11_list, 2)
α22_list <- round(α22_list, 2)
α33_list <- round(α33_list, 2)
α12_list <- round(α12_list, 2)
α13_list <- round(α13_list, 2)
σ11_list <- round(σ11_list, 2)
σ22_list <- round(σ22_list, 2)
σ33_list <- round(σ33_list, 2)
σ12_list <- round(σ12_list, 2)
σ13_list <- round(σ13_list, 2)
σ23_list <- round(σ23_list, 2)



###### CREATE TABLE OF RESULTS ######
coccolith_result_table_SE <- data.frame(
  Model.k = mapply(function(model, k) paste0(model, " (", k, ")"), models_list, k_list),
  DAICc.logL = mapply(function(deltaAICc, logL) paste0(deltaAICc, "/", logL), deltaAICc_list, logL_list),
  anc1.anc2.anc3 = mapply(function(anc1, anc2, anc3) paste0("(", anc1, ", ", anc2, ", ", anc3, ")"), anc_trait_list, anc_envi1_list, anc_envi2_list),
  opt1.opt2.opt3 = mapply(function(opt1, opt2, opt3) paste0("(", opt1, ", ", opt2, ", ", opt3, ")"), opt_trait_list, opt_envi1_list, opt_envi2_list),
  a11.a22.a33 = mapply(function(α11, α22, α33) paste0("(", α11, ", ", α22, ", ", α33, ")"), α11_list, α22_list, α33_list),
  a12.a13 = mapply(function(α12, α13) paste0("(", α12, ", ", α13, ")"), α12_list, α13_list),
  s11.s22.s33 = mapply(function(σ11, σ22, σ33) paste0("(", σ11, ", ", σ22, ", ", σ33, ")"), σ11_list, σ22_list, σ33_list),
  s12.s13.s23 = mapply(function(σ12, σ13, σ23) paste0("(", σ12, ", ", σ13, ", ", σ23, ")"), σ12_list, σ13_list, σ23_list)
)

#Export the table of results
write.csv(coccolith_result_table_SE, file = "../coccolith_results_table_SE.csv", row.names = FALSE)

