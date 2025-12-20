#######################################################
##                                                   ##
##  COCCOLITH SCRIPT - Model 2b (iteration 40 - 50)  ##
##                                                   ## 
#######################################################


#R version 4.2.1
.libPaths("/cluster/home/marionth/R")
library(doParallel) #doParallel version 1.0.17
library(dplyr) #dplyr version 1.1.0
library(evoTS) #evoTS version 1.0.3
library(foreach) #foreach version 1.5.2
library(MCSPCD) #MCSPCD version 0.0.0.9000
library(paleoTS) #paleoTS version 0.6.2
library(parallel) #parallel version 4.2.1
library(stats) #stats version 4.2.1
library(usethis) #usethis version 2.1.6


setwd("/cluster/home/marionth/Bornemann")



###### LOAD MORPHOLOGICAL DATA ######
bodylength_raw <- read.delim("/cluster/home/marionth/Bornemann/Biscutum_TotalcoccolithLength.txt")

bodylength_ts <- data.frame(
  mean = bodylength_raw[[2]],
  variance = bodylength_raw[[3]],
  samplesize = bodylength_raw[[1]],
  time = bodylength_raw[[4]]
)



###### LOAD ENVIRONMENTAL DATA ######
isotopes <- read.delim("/cluster/home/marionth/Bornemann/isotopes.txt")

oxygen_isotopes_ts <- data.frame(
  mean = isotopes[[4]],
  variance = abs(isotopes[[4]])*0.001,
  samplesize = 1,
  time = bodylength_raw[[4]]
)

carbon_isotopes_ts <- data.frame(
  mean = isotopes[[3]],
  variance = abs(isotopes[[3]])*0.001,
  samplesize = 1,
  time = bodylength_raw[[4]]
)



###### CREATING PALEOTS AND EVOTS OBJECTS ######
#Creating paleoTS object for the body length time series
bodylength_obj <- paleoTS::as.paleoTS(mm = bodylength_ts$mean, vv = bodylength_ts$variance, nn = bodylength_ts$samplesize, tt = bodylength_ts$time)
ln.bodylength_obj <- paleoTS::ln.paleoTS(bodylength_obj) #log-transformation of the data
ln.bodylength_obj$tt <- ln.bodylength_obj$tt/(max(ln.bodylength_obj$tt)) #transfer of the time vector into a length unit of 1 

#Creating paleoTS object for the oxygen isotope time series
oxygen_isotopes_obj <- paleoTS::as.paleoTS(mm = oxygen_isotopes_ts$mean, vv = oxygen_isotopes_ts$variance, nn = oxygen_isotopes_ts$samplesize, tt = oxygen_isotopes_ts$time, oldest="first")
oxygen_isotopes_obj$tt <- oxygen_isotopes_obj$tt/(max(oxygen_isotopes_obj$tt)) #transfer of the time vector into a length unit of 1

#Creating paleoTS object for the carbon isotope time series
carbon_isotopes_obj <- paleoTS::as.paleoTS(mm = carbon_isotopes_ts$mean, vv = carbon_isotopes_ts$variance, nn = carbon_isotopes_ts$samplesize, tt = carbon_isotopes_ts$time, oldest="first")
carbon_isotopes_obj$tt <- carbon_isotopes_obj$tt/(max(carbon_isotopes_obj$tt)) #transfer of the time vector into a length unit of 1

#Creating and plotting the multivariate evoTS object
ln.bodylength_isotopes_multiobj <- make.multivar.evoTS(ln.bodylength_obj, oxygen_isotopes_obj, carbon_isotopes_obj)



###### RUNNING Model_2b ######
#Parallelisation of the job 
registerDoParallel(52)
start_time <- Sys.time() 

#Model_2b: only body length is an OU in A matrix; full R matrix 
A = matrix(c(1,0,0,0,0,0,0,0,0), nrow = 3, byrow=TRUE)
R = matrix(c(1,1,1,1,1,1,1,1,1), nrow = 3, byrow=TRUE)

#Loop for 100 iterations 
coccolith_model_2b_it40 <- foreach(i=1:100) %dopar% { 
  result <- multiResultClass() 
  
#Set and store a different seed for each iteration
  seed <- (410+i) 
  set.seed(seed)
  result$result1 <- seed
  
#Run iterations using the parametrized matrices A and R
  result$result2 <- tryCatch({fit.multivariate.OU.user.defined(ln.bodylength_isotopes_multiobj, A.user = A, R.user = R, iterations = 1, hess = TRUE)
    
#Handle non-converging iterations
  }, error = function(e) {
    message(paste("Error in iteration", i, ":", e$message))
    return(NA) 
  })
  
#Record the running time for each iteration
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  print(elapsed_time)
  result$result3 <- elapsed_time
  
  return(result)
}

print ("DONE")
end_time <- Sys.time(); end_time - start_time

save.image(file = 'coccolith_allruns_model_2b_it40.RData')



###### EXTRACTING AND SAVING THE ITERATION WITH BEST PARAMETER ESTIMATES ######
#Filtering 10 iterations which converged
coccolith_model_2b_it40_filtered <- Filter(function(iteration) !is.na(iteration$result2[[1]]), coccolith_model_2b_it40)
coccolith_model_2b_it40_filtered <- Filter(function(iteration) !is.null(iteration), coccolith_model_2b_it40_filtered)
coccolith_model_2b_it40_results <- coccolith_model_2b_it40_filtered[1:10]

#Extracting the AICc values of each iterations
AICc <- matrix(ncol=2, nrow=10)
colnames(AICc) <- c("iteration", "AIC")

for (j in 1:10) {
  AICc[j,1] <- j 
  AICc[j,2] <- coccolith_model_2b_it40_results[[j]]$result2$AICc 
}

#Finding the best AICc
AICc_min <- which.min(AICc[,2])

#Finding the best iteration  
coccolith_model_2b_best_it40 <- coccolith_model_2b_it40_results[[AICc_min]]

#calculating the correlation coefficients 
Rmatrix_best_coccolith_model_2b_it40 <- coccolith_model_2b_best_it40$result2$R
cormatrix_best_coccolith_model_2b_it40 <- cov2cor(Rmatrix_best_coccolith_model_2b_it40)

#Saving the best iteration
print(coccolith_model_2b_best_it40)
print(cormatrix_best_coccolith_model_2b_it40)
save(coccolith_model_2b_best_it40, cormatrix_best_coccolith_model_2b_it40, file = 'coccolith_bestrun_model_2b_it40.RData')