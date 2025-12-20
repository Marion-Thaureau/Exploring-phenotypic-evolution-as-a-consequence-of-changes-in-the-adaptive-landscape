###########################################################
##                                                       ##
##  FORAMINIFERA SCRIPT - Model 6a (iteration 50 - 100)  ##
##                                                       ## 
###########################################################


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


setwd("/cluster/home/marionth/Hodell")



###### LOAD MORPHOLOGICAL AND ENVIRONMENTAL DATA ######
bodylength_raw <- read.table("Hodell_and_Vayavananda_1993_Fohsella_length.txt",header=T)
oxygen_isotopes_raw <- read.csv("table_d18O.csv")
carbon_isotopes_raw <- read.csv("table_d13C.csv")

#removing the NA rows
dataset_raw = cbind(bodylength_raw, oxygen_isotopes_raw, carbon_isotopes_raw)
dataset = dataset_raw %>% filter(!is.na(d13C_mean)) 

#creating the body length time series (already log-transformed)
bodylength_ts <- data.frame(
  mean = dataset$trait_mean,
  variance = dataset$trait_var,
  samplesize = dataset$N,
  time = dataset$age_MY
)

#creating the oxygen isotope time series
oxygen_isotopes_ts <- data.frame(
  mean = dataset$d18O_mean,
  variance = abs(dataset$d18O_mean)*0.001,
  samplesize = 1,
  time = dataset$age_MY
)

#creating the carbon isotope time series
carbon_isotopes_ts <- data.frame(
  mean = dataset$d13C_mean,
  variance = abs(dataset$d13C_mean)*0.001,
  samplesize = 1,
  time = dataset$age_MY
)



###### CREATING PALEOTS AND EVOTS OBJECTS ######
#Creating paleoTS object for the body length time series
ln.bodylength_obj <- paleoTS::as.paleoTS(mm = bodylength_ts$mean, vv = bodylength_ts$variance, nn = bodylength_ts$samplesize, tt = bodylength_ts$time)
ln.bodylength_obj$tt <- ln.bodylength_obj$tt/(max(ln.bodylength_obj$tt)) #transfer of the time vector into a length unit of 1 

#Creating paleoTS object for the oxygen isotope time series
oxygen_isotopes_obj <- paleoTS::as.paleoTS(mm = oxygen_isotopes_ts$mean, vv = oxygen_isotopes_ts$variance, nn = oxygen_isotopes_ts$samplesize, tt = oxygen_isotopes_ts$time, oldest="first")
oxygen_isotopes_obj$tt <- oxygen_isotopes_obj$tt/(max(oxygen_isotopes_obj$tt)) #transfer of the time vector into a length unit of 1

#Creating paleoTS object for the carbon isotope time series
carbon_isotopes_obj <- paleoTS::as.paleoTS(mm = carbon_isotopes_ts$mean, vv = carbon_isotopes_ts$variance, nn = carbon_isotopes_ts$samplesize, tt = carbon_isotopes_ts$time, oldest="first")
carbon_isotopes_obj$tt <- carbon_isotopes_obj$tt/(max(carbon_isotopes_obj$tt)) #transfer of the time vector into a length unit of 1

#Creating and plotting the multivariate evoTS object
ln.bodylength_isotopes_multiobj <- make.multivar.evoTS(ln.bodylength_obj, oxygen_isotopes_obj, carbon_isotopes_obj)



###### RUNNING Model_6a ######
#Parallelisation of the job 
registerDoParallel(52)
start_time <- Sys.time() 

#Model_6a: both isotopes influence the primary optimum of the body length in A matrix; diagonal R matrix 
A = matrix(c(1,1,1,0,1,0,0,0,1), nrow = 3, byrow=TRUE)
R = matrix(c(1,0,0,0,1,0,0,0,1), nrow = 3, byrow=TRUE)

#Loop for 65 iterations 
foraminifera_model_6a_it50 <- foreach(i=1:104) %dopar% { 
  result <- multiResultClass() 
  
#Set and store a different seed for each iteration
  seed <- (114+i) 
  set.seed(seed)
  result$result1 <- seed
  
#Run iterations using the parametrized matrices A and R
  result$result2 <- tryCatch({fit.multivariate.OU.user.defined(ln.bodylength_isotopes_multiobj, A.user = A, R.user = R, iterations = 1, hess = TRUE,
                                                               user.init.diag.A = c(18.18, 8.65, 30.57),
                                                               user.init.upper.diag.A = NULL,
                                                               user.init.diag.R = c(0.39, 1.27, 1.62),
                                                               user.init.off.diag.R = NULL)

    
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

save.image(file = 'foraminifera_allruns_model_6a_it50_sp2.RData')



###### EXTRACTING AND SAVING THE ITERATION WITH BEST PARAMETER ESTIMATES ######
#Filtering 50 iterations which converged
foraminifera_model_6a_it50_filtered <- Filter(function(iteration) !is.na(iteration$result2[[1]]), foraminifera_model_6a_it50)
foraminifera_model_6a_it50_filtered <- Filter(function(iteration) !is.null(iteration), foraminifera_model_6a_it50_filtered)
foraminifera_model_6a_it50_results <- foraminifera_model_6a_it50_filtered[1:50]

#Extracting the AICc values of each iterations
AICc <- matrix(ncol=2, nrow=50)
colnames(AICc) <- c("iteration", "AIC")

for (j in 1:50) {
  AICc[j,1] <- j 
  AICc[j,2] <- foraminifera_model_6a_it50_results[[j]]$result2$AICc 
}

#Finding the best AICc
AICc_min <- which.min(AICc[,2])

#Finding the best iteration  
foraminifera_model_6a_best_it50 <- foraminifera_model_6a_it50_results[[AICc_min]]

#Saving the best iteration
print(foraminifera_model_6a_best_it50)
save(foraminifera_model_6a_best_it50, file = 'foraminifera_bestrun_model_6a_it50_sp2.RData')