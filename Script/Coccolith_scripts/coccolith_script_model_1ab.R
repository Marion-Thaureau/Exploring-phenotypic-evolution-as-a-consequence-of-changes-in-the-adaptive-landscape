################################################
##                                            ##
##  COCCOLITH SCRIPT - Model 1a and Model 1b  ##
##                                            ## 
################################################


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
  variance = abs(isotopes[[4]]*0.001),
  samplesize = 1,
  time = isotopes[[5]]
)

carbon_isotopes_ts <- data.frame(
  mean = isotopes[[3]],
  variance = isotopes[[3]]*0.001,
  samplesize = 1,
  time = isotopes[[5]]
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



###### RUNNING Model_1a AND Model_1b ######
#Model_1a: empty A matrix; diagonal R matrix 
coccolith_model_1a <- fit.multivariate.URW(ln.bodylength_isotopes_multiobj, R = "diag", r = "fixed", iterations = 100, iter.sd = 0.000001, hess = TRUE)

#Model_1b: empty A matrix; full R matrix
coccolith_model_1b <- fit.multivariate.URW(ln.bodylength_isotopes_multiobj, R = "symmetric", r = "fixed", iterations = 100, iter.sd = 0.000001, hess = TRUE)

save.image(file <- 'coccolith_allruns_model_1ab.RData')

#calculating the correlation coefficients 
Rmatrix_best_coccolith_model_1b <- coccolith_model_1b$R
cormatrix_best_coccolith_model_1b <- cov2cor(Rmatrix_best_coccolith_model_1b)



###### EXTRACTING AND SAVING THE ITERATION WITH BEST PARAMETER ESTIMATES ######
print(coccolith_model_1a)
save(coccolith_model_1a, file = 'coccolith_bestrun_model_1a.RData')

print(coccolith_model_1b)
print(cormatrix_best_coccolith_model_1b)
save(coccolith_model_1b, cormatrix_best_coccolith_model_1b, file = 'coccolith_bestrun_model_1b.RData')