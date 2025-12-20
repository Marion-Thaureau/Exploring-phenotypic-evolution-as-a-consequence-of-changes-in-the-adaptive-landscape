#############################################
##                                         ##
##  SALMON SCRIPT - Model 1a and Model 1b  ##
##                                         ## 
#############################################


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


setwd("/cluster/home/marionth/Jensen")



###### LOAD MORPHOLOGICAL DATA ######
bodymass_raw <- read.delim("/cluster/home/marionth/Jensen/Atlantic_salmon_size.txt")

#Cleaning the data and log-transforming the data
bodymass_raw <- bodymass_raw %>% filter(!is.na(Year.of.catch)) #filter out NA rows
ln.bodymass_raw <- within(bodymass_raw, Mass <- log(Mass))

#Calculating the mean and variance of the observations for each year
year <- unique(ln.bodymass_raw$Year.of.catch)
ln.bodymass <- data.frame()

for (i in 1:length(year)) {
  year_subset <- subset(ln.bodymass_raw, ln.bodymass_raw$Year.of.catch == year[i]) #Create subsets based on the year
  year_stats <- data.frame(samplesize = nrow(year_subset), #calculate the mean and the variance for each subset
                        mean = mean(year_subset$Mass, na.rm = TRUE),  
                        variance = var(year_subset$Mass, na.rm = TRUE),  
                        time = year[i])
  ln.bodymass <- rbind(ln.bodymass, year_stats)
}

ln.bodymass <- ln.bodymass[order(ln.bodymass$time), ] #ordering chronologically the dataset



###### LOAD ENVIRONMENTAL DATA ######
discharge_raw <- read.delim("/cluster/home/marionth/Jensen/Discharge.txt")
discharge <- subset(discharge_raw, select = -DC_year) #we keep the discharge measured between June and September, not the average over the complete year

#Keeping the years common to both bodymass record and discharge
years_ts <- intersect(ln.bodymass$time, discharge$Year)
ln.bodymass_ts <- ln.bodymass[ln.bodymass$time %in% years_ts, ]
discharge_ts <- discharge[discharge$Year %in% years_ts, ]

#Implementation of variance and sample size 
discharge_ts$variance <- discharge_ts$DC_Jun_Sept*0.05
discharge_ts$samplesize <- 1


  
###### CREATING PALEOTS AND EVOTS OBJECTS ######
#Creating paleoTS object for the bodymass time series
ln.bodymass_obj <- paleoTS::as.paleoTS(mm = ln.bodymass_ts$mean, vv = ln.bodymass_ts$variance, nn = ln.bodymass_ts$samplesize, tt = ln.bodymass_ts$time, oldest = "first")
ln.bodymass_obj$tt <- ln.bodymass_obj$tt/(max(ln.bodymass_obj$tt)) #transfer of the time vector into a length unit of 1 

#Creating paleoTS object for the waterflow time series
discharge_obj <- paleoTS::as.paleoTS(mm = discharge_ts$DC_Jun_Sept, vv = discharge_ts$variance, nn = discharge_ts$samplesize, tt = discharge_ts$Year, oldest = "first")
ln.discharge_obj <- paleoTS::ln.paleoTS(discharge_obj) #log-transformation of the data
ln.discharge_obj$tt <- ln.discharge_obj$tt/(max(ln.discharge_obj$tt)) #transfer of the time vector into a length unit of 1 

#Creating and plotting the multivariate evoTS object
ln.bodymass_ln.discharge_multiobj <- make.multivar.evoTS(ln.bodymass_obj,ln.discharge_obj)



###### RUNNING Model_1a AND Model_1b ######
#Model_1a: empty A matrix; diagonal R matrix 
salmon_model_1a <- fit.multivariate.URW(ln.bodymass_ln.discharge_multiobj, R = "diag", r = "fixed", iterations = 100, iter.sd = 0.000001, hess = TRUE)

#Model_1b: empty A matrix; full R matrix
salmon_model_1b <- fit.multivariate.URW(ln.bodymass_ln.discharge_multiobj, R = "symmetric", r = "fixed", iterations = 100, iter.sd = 0.000001, hess = TRUE)

save.image(file <- 'salmon_allruns_model_1ab.RData')

#calculating the correlation coefficients 
Rmatrix_best_salmon_model_1b <- salmon_model_1b$R
cormatrix_best_salmon_model_1b <- cov2cor(Rmatrix_best_salmon_model_1b)



###### EXTRACTING AND SAVING THE ITERATION WITH BEST PARAMETER ESTIMATES ######
print(salmon_model_1a)
save(salmon_model_1a, file = 'salmon_bestrun_model_1a.RData')

print(salmon_model_1b)
print(cormatrix_best_salmon_model_1b)
save(salmon_model_1b, cormatrix_best_salmon_model_1b, file = 'salmon_bestrun_model_1b.RData')
