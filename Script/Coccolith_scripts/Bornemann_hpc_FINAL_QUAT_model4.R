####################
# BORNEMANN SCRIPT #
####################


.libPaths("/cluster/home/marionth/R")
library(usethis)
library(evoTS)
library(paleoTS)
library(stats)
library(dplyr)
library(parallel)
library(foreach)
library(doParallel)
library(MCSPCD)


setwd("/cluster/home/marionth/Bornemann")




# LOAD MORPHOLOGICAL DATA
## Coccolith ellipsicity

data1 <- read.delim("/cluster/home/marionth/Bornemann/Biscutum_EllipsicitePlacolith.txt")
summary(data1)

sample_size=data1[[1]]
ellipsicite_means=data1[[2]]
ellipsicite_variance=data1[[3]]
time_vector=data1[[4]]

ellipsicite=paleoTS::as.paleoTS(mm=ellipsicite_means, vv=ellipsicite_variance, nn=sample_size, tt=time_vector) # faire un paleoTS object
ln.ellipsicite=paleoTS::ln.paleoTS(ellipsicite) # creer une approximative log-transformation des donnees si ce n est pas deja fait
ln.ellipsicite$tt<-ln.ellipsicite$tt/(max(ln.ellipsicite$tt)) # transferer time vector en une unite de longueur

plotevoTS(ln.ellipsicite) # plotter l evolution du trait au cours du temps



## Size of central unit

data2 <- read.delim("/cluster/home/marionth/Bornemann/Biscutum_LengthRcentralunit.txt")
summary(data2)

lengthcentral_means=data2[[2]]
lengthcentral_variance=data2[[3]]

lengthcentral=paleoTS::as.paleoTS(mm=lengthcentral_means, vv=lengthcentral_variance, nn=sample_size, tt=time_vector) # faire un paleoTS object
ln.lengthcentral=paleoTS::ln.paleoTS(lengthcentral) # creer une approximative log-transformation des donnees si ce n est pas deja fait
ln.lengthcentral$tt<-ln.lengthcentral$tt/(max(ln.lengthcentral$tt)) # transferer time vector en une unite de longueur

plotevoTS(ln.lengthcentral) # plotter l evolution du trait au cours du temps




## Total length of the coccolith

data3 <- read.delim("/cluster/home/marionth/Bornemann/Biscutum_TotalcoccolithLength.txt")
summary(data3)

totallength_means=data3[[2]]
totallength_variance=data3[[3]]

totallength=paleoTS::as.paleoTS(mm=totallength_means, vv=totallength_variance, nn=sample_size, tt=time_vector) # faire un paleoTS object
ln.totallength=paleoTS::ln.paleoTS(totallength) # creer une approximative log-transformation des donnees si ce n est pas deja fait
ln.totallength$tt<-ln.totallength$tt/(max(ln.totallength$tt)) # transferer time vector en une unite de longueur

plotevoTS(ln.totallength) # plotter l evolution du trait au cours du temps





# LOAD ENVIRONMENTAL DATA




isotopes <- read.delim("/cluster/home/marionth/Bornemann/isotopes.txt")

geochemistrysample_size=isotopes[[1]]
geochemistrytime_vector=isotopes[[7]]
geochemistryvariance=isotopes[[2]]
d13Ccarb_means=isotopes[[5]]
d18Ocarb_means=isotopes[[6]]



d13Ccarb=paleoTS::as.paleoTS(mm=d13Ccarb_means, vv=geochemistryvariance+0.0001, nn=geochemistrysample_size, tt=geochemistrytime_vector, oldest="first")
d13Ccarb$tt=d13Ccarb$tt/(max(d13Ccarb$tt))

d18Ocarb=paleoTS::as.paleoTS(mm=d18Ocarb_means, vv=geochemistryvariance+0.0001, nn=geochemistrysample_size, tt=geochemistrytime_vector, oldest="first")
d18Ocarb$tt=d18Ocarb$tt/(max(d18Ocarb$tt))






multiENVI=make.multivar.evoTS(ln.totallength, d18Ocarb, d13Ccarb) 
plotevoTS.multivariate(multiENVI, y_min=-3, y_max=3, x.label="relative time", y.label="Trait and proxy d18O and d13C means", pch=c(20,20))





#number of cores that I want to use (also defined in the batch srcipt)
registerDoParallel(52)

#To calculate the time needed to process each iteration
start_time <- Sys.time()



##NO COEVOLUTION
A = matrix(c(1,0,1,0,1,0,0,0,1), nrow = 3, byrow=TRUE)
R = matrix(c(1,0,0,0,1,0,0,0,1), nrow = 3, byrow=TRUE) 




#Loop foreach 
OUOU.model1 <- foreach(i=1:100) %dopar% { #calculation for each iteration 
  result <- multiResultClass() #create a table to store the results
  
  #set and register a different seed for each loop
  seed <- (10+i) 
  set.seed(seed)
  
  result$result1 <- seed
  result$result2 <- fit.multivariate.OU.user.defined(multiENVI, A.user=A,R.user=R, iterations = 1, hess = TRUE, 
                                                     user.init.diag.A = c(186.618,41.80156,1.14631),
                                                     user.init.upper.diag.A = c(25.08519))
  
  end_time <- Sys.time(); end_time - start_time 
  print( end_time - start_time ) 
  result$result3 <- end_time - start_time 
  
  return(result)
  
}


print ("DONE")
end_time <- Sys.time(); end_time - start_time





save.image(file='Bornemann_FINAL_QUAT_model4.RData')





#creation of AIC matrix
AICvalues<-matrix(ncol=2,nrow=100)
colnames(AICvalues) <- c("iteration", "AIC")


#Loop to fill the AIC matrix
for (j in 1:100) {
  AICvalues[j,1] = j #iteration number
  AICvalues[j,2] = OUOU.model1[[j]]$result2$AIC #AIC value
}


#Finding the best AIC
BestAICiteration = which.min(AICvalues[,2])

#Finding the best iteration associated to the best AIC
BESTRESULT=OUOU.model1[[BestAICiteration]]
BESTRESULT
