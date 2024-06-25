###### PACKAGES AND FILE DIRECTORY ######

.libPaths("/cluster/home/marionth/R")
library(evoTS)
library(paleoTS)
library(stats)
library(dplyr)
library(parallel)
library(foreach)
library(doParallel)
library(MCSPCD)


setwd("/cluster/home/marionth/Hodell")
getwd()




###### LENGTH WITH d18O PROXY ITERATIONS ######



data1=read.table("Hodell_and_Vayavananda_1993_Fohsella_length.txt",header=T)
data13C <- read.csv("/cluster/home/marionth/Hodell/table_d13C.csv")
datalength13C = cbind(data1,data13C)
data13Cbis = datalength13C %>% filter(!is.na(Mean_isotopsd13C)) #enlever toutes les lignes contenant NA dans la colonne "Mean" pour data13C et datalength




sample_sized13C=data13Cbis[[1]]
length_meansd13C=data13Cbis[[2]]
length_varianced13C=data13Cbis[[3]]
time_vectord13C=data13Cbis[[4]]


d13C_means=data13Cbis[[5]]
d13C_variance=data13Cbis[[6]]
d13C_samplesize=data13Cbis[[7]]
d13C_timevector=data13Cbis[[8]]


datad13C=paleoTS::as.paleoTS(mm=d13C_means, vv=d13C_variance, nn=d13C_samplesize, tt=d13C_timevector)
#pas log transformed
datad13C$tt=datad13C$tt/(max(datad13C$tt)) #rescale le temps

ln.lengthd13C=paleoTS::as.paleoTS(mm=length_meansd13C, vv=length_varianced13C, nn=sample_sized13C, tt=time_vectord13C)
ln.lengthd13C$tt=ln.lengthd13C$tt/(max(ln.lengthd13C$tt))



data18O <- read.csv("/cluster/home/marionth/Hodell/table_d18O.csv")
datalength18O = cbind(data1,data18O)
data18Obis = datalength18O %>% filter(!is.na(Mean_isotopsd18O)) #enlever toutes les lignes contenant NA dans la colonne "Mean" pour data13C et datalength


d18O_means=data18Obis[[5]]
d18O_variance=data18Obis[[6]]
d18O_samplesize=data18Obis[[7]]
d18O_timevector=data18Obis[[8]]


datad18Oi=paleoTS::as.paleoTS(mm=d18O_means, vv=d18O_variance, nn=d18O_samplesize, tt=d18O_timevector) #create two paleoTS object
#pas log transformed
datad18Oi$tt=datad18Oi$tt/(max(datad18Oi$tt))  #rescale du temps   


multiENVI=make.multivar.evoTS(ln.lengthd13C, datad18Oi, datad13C) 
plotevoTS.multivariate(multiENVI, y_min=-3, y_max=3, x.label="relative time", y.label="Trait and proxy d18O and d13C means", pch=c(20,20))





#number of cores that I want to use (also defined in the batch srcipt)
registerDoParallel(52)

#To calculate the time needed to process each iteration
start_time <- Sys.time()



##NO COEVOLUTION
A = matrix(c(1,1,0,0,1,0,0,0,1), nrow = 3, byrow=TRUE)
R = matrix(c(1,0,0,0,1,0,0,0,1), nrow = 3, byrow=TRUE) 




#Loop foreach 
OUOU.model1 <- foreach(i=1:100) %dopar% { #calculation for each iteration 
  result <- multiResultClass() #create a table to store the results
  
  #set and register a different seed for each loop
  seed <- (10+i) 
  set.seed(seed)
  
  result$result1 <- seed
  result$result2 <- fit.multivariate.OU.user.defined(multiENVI, A.user=A,R.user=R, iterations = 1, hess = TRUE, 
                                                     user.init.diag.A = c(29.67556,0.2847143,17.772),
						     user.init.upper.diag.A = c(-25.6442576))
  
  end_time <- Sys.time(); end_time - start_time 
  print( end_time - start_time ) 
  result$result3 <- end_time - start_time 
  
  return(result)
  
}


print ("DONE")
end_time <- Sys.time(); end_time - start_time





save.image(file='Hodell_quat_M6.RData')





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
