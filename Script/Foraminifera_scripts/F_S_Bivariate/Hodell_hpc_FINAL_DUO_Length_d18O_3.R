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
data18O <- read.csv("/cluster/home/marionth/Hodell/table_d18O.csv")
datalength18O = cbind(data1,data18O)
data18Obis = datalength18O %>% filter(!is.na(Mean_isotopsd18O)) #enlever toutes les lignes contenant NA dans la colonne "Mean" pour data13C et datalength

sample_sized18O=data18Obis[[1]]
length_meansd18O=data18Obis[[2]]
length_varianced18O=data18Obis[[3]]
time_vectord18O=data18Obis[[4]]
d18O_means=data18Obis[[5]]
d18O_variance=data18Obis[[6]]
d18O_samplesize=data18Obis[[7]]
d18O_timevector=data18Obis[[8]]

ln.lengthd18O=paleoTS::as.paleoTS(mm=length_meansd18O, vv=length_varianced18O, nn=sample_sized18O, tt=time_vectord18O) 
ln.lengthd18O$tt=ln.lengthd18O$tt/(max(ln.lengthd18O$tt))
datad18Oi=paleoTS::as.paleoTS(mm=d18O_means, vv=d18O_variance, nn=d18O_samplesize, tt=d18O_timevector) #create two paleoTS object
datad18Oi$tt=datad18Oi$tt/(max(datad18Oi$tt))  #rescale du temps   


multid18Oi=make.multivar.evoTS(ln.lengthd18O, datad18Oi) 
plotevoTS.multivariate(multid18Oi, y_min=-3, y_max=3, x.label="relative time", y.label="Trait and proxy d13C means", pch=c(20,20))







#number of cores that I want to use (also defined in the batch script)
registerDoParallel(15)

#To calculate the time needed to process each iteration
start_time <- Sys.time()





#Loop foreach 
OUOUmodel.uptrid18O <- foreach(i=1:100) %dopar% { #calculation for each iteration 
  result <- multiResultClass() #create a table to store the results
  
  #set and register a different seed for each loop
  seed <- (10+i) 
  set.seed(seed)
  
  result$result1 <- seed
  result$result2 <- fit.multivariate.OU(multid18Oi,A.matrix ="upper.tri", R.matrix = "diag", trace=TRUE, iterations=1) 

  
  end_time <- Sys.time(); end_time - start_time 
  print( end_time - start_time ) 
  result$result3 <- end_time - start_time 
  
  return(result)
  
}


print ("DONE")
end_time <- Sys.time(); end_time - start_time


#creation of AIC matrix
AICvalues<-matrix(ncol=2,nrow=100)
colnames(AICvalues) <- c("iteration", "AIC")


#Loop to fill the AIC matrix
for (j in 1:100) {
  AICvalues[j,1] = j #iteration number
  AICvalues[j,2] = OUOUmodel.uptrid18O[[j]]$result2$AIC #AIC value
}


#Finding the best AIC
BestAICiteration = which.min(AICvalues[,2])

#Finding the best iteration associated to the best AIC
BESTRESULT=OUOUmodel.uptrid18O[[BestAICiteration]]
BESTRESULT



save.image(file='Hodell_Length_d18O_3.RData')
