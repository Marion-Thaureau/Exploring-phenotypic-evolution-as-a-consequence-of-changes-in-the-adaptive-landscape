.libPaths("/cluster/home/marionth/R")
library(evoTS)
library(paleoTS)
library(stats)
library(dplyr)
library(parallel)
library(foreach)
library(doParallel)
library(MCSPCD)


setwd("/cluster/home/marionth/Jensen")



###### LOAD MORPHOLOGICAL DATA ######



data1=read.delim("/cluster/home/marionth/Jensen/Atlantic_salmon_size.txt")
data1 = data1 %>% filter(!is.na(Year.of.catch))

year=unique(data1$Year.of.catch)  #garder qu un exemplaire de chaque date
data2=matrix(nrow=length(year), ncol=4)
colnames(data2)<-c("nn"," mm", "vv", "tt")

for (i in 1:length(year)){
  subdat=subset(data1, data1$Year.of.catch == year[i])  # un subset par annee
  
  data2[i,1] = nrow(subdat)
  data2[i,2] = mean(subdat$Mass)
  data2[i,3] = var(subdat$Mass)
  data2[i,4] = year[i]
}

data2=data2[order(year),]



###### FOR FROM JUN TO SEPT WATER FLOW VALUES ######
###### LOAD ENVIRONMENTAL DATA AND RESIZE ######



Discharge <- read.delim("/cluster/home/marionth/Jensen/Discharge.txt")
Envi = Discharge[-c(1:9,28:35,87:88),]
Morpho = data2[-c(1:2),]  # garder que les lignes ou on a les donnees pour l environment et le morpho



###### CREATING PALEOTS AND EVOTS OBJECTS ######



dataEnvi=paleoTS::as.paleoTS(mm=Envi[,3], vv=c(rep(0.000001,length(Envi[,1]))), nn=c(rep(1000,length(Envi[,1]))), tt=Envi[,1], oldest = "first")
ln.dataEnvi=paleoTS::ln.paleoTS(dataEnvi)
dataEnvi$tt=dataEnvi$tt/(max(dataEnvi$tt))
ln.dataEnvi$tt=ln.dataEnvi$tt/(max(ln.dataEnvi$tt))

dataMorpho=paleoTS::as.paleoTS(mm=Morpho[,2], vv=Morpho[,3], nn=Morpho[,1], tt=Morpho[,4], oldest = "first")
ln.dataMorpho=paleoTS::ln.paleoTS(dataMorpho)
dataMorpho$tt=dataMorpho$tt/(max(dataMorpho$tt))
ln.dataMorpho$tt=ln.dataMorpho$tt/(max(ln.dataMorpho$tt))





multilnmorpholnenvi=make.multivar.evoTS(ln.dataMorpho, ln.dataEnvi) 
plotevoTS.multivariate(multilnmorpholnenvi, y_min=0, y_max=20, x.label="relative time", y.label="Salmon mass and water flow", pch=c(20,20))





###### MULTIVARIATE MODELS ######




#number of cores that I want to use (also defined in the batch script)
registerDoParallel(15)

#To calculate the time needed to process each iteration
start_time <- Sys.time()



A = matrix(c(1,0,0,1), nrow = 2, byrow=TRUE) #1 for parametrised and 0 for not parametrised
R = matrix(c(1,1,1,1), nrow = 2, byrow=TRUE) #5 traits donc 25 cases dans la matrice





#Loop foreach 
OUOUmodel.upperrowdiagENVIjs <- foreach(i=1:100) %dopar% { #calculation for each iteration 
  result <- multiResultClass() #create a table to store the results
  
  #set and register a different seed for each loop
  seed <- (10+i) 
  set.seed(seed)
  
  result$result1 <- seed
  result$result2 <- fit.multivariate.OU.user.defined(multilnmorpholnenvi,A.user=A,R.user=R, iterations=1, iter.sd=0.000001)
  
  
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
  AICvalues[j,2] = OUOUmodel.upperrowdiagENVIjs[[j]]$result2$AIC #AIC value
}


#Finding the best AIC
BestAICiteration = which.min(AICvalues[,2])

#Finding the best iteration associated to the best AIC
BESTRESULT=OUOUmodel.upperrowdiagENVIjs[[BestAICiteration]]
BESTRESULT




save.image(file='Jensen_FINAL_junsept_5.RData')
