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



##FIT D UNE UNBIASED RANDOM WALK
modelURWlndiag=fit.multivariate.URW(multilnmorpholnenvi, R="diag", r="fixed",trace=TRUE, iterations=100, iter.sd=0.000001)




##ACCELERATED AND DECELERATED MODELS
modelaccellndiag=fit.multivariate.URW(multilnmorpholnenvi, R="diag", r="accel",trace=TRUE, iterations=100, iter.sd=0.000001)


modeldecellndiag=fit.multivariate.URW(multilnmorpholnenvi, R="diag", r="decel",trace=TRUE, iterations=100, iter.sd=0.000001)








save.image(file='Jensen_FINAL_junsept_1.RData')




modelURWlndiag


modelaccellndiag

modeldecellndiag
