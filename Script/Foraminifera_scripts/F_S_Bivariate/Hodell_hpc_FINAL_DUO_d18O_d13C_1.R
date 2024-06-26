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


d13C_means=data13Cbis[[5]]
d13C_variance=data13Cbis[[6]]
d13C_samplesize=data13Cbis[[7]]
d13C_timevector=data13Cbis[[8]]


datad13C=paleoTS::as.paleoTS(mm=d13C_means, vv=d13C_variance, nn=d13C_samplesize, tt=d13C_timevector)
datad13C$tt=datad13C$tt/(max(datad13C$tt)) #rescale le temps






data18O <- read.csv("/cluster/home/marionth/Hodell/table_d18O.csv")
datalength18O = cbind(data1,data18O)
data18Obis = datalength18O %>% filter(!is.na(Mean_isotopsd18O)) #enlever toutes les lignes contenant NA dans la colonne "Mean" pour data13C et datalength


d18O_means=data18Obis[[5]]
d18O_variance=data18Obis[[6]]
d18O_samplesize=data18Obis[[7]]
d18O_timevector=data18Obis[[8]]


datad18Oi=paleoTS::as.paleoTS(mm=d18O_means, vv=d18O_variance, nn=d18O_samplesize, tt=d18O_timevector) #create two paleoTS object
datad18Oi$tt=datad18Oi$tt/(max(datad18Oi$tt))  #rescale du temps   


multiENVI=make.multivar.evoTS(datad18Oi, datad13C) 
plotevoTS.multivariate(multiENVI, y_min=-3, y_max=3, x.label="relative time", y.label="Trait and proxy d13C means", pch=c(20,20))
 


##FIT D UNE UNBIASED RANDOM WALK
modelURWENVIdiag=fit.multivariate.URW(multiENVI, R="diag", r="fixed",iterations=100,trace=TRUE)





##ACCELERATED AND DECELERATED MODELS
modelaccelENVIdiag=fit.multivariate.URW(multiENVI, R="diag", r="accel",iterations=100,trace=TRUE)




modeldecelENVIdiag=fit.multivariate.URW(multiENVI, R="diag", r="decel",iterations=100,trace=TRUE)




save.image(file='Hodell_ENVI.RData')



modelURWENVIdiag
modelaccelENVIdiag
modeldecelENVIdiag