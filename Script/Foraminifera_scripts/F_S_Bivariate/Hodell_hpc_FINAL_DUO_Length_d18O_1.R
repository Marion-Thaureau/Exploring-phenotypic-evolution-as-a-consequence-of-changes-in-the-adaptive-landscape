###### PACKAGES AND FILE DIRECTORY ######

.libPaths("/cluster/home/marionth/R")
library(evoTS)
library(paleoTS)
library(stats)
library(dplyr)


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






##FIT D UNE UNBIASED RANDOM WALK
modelURWd18Odiag=fit.multivariate.URW(multid18Oi, R="diag", r="fixed",iterations=100,trace=TRUE)







##ACCELERATED AND DECELERATED MODELS
modelacceld18Odiag=fit.multivariate.URW(multid18Oi, R="diag", r="accel",iterations=100,trace=TRUE)




modeldeceld18Odiag=fit.multivariate.URW(multid18Oi, R="diag", r="decel",iterations=100,trace=TRUE)







save.image(file='Hodell_Length_d18O_1.RData')





modelURWd18Odiag

modelacceld18Odiag
modeldeceld18Odiag
