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

ln.lengthd13C=paleoTS::as.paleoTS(mm=length_meansd13C, vv=length_varianced13C, nn=sample_sized13C, tt=time_vectord13C) 
ln.lengthd13C$tt=ln.lengthd13C$tt/(max(ln.lengthd13C$tt))
datad13C=paleoTS::as.paleoTS(mm=d13C_means, vv=d13C_variance, nn=d13C_samplesize, tt=d13C_timevector)
datad13C$tt=datad13C$tt/(max(datad13C$tt)) #rescale le temps


multid13Ci=make.multivar.evoTS(ln.lengthd13C, datad13C) 
plotevoTS.multivariate(multid13Ci, y_min=-3, y_max=3, x.label="relative time", y.label="Trait and proxy d13C means", pch=c(20,20))





## WITH INDEPENDANT EVOLUTION (1)
OUOUmodel.indepd13C=fit.multivariate.OU(multid13Ci,A.matrix ="diag", R.matrix = "diag", trace=TRUE, iterations=100, iter.sd=0.000001) 
# AICc = 







save.image(file='Hodell_Length_d13C_2sd.RData')





OUOUmodel.indepd13C
