####################
# BORNEMANN SCRIPT #
####################


.libPaths("/cluster/home/marionth/R")
library(usethis)
library(evoTS)
library(paleoTS)



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



## Length ratio of R unit to placolith

data2 <- read.delim("/cluster/home/marionth/Bornemann/Biscutum_RatioRcentralunitTotallength.txt")
summary(data2)

lengthratio_means=data2[[2]]
lengthratio_variance=data2[[3]]

lengthratio=paleoTS::as.paleoTS(mm=lengthratio_means, vv=lengthratio_variance, nn=sample_size, tt=time_vector) # faire un paleoTS object
ln.lengthratio=paleoTS::ln.paleoTS(lengthratio) # creer une approximative log-transformation des donnees si ce n est pas deja fait
ln.lengthratio$tt<-ln.lengthratio$tt/(max(ln.lengthratio$tt)) # transferer time vector en une unite de longueur

plotevoTS(ln.lengthratio) # plotter l evolution du trait au cours du temps



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




multienvi=make.multivar.evoTS(d18Ocarb, d13Ccarb)





## WITH ENVIRONMENTAL FACTOR AFFECTING THE OPTIMUM OF PHENOTYPIC TRAITS FOR NUMBER OF ENVIRONMENTAL VARIABLES REDUCED (2)
A = matrix(c(0,1,0,1), nrow = 2, byrow=TRUE) #1 for parametrised and 0 for not parametrised
R = matrix(c(1,0,0,1), nrow = 2, byrow=TRUE) #5 traits donc 25 cases dans la matrice
OUOUmodel.fulldiagENVI=fit.multivariate.OU.user.defined(multienvi,A.user=A,R.user=R, iterations=100)
# AICc = 







save.image(file='Bornemann_FINAL_d18O_d13C_8.RData')


OUOUmodel.fulldiagENVI

