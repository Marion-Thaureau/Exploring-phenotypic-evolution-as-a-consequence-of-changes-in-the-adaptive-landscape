#R version 4.2.1
.libPaths("/cluster/home/marionth/R")
library(evoTS)

setwd("/cluster/home/marionth/Detrending")


nrep <- 100
winners <- character(nrep)   # store best model name each run
t1 <- seq(1, 30, by = 1)
tt <- cbind(t1, t1)
vv <- cbind(rep(1, 30), rep(1, 30))
nn <- cbind(rep(30, 30), rep(30, 30))
set.seed(123)
i <- 0          ###
attempt <- 0    ###
while (i < nrep) {   ###
  attempt <- attempt + 1   ###
  resid1 <- rnorm(30, mean = 0, sd = 1.5)
  x1 <- (1 + .5 * t1 + resid1) #variable 1 is assosciated with time, but not directly affecting variable 2
  resid2 <- rnorm(30, mean = 0, sd = 2)
  x2 <- (1 + .3 * t1 + resid2) #variable 2 is assosciated with time, but not directly affecting variable 1
  xx <- cbind(x1, x2)

  x <- list(xx = xx, vv = vv, nn = nn, tt = tt)
  plotevoTS.multivariate(x)
  plot(x2,x1)

  x$tt[,1]<-x$tt[,1]/max(x$tt[,1])
  x$tt[,2]<-x$tt[,1]

  # ---- Fit models ----
  m1<-fit.multivariate.URW(x, R="diag")
  print(m1$converge)
  if (m1$converge != "Model converged successfully") { message("Attempt ", attempt, ": m1 did not converge, re-simulating."); next }   ###
  m2<-fit.multivariate.URW(x)
  print(m2$converge)
  if (m2$converge != "Model converged successfully") { message("Attempt ", attempt, ": m2 did not converge, re-simulating."); next }   ###
  m3 <- fit.multivariate.OU(x, A.matrix="diag", R.matrix="diag")
  print(m3$converge)
  if (m3$converge != "Model converged successfully") { message("Attempt ", attempt, ": m3 did not converge, re-simulating."); next }   ###
  m4 <- fit.multivariate.OU(x, A.matrix="diag", R.matrix="symmetric")
  print(m4$converge)
  if (m4$converge != "Model converged successfully") { message("Attempt ", attempt, ": m4 did not converge, re-simulating."); next }   ###
  m5 <- fit.multivariate.OU(x, A.matrix="upper.tri", R.matrix="diag")
  print(m5$converge)
  if (m5$converge != "Model converged successfully") { message("Attempt ", attempt, ": m5 did not converge, re-simulating."); next }   ###
  m6 <- fit.multivariate.OU(x, A.matrix="upper.tri", R.matrix="symmetric")
  print(m6$converge)
  if (m6$converge != "Model converged successfully") { message("Attempt ", attempt, ": m6 did not converge, re-simulating."); next }   ###

  i <- i + 1   ###

  # ---- Extract AICcs ----
  aics <- c(
    m1 = m1$AICc,
    m2 = m2$AICc,
    m3 = m3$AICc,
    m4 = m4$AICc,
    m5 = m5$AICc,
    m6 = m6$AICc
  )

  print(i)
  print(aics)

  # ---- Store best model ----
  winners_s1[i] <- names(which.min(aics))
}

# ---- Count how often each model wins ----
table(winners_s1)

save(winners_s1, file = "reviewer1_revised.RData")
