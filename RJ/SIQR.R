# This script includes codes used to produce real data analysis and simulation study
# in the R journal manuscript, titled
# SIQR: An R Package for Single-index Quantile Regression
#
# Authors: Tianhai Zu and Yan Yu.
# Date: Dec 7 2020
#
# The SIQR package proposed in the manuscript is available at github repository:
# https://github.com/zzz1990771/siqr .
# It can be installed through devtools:
#
#   devtools::install_github("zzz1990771/siqr")
#
# Or it can be installed manually:
#
#   install.packages("quantreg","KernSmooth")
#   install.packages(pkgs= "https://github.com/zzz1990771/siqr/raw/master/siqr_0.1.0.zip", repos = NULL, type = "win.binary")
#
# This package is built under R version 3.6.3.
#


library(siqr)
### Boston Housing data

#load data from MASS
library(MASS)
#help(Boston)
#data transformation
medv<- Boston$medv
RM <- Boston$rm
logTAX <- log(Boston$tax)
PTRATIO <- Boston$ptratio
logLSTAT <- log(Boston$lstat)
X <- cbind(RM,logTAX,PTRATIO,logLSTAT)
y0<-medv - mean(medv)

#initials
beta0 <- NULL

#fit the SIQR for different tau
tau.vec <- c(0.25,0.5,0.75)
est.coefficient <- matrix(NA, nrow = length(tau.vec), ncol = 5)
est.coefficient[,1] <- tau.vec
for (i in 1:length(tau.vec)){
  est <- siqr(y0,X,beta.inital = beta0, tau=tau.vec[i],maxiter = 20,tol = 1e-6)
  est.coefficient[i,2:5] <- est$beta
}
colnames(est.coefficient) <- c("quantile tau",colnames(X))
#print out estimated coefficients
est.coefficient

#The estimated 0.25 quantiles
#and their 95% point-wise confidence bounds
est.tau25 <- siqr(y0,X,beta.inital = NULL, tau=0.25)
plot.siqr(est.tau25,bootstrap.interval = TRUE)

#The estimated 0.50 quantiles
#and their 95% point-wise confidence bounds
est.tau50 <- siqr(y0,X,beta.inital = NULL, tau=0.5)
plot.siqr(est.tau50,bootstrap.interval = TRUE)

#The estimated 0.75 quantiles
#and their 95% point-wise confidence bounds
est.tau75 <- siqr(y0,X,beta.inital = NULL, tau=0.75)
plot.siqr(est.tau75,bootstrap.interval = TRUE)


### Simulation

#### Setting 1
n <- 400
beta0 <- c(1, 1, 1)/sqrt(3)
#the replication number set to very low for demo purpose.
n.sim <- 5
data <- generate.data(n, true.theta=beta0, setting = "setting1",ncopy = n.sim)

#for demo only, we used large n.sim in paper
library(foreach)
sim.results.50 <- foreach(m = 1:n.sim,.combine = "rbind") %do% {
  X <- data$X
  Y <- data$Y[[m]]
  est <- siqr(Y, X, beta.inital = c(2,1,0), tau=0.5,maxiter = 30,tol = 1e-8)
  return(est$beta)
}

sim.results.25 <- foreach(m = 1:n.sim,.combine = "rbind") %do% {
  X <- data$X
  Y <- data$Y[[m]]
  est <- siqr(Y, X, beta.inital = c(2,1,0), tau=0.25,maxiter = 30,tol = 1e-8)
  return(est$beta)
}
sim.results.75 <- foreach(m = 1:n.sim,.combine = "rbind") %do% {
  X <- data$X
  Y <- data$Y[[m]]
  est <- siqr(Y, X, beta.inital = c(2,1,0), tau=0.75,maxiter = 30,tol = 1e-8)
  return(est$beta)
}

# obtain a box plot ofestimated single-index coefficients for each tau
# note that you may get some warning message from boxplot when using low
# number of replication
boxplot(data.frame((sim.results.25)), outline=T,notch=T,range=1,main = "Boxplots of Coefficient Estimates, tau = 0.25",horizontal = F)

boxplot(data.frame((sim.results.50)), outline=T,notch=T,range=1,main = "Boxplots of Coefficient , tau = 0.50",horizontal = F)

boxplot(data.frame((sim.results.75)), outline=T,notch=T,range=1,main = "Boxplots of Coefficient Estimates, tau = 0.75",horizontal = F)

#### Setting 3
n <- 400
beta0 <- c(1, 2)/sqrt(5)
#the replication number set to very low for demo purpose.
n.sim <- 5
tau <- 0.5
data <- generate.data(n, true.theta=beta0, setting = "setting3",ncopy = n.sim)

sim.results <- foreach(m = 1:n.sim,.combine = "rbind") %do% {
  X <- data$X
  Y <- data$Y[[m]]
  est <- siqr(Y, X, beta.inital = NULL, tau=tau,maxiter = 30,tol = 1e-8)
  est$beta
}

#print estimated coefficients
est.mean <- c(tau,apply(sim.results,2,mean))
names(est.mean) <- c("tau","beta1.hat","beta2.hat")
est.mean

#print the simulation standard error of estimated coefficients
est.se <- c(tau,apply(sim.results,2,sd))
names(est.se) <- c("tau","beta1.se.hat","beta1.se.hat")
est.se


# obtain a box plot ofestimated single-index coefficients for tau = 0.5
boxplot(data.frame((sim.results)), outline=T,notch=T,range=1,main = "Boxplots of Coefficient Estimates, Example 2",horizontal = F)


# plot the estimated quantiles and their 95% point-wise confidence bounds
est.sim.05 <- siqr(data$Y[[1]],data$X,beta.inital = NULL, tau=0.5)
plot.siqr(est.sim.05,bootstrap.interval = TRUE)




