#real data boston housing example

#load data from MASS
library(MASS)
#help(Boston)
medv<- Boston$medv
RM <- Boston$rm
logTAX <- log(Boston$tax)
PTRATIO <- Boston$ptratio
logLSTAT <- log(Boston$lstat)
X <- cbind(RM,logTAX,PTRATIO,logLSTAT)
y0<-scale(medv)

est<-siqr(y0,X,p=0.5,maxiter = 60)





