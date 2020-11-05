#real data boston housing example

#load data from MASS
library(MASS)
#help(Boston)
medv<- Boston$medv
RM <- Boston$rm
logTAX <- log(Boston$tax/10000)
PTRATIO <- Boston$ptratio
logLSTAT <- log(Boston$lstat)


X <- cbind(RM,logTAX,PTRATIO,logLSTAT)
y0<-medv - mean(medv)

#result is not the same with Wu 2010 as initial was not normalized in Wu 2010
gamma0<-c(1,-1,0,-1);
est<-siqr(y0,X,gamma.inital = gamma0, p=0.5,maxiter = 20,tol = 1e-6)




#asyp se



#Bootstrap

# input
n <- 200
data <- generate_data(n)
X <- data$X
Y <- data$Y

gamma0 <- c(1,2,3)
B <- 100
#residual bootstrap
fit <- siqr(Y, X, gamma.inital = NULL, p=0.5,maxiter = 30,tol = 1e-8)
res <- Y - fit$yhat

gammaB.matrix<-matrix(nrow=B,ncol=NCOL(X));
b=0
while(b<=B){
  B.index <- sample(n,replace=T)
  X_boot <- X
  Y_boot <- fit$yhat + res[B.index]
  gammaB.matrix[b,] <- siqr(Y_boot, X_boot, gamma.inital = gamma0, p=0.5,maxiter = 30,tol = 1e-8)$gamma
  b=b+1
}
apply(gammaB.matrix,2,function(x){sd(x,na.rm=TRUE)})



#XY bootstrap
gammaB.matrix<-matrix(nrow=B,ncol=NCOL(X));
b=0
while(b<=B){
  tryCatch({
    B.index <- sample(n,replace=T)
    X_boot <- X[B.index,]
    Y_boot <- Y[B.index]
    gammaB.matrix[b,] <- siqr(Y_boot, X_boot, gamma.inital = gamma0, p=0.5,maxiter = 30,tol = 1e-8)$gamma
    b=b+1
  },error = function(e) {
    print(e)
  })

}
apply(gammaB.matrix,2,sd)

