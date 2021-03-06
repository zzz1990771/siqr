## ----eval = TRUE--------------------------------------------------------------
# install.packages("quantreg","KernSmooth")
# install.packages(pkgs= "https://github.com/zzz1990771/siqr/raw/master/siqr_0.1.0.zip", repos = NULL, type = "win.binary")
library(siqr)

#load data from MASS
library(MASS)
#help(Boston)
medv<- Boston$medv
RM <- Boston$rm
logTAX <- log(Boston$tax)
PTRATIO <- Boston$ptratio
logLSTAT <- log(Boston$lstat)

X <- cbind(RM,logTAX,PTRATIO,logLSTAT)
y0<-medv - mean(medv)

beta0 <- NULL
tau.vec <- c(0.25,0.5,0.75)
est.coefficient <- matrix(NA, nrow = length(tau.vec), ncol = 5)
est.coefficient[,1] <- tau.vec
for (i in 1:length(tau.vec)){
  est <- siqr(y0,X,beta.initial = beta0, tau=tau.vec[i],maxiter = 20,tol = 1e-6)
  est.coefficient[i,2:5] <- est$beta
}
colnames(est.coefficient) <- c("quantile tau",colnames(X))
est.coefficient



## ----cache=TRUE, eval = TRUE--------------------------------------------------
est.tau25 <- siqr(y0,X,beta.initial = NULL, tau=0.25)
plot.siqr(est.tau25,bootstrap.interval = TRUE)


## ----cache=TRUE, eval = TRUE--------------------------------------------------
est.tau50 <- siqr(y0,X,beta.initial = NULL, tau=0.5)
plot.siqr(est.tau50,bootstrap.interval = TRUE)


## ----cache=TRUE, eval = TRUE--------------------------------------------------
est.tau75 <- siqr(y0,X,beta.initial = NULL, tau=0.75)
plot.siqr(est.tau75,bootstrap.interval = TRUE)


## ----eval=FALSE---------------------------------------------------------------
#> n <- 400
#> beta0 <- c(1, 1, 1)/sqrt(3)
#> n.sim <- 200
#> tau.vec <- c(0.25,0.5,0.75)
#> tau <- tau.vec[1]
#> 
#> data <- generate.data(n, true.theta=beta0, setting = "setting1",ncopy = n.sim)
#> 
#> #paralell
#> library(parallel)
#> library(foreach)
#> cl<- makeCluster(12)
#> doParallel::registerDoParallel(cl)
#> sim.results.50 <- foreach(m = 1:n.sim,.combine = "rbind") %dopar% {
#>   X <- data$X
#>   Y <- data$Y[[m]]
#>   est <- siqr(Y, X, beta.initial = c(2,1,0), tau=0.5,maxiter = 30,tol = 1e-8)
#>   if(est$flag.conv == 0){
#>     return(NULL)
#>   }else{
#>     return(est$beta)
#>   }
#> }
#> 
#> sim.results.25 <- foreach(m = 1:n.sim,.combine = "rbind") %dopar% {
#>   X <- data$X
#>   Y <- data$Y[[m]]
#>   est <- siqr(Y, X, beta.initial = c(2,1,0), tau=0.25,maxiter = 30,tol = 1e-8)
#>   if(est$flag.conv == 0){
#>     return(NULL)
#>   }else{
#>     return(est$beta)
#>   }
#> }
#> sim.results.75 <- foreach(m = 1:n.sim,.combine = "rbind") %dopar% {
#>   X <- data$X
#>   Y <- data$Y[[m]]
#>   est <- siqr(Y, X, beta.initial = c(2,1,0), tau=0.75,maxiter = 30,tol = 1e-8)
#>   if(est$flag.conv == 0){
#>     return(NULL)
#>   }else{
#>     return(est$beta)
#>   }
#> }
#> stopCluster(cl)


## -----------------------------------------------------------------------------
sim.results.25 <- readRDS("./sim1.results25.RDS")
sim.results.50 <- readRDS("./sim1.results50.RDS")
sim.results.75 <- readRDS("./sim1.results75.RDS")


## -----------------------------------------------------------------------------
boxplot(data.frame((sim.results.25)), outline=T,notch=T,range=1,main = "Boxplots of Coefficient Estimates, tau = 0.25",horizontal = F,
names=c(expression(hat(beta)[1]),expression(hat(beta)[2]),expression(hat(beta)[3])))

## -----------------------------------------------------------------------------
boxplot(data.frame((sim.results.50)), outline=T,notch=T,range=1,main = "Boxplots of Coefficient Estimates, tau = 0.25",horizontal = F,
names=c(expression(hat(beta)[1]),expression(hat(beta)[2]),expression(hat(beta)[3])))


## -----------------------------------------------------------------------------
boxplot(data.frame((sim.results.75)), outline=T,notch=T,range=1,main = "Boxplots of Coefficient Estimates, tau = 0.25",horizontal = F,
names=c(expression(hat(beta)[1]),expression(hat(beta)[2]),expression(hat(beta)[3])))

## ----cache=TRUE, eval = FALSE-------------------------------------------------
#> est.sim.05 <- siqr(data$Y[[1]],data$X,beta.initial = NULL, tau=0.5)
#> plot.siqr(est.sim.05,bootstrap.interval = TRUE)


## ----eval=FALSE---------------------------------------------------------------
#> n <- 400
#> beta0 <- c(1, 2)/sqrt(5)
#> n.sim <- 100
#> tau <- 0.5
#> 
#> data <- generate.data(n, true.theta=beta0, setting = "setting3",ncopy = n.sim)
#> 
#> #paralell
#> library(parallel)
#> library(foreach)
#> cl<- makeCluster(12)
#> doParallel::registerDoParallel(cl)
#> sim.results <- foreach(m = 1:n.sim,.combine = "rbind") %dopar% {
#>   X <- data$X
#>   Y <- data$Y[[m]]
#>   est <- siqr(Y, X, beta.initial = NULL, tau=tau,maxiter = 30,tol = 1e-8)
#>   est$beta
#> }


## -----------------------------------------------------------------------------
tau <- 0.5
sim.results <- readRDS("./sim.results.RDS")
est.mean <- c(tau,apply(sim.results,2,mean))
names(est.mean) <- c("tau","beta1.hat","beta2.hat")
est.mean


## -----------------------------------------------------------------------------
est.se <- c(tau,apply(sim.results,2,sd))
names(est.se) <- c("tau","beta1.se.hat","beta1.se.hat")
est.se


## -----------------------------------------------------------------------------
boxplot(data.frame((sim.results)), outline=T,notch=T,range=1,main = "Boxplots of Coefficient Estimates, Example 2",horizontal = F,
names=c(expression(hat(beta)[1]),expression(hat(beta)[2])))


## ----cache=TRUE, eval = TRUE--------------------------------------------------
n <- 400
beta0 <- c(1, 2)/sqrt(5)
n.sim <- 100
tau <- 0.5
data <- generate.data(n, true.theta=beta0, setting = "setting3",ncopy = 2)
est.sim.05 <- siqr(data$Y[[1]],data$X,beta.initial = NULL, tau=0.5)
plot.siqr(est.sim.05,bootstrap.interval = TRUE)


## ----eval = FALSE-------------------------------------------------------------
#> Sys.sleep(100)

