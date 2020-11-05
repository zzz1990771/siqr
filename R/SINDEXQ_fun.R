#main program for simulations
#contains: library, utility, data, estimation
#July 5, 2006; by Tracy Wu

##1.  invoke libraries

################################################################################
#library(stats);
#library(quantreg); #will call function lprq(), rq()
#library(lokern); # call function glkerns() to compute h_mean
#library(mda) #call for polyreg(), this is not the right polyreg; use own
#library(KernSmooth);
#################################################################################

##2.  utility functions


## Data generation function for simulation and demonstration
## A sine-bump setting has been employed.
generate_data <- function(n,true.theta=c(1, 1, 1)/sqrt(3),sigma=0.1,family="Gaussian",ncopy=1){
  #parameter setting
  c1 = sqrt(3)/2-1.645/sqrt(12) #0.3912
  c2 = sqrt(3)/2+1.645/sqrt(12)#1.3409
  #rho = 0.3

  X = matrix(runif(length(true.theta)*n), ncol=length(true.theta))
  true.theta = sign(true.theta[1])*true.theta/sqrt(sum(true.theta^2));
  U = X%*%true.theta
  si = sin( (U-c1)*pi/(c2 -c1) )

  if(family=="Gaussian"){
    y = si + rnorm(length(si),0,sigma)
  }else if(family=="Binomial"){
    py = exp(si)/(1+exp(si))
    y = rbinom(length(si), size=1, prob=py)
  }else if(family=="Poisson"){
    py = exp(si)
    y = rpois(length(si),py)
  }
  if(ncopy>1){
    ylist <- lapply(vector(mode = "list", length = ncopy),function(x){si + rnorm(length(si),0,sigma)})
    return(list("X" = X, "Y" = ylist,"single_index_values"=si))
  }else{
    return(list("X" = X, "Y" = y,"single_index_values"=si))
  }
}


#' A supporting function that return the local polynomial regression quantile.
#' This estimates the quantile and its derivative at the point x_0
#'
#' @param x covariate sequence; y - response sequence; they form a sample.
#' @param h bandwidth(scalar); tau - left-tail probability
#' @param x0 point at which the quantile is estimated
#'
#' @return x0  a scalar
#' @return fv  quantile est; dv - quantile derivative est


lprq0<-function (x, y, h, tau = 0.5,x0)  #used in step 1 of the algorithm
{
    require(quantreg)
    fv <- x0
    dv <- x0

    z <- x - x0
    wx <- dnorm(z/h)
    r <- rq(y ~ z, weights = wx, tau = tau, ci = FALSE)
    fv <- r$coef[1]
    dv <- r$coef[2]
    list(x0 = x0, fv = fv, dv = dv)
}


#' Main estimation function of single index quantile regression model.
#' a two step method.
#'
#'
#' @param y response vector;
#' @param X covariate matrix;
#' @param p left-tail probability (quantile index), scalar
#' @param gamma.inital starting value of gamma, the single index cooefficients
#' @param maxiter max iteration number
#' @param tol toleration for convergence
#'
#' @return gamma - index direction, unit norm, first nonneg component
#'          flag.conv  - whether the iterations converge
#' @note in step 2 of the proposed algorithm, we may consider random sampling a "subsample" of size,
#'       say 5n, of the augmented sample(with sample size n^2);
#'        do it several times and take average of the estimates(need to write a sampling step, but not in this utility function
siqr<-function (y, X, p=0.5, gamma.inital=NULL, maxiter=40, tol=1e-9, method = "own")
{
  require(quantreg)
  if(is.null(gamma.inital)){
    gamma.inital <- coef(rq(y~X,tau=p))[-1]
  }
  flag.conv<-0; #flag whether maximum iteration is achieved

  gamma.new<-gamma.inital; #starting value
  if(method == "Wu"){
  }else{
    gamma.new<-sign(gamma.new[1])*gamma.new/sqrt(sum(gamma.new^2));
  }
  #gamma.new<-sign(gamma.new[1])*gamma.new/sqrt(sum(gamma.new^2));

  n<-NROW(y); d<-NCOL(X);
  a<-rep(0,n); b<-rep(0,n); #h<-rep(0,n);

  iter<-1;
  gamma.old<-2*gamma.new;


while((iter < maxiter) & (sum((gamma.new-gamma.old)^2)>tol))
#while(iter < maxiter)
 {
 #print(iter)
 gamma.old<-gamma.new;
 iter<-iter+1;
 ####################################
 #  step 1: compute a_j,b_j; j=1:n  #
 ####################################
  a<-rep(0,n); b<-rep(0,n);#h<-rep(0,n);
  x<-rep(0,n);
     for(jj in 1:d)
     {x<-x+X[,jj]*gamma.old[jj]; #n-sequence, dim=null
       #x0<-x0+X[j,jj]*gamma.old[jj]; #scalar
       }
   hm<-KernSmooth::dpill(x, y);
   hp<-hm*(p*(1-p)/(dnorm(qnorm(p)))^2)^.2;
 x0<-0;
 for(j in 1:n)
  {
     x0<-x[j];
     fit<-lprq0(x, y, hp, p, x0)
     a[j]<-fit$fv;
     b[j]<-fit$dv;
   }

 #############################
 # step 2: compute gamma.new #
 #############################
 # here, let v_j=1/n;
 ynew<-rep(0,n^2);
 xnew<-rep(0,n^2*d);
 xnew<-matrix(xnew,ncol=d);

 for (i in 1:n)
  { for (j in 1:n)
    { ynew[(i-1)*n+j]<-y[i]-a[j];
      for(jj in 1:d){ xnew[(i-1)*n+j,jj]<-b[j]*(X[i,jj]-X[j,jj]);}
    }
  }

  xg<-rep(0,n^2); #x*gamma
  for(jj in 1:d)
      {xg<-xg+xnew[,jj]*gamma.old[jj]; #n-sequence, dim=null
       }
  xgh<-rep(0,n^2); #x*gamma/h
  for (i in 1:n)
  {   for (j in 1:n)
      {
        xgh[(i-1)*n+j]<-xg[(i-1)*n+j]/hp;
      }
   }
  wts<-dnorm(xgh);
  #fit<-rq(ynew ~0+ xnew, weights = wts, tau = p, method="fn") ; #pfn for very large problems
  fit<-rq(ynew ~0+ xnew, weights = wts, tau = p, ci = FALSE) ; #default: br method, for several thousand obs
          # 0, to exclude intercept
  gamma.new<-fit$coef;
  gamma.new<-sign(gamma.new[1])*gamma.new/sqrt(sum(gamma.new^2));   #normalize

} #end iterates over iter;

iter<-iter;

flag.conv<-(iter < maxiter) ;# = 1 if converge; =0 if not converge
#flag.conv<- 1- ((iter=maxiter)&(sum((gamma.new-gamma.old)^2)<tol))

gamma<-gamma.new;
names(gamma) <- colnames(X)

si <- X%*%gamma
hm <- KernSmooth::dpill(si,y);
hp <- hm*(p*(1-p)/(dnorm(qnorm(p)))^2)^.2;

yhat<-rep(0,n);
for (i in 1:length(y)){
  local_fit<-lprq0(si, y, hp, p, si[i]);
  yhat[i]<-local_fit$fv;
}

err<- y-yhat;
R<- sum(abs(err)+(2*p-1)*err)/n;

list(gamma=gamma,flag.conv=flag.conv,X=X,y=y,yhat=yhat,p=p,rqfit=fit,MSAE = R)
}


#' plot function of siqr
#'
#'
#' @param model.obj The SIQR model object
#'
#' @return None
plot.si <- function(model.obj, bootstrap_interval = FALSE){
  si <- model.obj$X%*%model.obj$gamma
  y <- model.obj$y
  plot(si,y,xlab = "Single Index", ylab = "Predicted Y",col="gray");
  lines(sort(si),model.obj$yhat[order(si)],lty=1,lwd=1.5,col="red");

  if(bootstrap_interval){
    p <- model.obj$p
    hm <- KernSmooth::dpill(si,y)
    hp <- hm*(p*(1-p)/(dnorm(qnorm(p)))^2)^.2

    #get residual
    res <- y-model.obj$yhat
    n <- length(res)

    #get bootstrap y_hat
    #v1
    # B=100
    # y_hat_B <- matrix(NA,length(y.B),B)
    # for(b in 1:B){
    #   #get residual bootstrap data
    #   bs.index<-sample(n,replace=T)
    #   res.B<-res[bs.index]
    #   y.B<-model.obj$yhat+res.B
    #   fit.B <- siqr(y.B, X, gamma.inital = gamma0, p=p,maxiter = 20,tol = 1e-6, method = "Wu")
    #   y_hat_B[,b] <- fit.B$yhat
    # }


    #v2
    B=100
    y_hat_B <- matrix(NA,length(y),B)
    for(b in 1:B){
      for(i in 1:length(y)){
        #get residual bootstrap data
        bs.index<-sample(n,replace=T)
        res.B<-res[bs.index]
        y.B<-model.obj$yhat+res.B
        fit.B <- lprq0(si, y.B, hp, tau=p, si[i])
        y_hat_B[i,b] <- fit.B$fv
      }
    }

    #get sd of bootstrap Y_hat
    se_yhat <- apply(y_hat_B,1,sd)
    #2*sd +/- original y_hat to form the interval
    yhat_B_025 <- model.obj$yhat - 2 * se_yhat
    yhat_B_975 <- model.obj$yhat + 2 * se_yhat
    #plot
    #plot.si(model.obj = model.obj)
    lines(sort(si),yhat_B_025[order(si)],lty=6,lwd=1.5,col="blue")
    lines(sort(si),yhat_B_975[order(si)],lty=6,lwd=1.5,col="blue")
  }
}


























#' Main estimation function of single index quantile regression model.
#' a two step method.
#' revised version of index.gamma(), in that random sampling is implemented in its step 2
#'
#'
#' @param y - response vector; X - covariate matrix;
#' @param p - left-tail probability (quantile index), scalar
#' @param gamma.inital - starting value of gamma
#' @param maxiter
#' @param tol
#' @param subsize
#' @param subrep
#'
#' @return gamma - index direction, unit norm, first nonneg component
#'          flag.conv  - whether the iterations converge
#' @note TBD
index.gamma1<-function (y, X, p, gamma.inital,maxiter,tol,subsize,subrep)
{
flag.conv<-0; #flag whether maximum iteration is achieved

gamma.new<-gamma.inital; #starting value
gamma.new<-sign(gamma.new[1])*gamma.new/sqrt(sum(gamma.new^2));

n<-NROW(y); d<-NCOL(X);
a<-rep(0,n); b<-rep(0,n); #h<-rep(0,n);
gamma.inter<-rep(0,d*subrep);

iter<-1;
gamma.old<-2*gamma.new;
#gamma.old<-sign(gamma.old[1])*gamma.old/sqrt(sum(gamma.old^2));

while((iter < maxiter) & (sum((gamma.new-gamma.old)^2)>tol))
#while(iter < maxiter)
 {
 gamma.old<-gamma.new;
 iter<-iter+1;
 ####################################
 #  step 1: compute a_j,b_j; j=1:n  #
 ####################################
  a<-rep(0,n); b<-rep(0,n);#h<-rep(0,n);
  x<-rep(0,n);
     for(jj in 1:d)
     {x<-x+X[,jj]*gamma.old[jj]; #n-sequence, dim=null
       #x0<-x0+X[j,jj]*gamma.old[jj]; #scalar
       }
   hm<-KernSmooth::dpill(x, y);
   hp<-hm*(p*(1-p)/(dnorm(qnorm(p)))^2)^.2;
 x0<-0;
 for(j in 1:n)
  {
     x0<-x[j];
     fit<-lprq0(x, y, hp, p, x0)
     a[j]<-fit$fv;
     b[j]<-fit$dv;
   }

 #############################
 # step 2: compute gamma.new #
 #############################
 # here, let v_j=1/n;
 # augmented "observations"
 ynew<-rep(0,n^2);
 xnew<-rep(0,n^2*d);
 xnew<-matrix(xnew,ncol=d);
 for (i in 1:n)
    { for (j in 1:n)
      { ynew[(i-1)*n+j]<-y[i]-a[j];
      for(jj in 1:d){ xnew[(i-1)*n+j,jj]<-b[j]*(X[i,jj]-X[j,jj]); }
     }
   }  #generated ynew and xnew
  xg<-rep(0,n^2); #x*gamma
  for(jj in 1:d)
      {xg<-xg+xnew[,jj]*gamma.old[jj]; #n-sequence, dim=null
       }
  xgh<-rep(0,n^2); #x*gamma/h
  for (i in 1:n)
  {   for (j in 1:n)
      {
        xgh[(i-1)*n+j]<-xg[(i-1)*n+j]/hp;
      }
   }
 wts<-dnorm(xgh); #generated wts, the weights
   ######################
   #  sampling step
   ######################
 gamma.inter<-rep(0,d*subrep);
 gamma.inter<-matrix(gamma.inter, ncol=subrep);  #to store intermediate gammas
 for (subiter in 1:subrep)
  {   id.sample<-  sample(seq(1,n^2), subsize) #default: w/o replacement
      subx<-xnew[id.sample,]
      suby<-ynew[id.sample]
      subw<-wts[id.sample]

      subfit<-rq(suby ~0+ subx, weights = subw, tau = p, ci = FALSE) ; #default: br method, for several thousand obs
          # 0, to exclude intercept
    subgamma<-subfit$coef;
    subgamma<-sign(subgamma)*subgamma/sqrt(sum(subgamma^2));   #normalize
    gamma.inter[,subiter]<-subgamma;
  }  #end iterates over subiter
 ###################################
 #take average of gamma estimates
 ###################################
 gamma.inter<-gamma.inter #dim: d by subrep
 gamma.new<-apply(gamma.inter,1,mean)

 gamma.new<- sign(gamma.new)*gamma.new/sqrt(sum(gamma.new^2))

} #end iterates over iter;

gamma.new<-gamma.new
iter<-iter;

flag.conv<-(iter < maxiter) ;# = 1 if converge; =0 if not converge
#flag.conv<- 1- ((iter=maxiter)&(sum((gamma.new-gamma.old)^2)<tol))

gamma<-gamma.new;
list(gamma=gamma,flag.conv=flag.conv)
}

####################
#end of index.gamma1
#####################







