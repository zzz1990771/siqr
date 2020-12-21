#' Data generation function for simulation and demonstration
#' There are three settings.
#'
#' @param n sample size
#' @param true.theta true single-index coefficients,
#' default is c(1,1,1)/sqrt(3) for setting 1 and c(1,2)/sqrt(5) for other settings
#' @param sigma the standard deviation of error term
#' @param setting chose from three settings
#' @param ncopy generates multiple copies of data for Monte Carlo simulations
#'
#' @return X  predictors
#' @return Y  response variables
#' @return single.index.values single index term
generate.data <- function(n,true.theta=NULL,sigma=0.1,setting="setting1",ncopy=1){

  if(setting == "setting1"){

    #parameter setting
    true.theta = if(is.null(true.theta)) c(1, 1, 1)/sqrt(3) else true.theta
    c1 = sqrt(3)/2-1.645/sqrt(12) #0.3912
    c2 = sqrt(3)/2+1.645/sqrt(12)#1.3409

    X = matrix(runif(length(true.theta)*n), ncol=length(true.theta))
    true.theta = sign(true.theta[1])*true.theta/sqrt(sum(true.theta^2));
    U = X%*%true.theta
    si = sin( (U-c1)*pi/(c2 -c1) )
    y = si + rnorm(length(si),0,sigma)
    if(ncopy>1){
      ylist <- lapply(vector(mode = "list", length = ncopy),function(x){si + rnorm(n,0,sigma)})
    }
  }else if(setting == "setting2"){

  }else if(setting == "setting3"){
    true.theta = if(is.null(true.theta)) c(1, 2)/sqrt(5) else true.theta
    X = matrix(rnorm(length(true.theta)*n), ncol=length(true.theta))
    U = X%*%true.theta
    si = 5*cos(U)+exp(-U^2)
    e = rexp(n,rate=.5)
    y = si+e;
    if(ncopy>1){
      ylist <- lapply(vector(mode = "list", length = ncopy),function(x){si + rexp(n,rate=.5)})
    }
  }

  if(ncopy>1){
    return(list("X" = X, "Y" = ylist,"single.index.values"=si))
  }else{
    return(list("X" = X, "Y" = y,"single.index.values"=si))
  }
}


#' A supporting function that return the local polynomial regression quantile.
#' This estimates the quantile and its derivative at the point x.0
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
#' @param tau left-tail probability (quantile index), scalar
#' @param beta.inital starting value of beta, the single index coefficients
#' @param maxiter max iteration number
#' @param tol toleration for convergence
#'
#' @return beta - the fitted single index coefficients with unit norm and first component being non negative
#'         flag.conv  - whether the iterations converge
#' @note in step 2 of the proposed algorithm, we may consider random sampling a "subsample" of size,
#'       say 5n, of the augmented sample(with sample size n^2);
#'        do it several times and take average of the estimates(need to write a sampling step, but not in this utility function
siqr<-function (y, X, tau=0.5, beta.inital=NULL, maxiter=40, tol=1e-9, method = "own")
{
  require(quantreg)
  if(is.null(beta.inital)){
    beta.inital <- coef(rq(y~X,tau=tau))[-1]
    beta.inital[1] <- abs(beta.inital[1])
  }
  flag.conv<-0; #flag whether maximum iteration is achieved

  beta.new<-beta.inital; #starting value
  if(method == "Wu"){
  }else{
    beta.new<-sign(beta.new[1])*beta.new/sqrt(sum(beta.new^2));
  }
  #beta.new<-sign(beta.new[1])*beta.new/sqrt(sum(beta.new^2));

  n<-NROW(y); d<-NCOL(X);
  a<-rep(0,n); b<-rep(0,n); #h<-rep(0,n);

  iter<-0;
  beta.old<-2*beta.new;


while((iter < maxiter) & (sum((beta.new-beta.old)^2)>tol))
#while(iter < maxiter)
 {
 #print(iter)
 #print(beta.new)

 beta.old<-beta.new;
 iter<-iter+1;
 ####################################
 #  step 1: compute a.j,b.j; j=1:n  #
 ####################################
  a<-rep(0,n); b<-rep(0,n);#h<-rep(0,n);
  x<-rep(0,n);
     for(jj in 1:d)
     {x<-x+X[,jj]*beta.old[jj]; #n-sequence, dim=null
       #x0<-x0+X[j,jj]*beta.old[jj]; #scalar
       }
   hm<-KernSmooth::dpill(x, y);
   hp<-hm*(tau*(1-tau)/(dnorm(qnorm(tau)))^2)^.2;
 x0<-0;
 for(j in 1:n)
  {
     x0<-x[j];
     fit<-lprq0(x, y, hp, tau, x0)
     a[j]<-fit$fv;
     b[j]<-fit$dv;
   }

 #############################
 # step 2: compute beta.new #
 #############################
 # here, let v.j=1/n;
 ynew<-rep(0,n^2);
 xnew<-rep(0,n^2*d);
 xnew<-matrix(xnew,ncol=d);

 for (i in 1:n)
  { for (j in 1:n)
    { ynew[(i-1)*n+j]<-y[i]-a[j];
      for(jj in 1:d){ xnew[(i-1)*n+j,jj]<-b[j]*(X[i,jj]-X[j,jj]);}
    }
  }

  xg<-rep(0,n^2); #x*beta
  for(jj in 1:d)
      {xg<-xg+xnew[,jj]*beta.old[jj]; #n-sequence, dim=null
       }
  xgh<-rep(0,n^2); #x*beta/h
  for (i in 1:n)
  {   for (j in 1:n)
      {
        xgh[(i-1)*n+j]<-xg[(i-1)*n+j]/hp;
      }
   }
  wts<-dnorm(xgh);
  #fit<-rq(ynew ~0+ xnew, weights = wts, tau = tau, method="fn") ; #pfn for very large problems
  fit<-rq(ynew ~0+ xnew, weights = wts, tau = tau, ci = FALSE) ; #default: br method, for several thousand obs
          # 0, to exclude intercept
  beta.new<-fit$coef;
  beta.new<-sign(beta.new[1])*beta.new/sqrt(sum(beta.new^2));   #normalize

} #end iterates over iter;


flag.conv<-(iter < maxiter)

beta<-beta.new;
names(beta) <- colnames(X)

si <- X%*%beta
hm <- KernSmooth::dpill(si,y);
hp <- hm*(tau*(1-tau)/(dnorm(qnorm(tau)))^2)^.2;

yhat<-rep(0,n);
for (i in 1:length(y)){
  local.fit<-lprq0(si, y, hp, tau, si[i]);
  yhat[i]<-local.fit$fv;
}

err<- y-yhat;
R<- sum(abs(err)+(2*tau-1)*err)/n;

list(beta=beta,flag.conv=flag.conv,X=X,y=y,yhat=yhat,tau=tau,rqfit=fit,MSAE = R)
}


#' plot function of siqr
#'
#'
#' @param model.obj The SIQR model object
#'
#' @return None
plot.siqr <- function(model.obj, bootstrap.interval = FALSE){
  si <- model.obj$X%*%model.obj$beta
  y <- model.obj$y
  plot(si,y,xlab = "Single Index", ylab = "Predicted Y",col="gray",main="Fitted Quantiles Plot");
  lines(sort(si),model.obj$yhat[order(si)],lty=1,lwd=1.5,col="red");

  if(bootstrap.interval){
    tau <- model.obj$tau
    hm <- KernSmooth::dpill(si,y)
    hp <- hm*(tau*(1-tau)/(dnorm(qnorm(tau)))^2)^.2

    #get residual
    res <- y-model.obj$yhat
    n <- length(res)

    #get bootstrap y.hat
    #v1
    # B=100
    # y.hat.B <- matrix(NA,length(y.B),B)
    # for(b in 1:B){
    #   #get residual bootstrap data
    #   bs.index<-sample(n,replace=T)
    #   res.B<-res[bs.index]
    #   y.B<-model.obj$yhat+res.B
    #   fit.B <- siqr(y.B, X, beta.inital = beta0, tau=tau,maxiter = 20,tol = 1e-6, method = "Wu")
    #   y.hat.B[,b] <- fit.B$yhat
    # }


    #v2
    B=100
    y.hat.B <- matrix(NA,length(y),B)
    for(b in 1:B){
      for(i in 1:length(y)){
        #get residual bootstrap data
        bs.index<-sample(n,replace=T)
        res.B<-res[bs.index]
        y.B<-model.obj$yhat+res.B
        fit.B <- lprq0(si, y.B, hp, tau=tau, si[i])
        y.hat.B[i,b] <- fit.B$fv
      }
    }

    #get sd of bootstrap Y.hat
    se.yhat <- apply(y.hat.B,1,sd)
    #2*sd +/- original y.hat to form the interval
    yhat.B.025 <- model.obj$yhat - 2 * se.yhat
    yhat.B.975 <- model.obj$yhat + 2 * se.yhat
    #plot
    #plot.si(model.obj = model.obj)
    lines(sort(si),yhat.B.025[order(si)],lty=6,lwd=1.5,col="blue")
    lines(sort(si),yhat.B.975[order(si)],lty=6,lwd=1.5,col="blue")
  }
}


#' Function to print summary
#'
#' @param siqr.object the single index quantile regression model object
#'
#' @return the summarized information object
print.summary.siqr <- summary.siqr <- function(object, digits = max(5, getOption("digits") - 3),
                              signif.stars = getOption("show.signif.stars"), ...)

{
  cat("Tau: ",object$tau,"\n",sep="")

  if (length(object$beta)>0)
  { cat("\nsingle index coefficients:\n")
    printCoefmat(data.frame(Coefficients=est$beta), digits = digits, signif.stars = signif.stars, na.print = "NA", ...)
  }
  cat("\n")
  cat("Model MSAE: ",object$MSAE,"\n",sep="")
  cat("Model convergence: ",object$flag.conv,"\n",sep="")
  invisible(object)
}

