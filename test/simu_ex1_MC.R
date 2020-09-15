############################################################################################
# Simulation example 1: sine-bump. 
# Monte Carlo variance study. 
# updated 12/20/2009, by Yan Yu
# Compare asymptotic standard error estmate, bootstrap standard error estmate, with 
# the Monte Carlo standard error estmate (approximation to the truth) 
# Step 1: Source R functions from SimulationStudy.R
# Step 2: Run simu_x.R to fix xx design matrix for parrallel computing
# Step 3: simu_ex1_MC.R
##############################################################################################
# The bootstrap s.e. of gammahat is calculated by sampling the data with replacement B times;
# asymptotic s.e. and average bootstrap s.e. are compared with the Monte Carlo s.e. 
##############################################################################################

library(MASS); #ginv function;
## fix xx design matrix for parrallel computing
load("SimuXX.Rdata");

n<-200
#x1<-runif(n, min=0, max=1); x2<-runif(n, min=0, max=1); x3<-runif(n, min=0, max=1);
#xx<-cbind(x1,x2,x3);
#save(xx,file="SimuXX.Rdata");

nsimu<-100;
B<-5; #number of bootstrap samples;
p<-0.9;  
n<-200
sigma<-.1;
maxiter<-30; #maximum number of iteration 
crit<-1e-8;  

gamma.true<-c(1,1,1);
gamma.true<-sign(gamma.true[1])*gamma.true/sqrt(sum(gamma.true^2));
d<-length(gamma.true);
gamma0<-c(1,2,3);  #starting value; do not confuse with true value 

A<-sqrt(3)/2-1.645/sqrt(12); C<-sqrt(3)/2+1.645/sqrt(12);
index<-xx%*%gamma.true;
mu<-sin(pi*(index-A)/(C-A));

##Asy. s.e.
#linear regression to find E(x|x\trans\gamma);use unconditional to replace conditional expectation;
mean.xx<-apply(xx,2,mean);
xc<-xx-mean.xx;
partial.g0<-cos(pi*(index-A)/(C-A))*pi/(C-A); #n by 1 vector;
zp<-qnorm(p);
pdf.y<-exp(-(0.1*zp)^2/(2*sigma^2))/(sigma*sqrt(pi*2));
xc.1<-matrix(nrow=n,ncol=d);
xc.2<-matrix(nrow=n,ncol=d);
i=1
for (i in 1:n)
{ xc.1[i,]<-partial.g0[i]*xc[i,];
  xc.2[i,]<-sqrt(pdf.y)*partial.g0[i]*xc[i,];
}
c0<-t(xc.1)%*%(xc.1)/n;    #matrix product;
c1<-t(xc.2)%*%(xc.2)/n;   
var<-p*(1-p)*ginv(c1)%*%c0%*%ginv(c1)/n;
asyse<-sqrt(diag(var));

#initialization
gamma.est<- matrix(rep(0,d*nsimu),nrow=d);
conv<-rep(0,nsimu);
bootse<-matrix(nrow=d,ncol=nsimu);

t.start<-date()
simu=1;
for (simu in 1:nsimu)
{
e<-rnorm(n,mean=0,sd=sigma);
y<-mu+e; 

est<-index.gamma(y,xx,p,gamma0,maxiter,crit);
gamma.est[,simu]<-est$gamma;
conv[simu]<-est$flag.conv;

## Bootstrap s.e.
gammaB.matrix<-matrix(nrow=d,ncol=B);
#gammaB.matrix<-array(rep(0,d*B*nsimu),c(d,B,nsimu))

##Gooijer Bootstrap se
gammaB.matrix<-matrix(nrow=d,ncol=B);

for (bs in 1:B) {
indexhat<-xx%*%est$gamma;
muhat<-sin(pi*(indexhat-A)/(C-A));
r <- y - muhat;
bs.index<-sample(n,replace=T);
rB <- r[bs.index];
yB<- muhat + rB;
xxB<-xx;
gammaB<-index.gamma(yB, xxB, p, gamma0,maxiter,crit);
gammaB.matrix[,bs]<-gammaB$gamma;
gammaB.std<-apply(gammaB.matrix,1,sd);
}
bootse[,simu]<-gammaB.std;

} #end MC loop
abootse <- apply(bootse,1,mean)


##Monte Carlo se
gamma.avg<-apply(gamma.est,1,mean)
mcse<-apply(gamma.est,1,sd)
##Relative distance measure of asy. se. and boot se. comparing with MC se.
asyd2 <- (sqrt(sum((asyse-mcse)^2)))/(sqrt(sum(mcse^2))) 
asyd1 <- (sum(abs(asyse-mcse)))/(sum(abs(mcse)))
bootd2 <- (sqrt(sum((abootse-mcse)^2)))/(sqrt(sum(mcse^2))) 
bootd1 <- (sum(abs(abootse-mcse))) /(sum(abs(mcse)))

save(p,nsimu,B,xx,conv,gamma.true,gamma.est,gamma.avg,asyse,bootse,mcse,gammaB.matrix,file="Simu.Rdata");





