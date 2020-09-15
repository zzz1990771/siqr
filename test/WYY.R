

###########################################################
#
#  simulations
#
###########################################################

#################################################################
# simulation example 1: sine-bump. generate data 1 (homoscedastic, normal)
#################################################################

#true para: (1,1,1)/sqrt(3)

t.start<-date()
#nsimu<-50; #number of replications;
#usually each 50 simulations generates a result table to be combined later
nsimu<-50;
n<-200
subsize=20*n;
subrep=2;
sigma=0.1;
maxiter<-30; #25
crit<-1e-8;  #9
p<-0.7;  # p=0.5 is best for estimation of gamma

gamma.est<- matrix(rep(0,3*nsimu),nrow=3);
conv<-rep(0,nsimu);

for (simu in 1:nsimu)
{
  data <- generate_data(n)
  xx   <- data$X
  y    <- data$Y

 #initial value
 gamma0<-c(1,2,0);  #starting value; do not confuse with true value
 est<-siqr(y,xx,p,gamma0,maxiter,crit);
 gamma.est[,simu]<-est$gamma;
 conv[simu]<-est$flag.conv;
}
gamma.est<-gamma.est;
gamma.est
conv<-conv
t.end<-date()

t.start
t.end

##monte carlo study
gamma.avg<-apply(gamma.est,1,mean)
gamma.std<-apply(gamma.est,1,sd)


################################################################
# simulation example 2: heteroscedastic, normal
# true model: (#Y=sin(0.75eta)+1+{ 0.3/sqrt{sin(0.75eta)+1} }Z )
#            y=10*sin(0.75*index)+sqrt(sin(index)+1)*z;
# true para:  (1,2)/sqrt(5)
# updated: July 19, 2006
################################################################
t.start<-date()
nsimu<-3;    #usually each 50 simulations generates a result table to be combined later
gamma.est<- rep(0,2*nsimu);
gamma.est<- matrix(gamma.est,nrow=2);
conv<-rep(0,nsimu);
n<-400; #sample size is large!!
subsize<-20*n;
subrep<-2;
#gamma.true<-c(1,2); #true value
gamma.true<-c(1,-2); #true value
gamma.true<-sign(gamma.true[1])*gamma.true/sqrt(sum(gamma.true^2));
gamma0<-c(1,-1)      #starting value; do not confuse with true value
gamma0<-sign(gamma0[1])*gamma0/sqrt(sum(gamma0^2));
gamma.true
gamma0
p<-0.5;
maxiter<-40; #25
crit<-1e-9;  #9

for (simu in 1:nsimu)
{
#generate data
z<-rnorm(n*1.2);
x1<-rnorm(n*1.2,mean=0,sd=0.25);
x2<-rnorm(n*1.2,mean=0,sd=0.25);
xx<-cbind(x1,x2);
index<-xx%*%gamma.true;
#y<-sin(0.75*index)+1+0.1*sqrt(sin(0.75*index)+1)*z; ###true model
y<-10*sin(0.75*index)+sqrt(sin(index)+1)*z;
# may check scatter plot of (y,index)
data<-cbind(index,xx,y)
ix<-sort(index,index=T)$ix
id<-ix[ceiling(0.1*n):(ceiling(0.1*n)+n-1)]
subdata<-data[id,]  #only use interior data
index<-subdata[,1];
xx<-subdata[,2:3]
y<-subdata[,4]

#start estimation
#est<-siqr(y,xx,p,gamma0,maxiter,crit);
est<-index.gamma1(y,xx,p,gamma0,maxiter,crit,subsize,subrep);
gamma.est[,simu]<-est$gamma;
conv[simu]<-est$flag.conv;
}  #end iterates over simu

gamma.est<-gamma.est;
gamma.est
conv<-conv
t.end<-date()

gamma.avg<-apply(gamma.est,1,mean)
gamma.std<-apply(gamma.est,1,sd)


##############################
#  plot the histograms
##############################
#id1<-sort(gamma.total[1,],index=T)$ix
#id<-id1[8:107]
#gamma.total2<-gamma.total[,id] #exclude outliers

gamma.total<-gamma.total0[,c(-2,-4,-5,-6,-8,-10,-19,-21,-22)]
m<-apply(gamma.total,1,mean)
s<-apply(gamma.total,1,sd)
apply(gamma.total,1,median)

par(mfcol=c(2,1))
mfg=c(1,1,2,1)
#draw histogram of gamma and overlay histogram with normal distribution with estimated mean, sd
grid1<-seq(-0.38,0.58,0.001)
hist(gamma.total[1,],freq=FALSE,col="lightblue",ylim=c(0,15), xlab="gamma1",main="Histograms of estimates: heteroscedastic");
# true=0.4472136
lines(grid1,dnorm(grid1, mean=m[1], sd= s[1]))
lines(rep(gamma.true[1],400),seq(0,15,length.out=400),lwd=4 )

mfg=c(2,1,2,1)
grid2<-seq(0.82,0.94,0.001)
#hist(gamma.total[2,])
#hist(gamma.total2[2,],br=c(0.855,0.865,0.875,0.885,0.895,0.905,0.915,0.925),freq=TRUE,main="boxplot for $\gamma_2$")
hist(gamma.total[2,],freq=FALSE,col="lightblue",ylim=c(0,32),xlab="gamma2",main=" ")
# true=0.8944272
lines(grid2,dnorm(grid2, mean=m[2], sd=s[2]));
lines(rep(gamma.true[2],400),seq(0,32,length.out=400),lwd=4 )

dev.off()
########################################################


##################################################
# simulation example 3; Exponentional error
#
# July 19, 2006
##################################################

t.start<-date()
nsimu<-60;    #usually each 50 simulations generates a result table to be combined later
gamma.est<- rep(0,2*nsimu);
gamma.est<- matrix(gamma.est,nrow=2);
conv<-rep(0,nsimu);
n<-400;
subsize=30*n;
subrep=1;
gamma.true<-c(1,2); #true value
gamma.true<-sign(gamma.true[1])*gamma.true/sqrt(sum(gamma.true^2));
gamma0<-c(1,1)      #starting value; do not confuse with true value
gamma0<-sign(gamma0[1])*gamma0/sqrt(sum(gamma0^2));
gamma.true
gamma0
p<-0.5;
maxiter<-40; #25
crit<-1e-9;  #9


for (simu in 1:nsimu)
{
#generate data
n<-400;
n1<-n*1.2
x1<-rnorm(n1);
x2<-rnorm(n1);
xx<-cbind(x1,x2);
e<-rexp(n1,rate=.5) ; #exp(1)
gamma.true<-c(1,2);
gamma.true<-sign(gamma.true[1])*gamma.true/sqrt(sum(gamma.true^2));
index<-xx%*%gamma.true;
y<- 5*cos(index)+exp(-index^2)+e;

data<-cbind(index,xx,y)
ix<-sort(index,index=T)$ix
id<-ix[ceiling(0.1*n):(ceiling(0.1*n)+n-1)]
subdata<-data[id,]
##############only use interior data to avoid sparse data at boundary
index<-subdata[,1];
xx<-subdata[,2:3]
y<-subdata[,4]
#  plot(index,y, main="$y$ vs index");

#start estimation
#est<-siqr(y,xx,p,gamma0,maxiter,crit);
est<-index.gamma1(y,xx,p,gamma0,maxiter,crit,subsize,subrep);
gamma.est[,simu]<-est$gamma;
conv[simu]<-est$flag.conv;
}  #end iterates over simu

gamma.est<-gamma.est;
gamma.est
conv<-conv
t.end<-date()
#################################
# end of simulations
#################################

gamma.avg<-apply(gamma.est,1,mean)
gamma.std<-apply(gamma.est,1,sd)

gamma.true<-c(1,2)/sqrt(5)
gamma.total0<-cbind(gamma.estA,gamma.estB) # 110 results
gamma.total<-gamma.total0[,c(-32,-33,-75,-80,-90,-108,-109)]
m<-apply(gamma.total,1,mean)
s<-apply(gamma.total,1,sd)
apply(gamma.total,1,median)

par(mfcol=c(2,1))
mfg=c(1,1,2,1)
#draw histogram of gamma and overlay histogram with normal distribution with estimated mean, sd
grid1<-seq(-0.38,0.58,0.001)
hist(gamma.total[1,],freq=FALSE,col="lightblue",ylim=c(0,15), xlab="gamma1",main="Histograms of estimates: exponential");
# true=0.4472136
lines(grid1,dnorm(grid1, mean=m[1], sd=s[1]))
lines(rep(gamma.true[1],400),seq(0,15,length.out=400),lwd=4 )

mfg=c(2,1,2,1)
grid2<-seq(0.82,0.94,0.001)
#hist(gamma.total[2,])
#hist(gamma.total2[2,],br=c(0.855,0.865,0.875,0.885,0.895,0.905,0.915,0.925),freq=TRUE,main="boxplot for $\gamma_2$")
hist(gamma.total[2,],freq=FALSE,col="lightblue",ylim=c(0,32),xlab="gamma2",main=" ")
# true=0.8944272
lines(grid2,dnorm(grid2, mean=m[2], sd= s[2]));
lines(rep(gamma.true[2],400),seq(0,32,length.out=400),lwd=4 )

dev.off()


########################################################
# end of simulation example 3, last edit: July 19,2006
########################################################


#############################################################################
#   boston housing example
#    July 20, 2006
#############################################################################

library(MASS)
#help(Boston)
y<-Boston[,14]
rm<-Boston[,6]
lstat<-Boston[,13]
dis<-Boston[,8]

########################################
##  preliminary
#standize the covariates
y0<-(y-mean(y))/sd(y)
rm0<- (rm-mean(rm))/sd(rm)
lstat0<-(lstat-mean(lstat))/sd(lstat)
loglstat<-log(lstat)
dis0<-(dis-mean(dis))/sd(dis)
logdis<-log(dis)
#mean(rm0)   sd(rm0)
###########################################
# plot before transformation
#
par(mfcol=c(1,3))
mfg=c(1,1,1,3)
plot(rm0,y0)
mfg=c(1,2,1,3)
plot(lstat0,y0)
mfg=c(1,3,1,3)
plot(dis0,y0)

# or
#par(mfcol=c(1,3));mfg=c(1,1,1,3);plot(y0,rm0); mfg=c(1,2,1,3);plot(y0,lstat0); mfg=c(1,3,1,3); plot(y0,dis0)}

############################################
# plot after transformation
#
par(mfcol=c(1,3))
mfg=c(1,1,1,3)
plot(rm0,y0)
mfg=c(1,2,1,3)
plot(lstat0,y0)
mfg=c(1,3,1,3)
plot(logdis,y0)
#or
#par(mfcol=c(1,3));mfg=c(1,1,1,3);plot(y0,rm0);mfg=c(1,2,1,3);plot(y0,lstat0);mfg=c(1,3,1,3);plot(y0,logdis)
########################################3

yfit<-lm(y0~rm0+lstat0+dis0)

#gamma0<-c(1,-2,.5)  #starting value; can do usual regression to determine the sign first

coef<-yfit$coef
gamma0<-coef[2:4]
gamma0<-sign(gamma0[1])*gamma0/sqrt(sum(gamma0^2));
gamma0

#xx<-cbind(rm0,lstat0,logdis)
n<-length(y)
subsize=30*n;
subrep=3;

###########################################
p<-0.5;
maxiter<-40; #25
crit<-1e-9;  #9
#est<-index.gamma1(y0,xx,p,gamma0,maxiter,crit,subsize,subrep);
est<-siqr(y0,xx,p,gamma0,maxiter,crit);
gamma.est<-est$gamma;
conv<-est$flag.conv;

##########################################
# try p=0.1, 0.25, 0.5, 0.75, 0.90
##########################################
p<-c(0.1,0.25,0.5,0.75,0.9)
np<-length(p)
gamma.est<- rep(0,3*np);
gamma.est<- matrix(gamma.est,nrow=3);
conv<-rep(0,np)
for (i in 1:np)
{
est<-siqr(y,xx,p[i],gamma0,maxiter,crit);
gamma.est[,i]<-est$gamma;
conv[i]<-est$flag.conv
}

##########use interior index points
xx<-cbind(rm0,lstat0,dis0)
index.est<-xx%*%gamma.est;
data<-cbind(index.est,xx,y)
ix<-sort(index,index=T)$ix
id<-ix[ceiling(0.02*n):n-ceiling(0.02*n)]
subdata<-data[id,]

index1<-subdata[,1];
xx1<-subdata[,2:4]
y1<-subdata[,5]
####plot y vs index
index2<-rep(0,length(index1))
index2<-index1

#hm<-dpill(index1,y1);
#hp<-hm*(p*(1-p)/(dnorm(qnorm(p)))^2)^.2;
newfit<-lprq(index1, y1, hp, tau = .5)

plot(index1,y1)
plot(index,y)

lines(newfit$xx,newfit$fv)







