rm(list=ls())
set.seed(2023-6-19)
alp<-c(0.5,0.4,0.1);mu<-c(-1,12,30);sig<-c(2,3,4);beta<-1
sig<-sig*beta
#try what can be found?
y<-c();n<-1000;#sample size
m1<-30#maximum steps
m<-3#three components
for(i in 1:n)
{u<-runif(1)
if(u<alp[1]) y[i]<-rnorm(1,mu[1],sig[1])  
else if (u<alp[1]+alp[2]) y[i]<-rnorm(1,mu[2],sig[2])   
else y[i]<-rnorm(1,mu[3],sig[3]) }   #geneatiing data
par(mfrow=c(1,1))
hist(y)
alp0<-matrix(0,m1,m);p<-matrix(0,n,m)
mu0<-matrix(0,m1,m);sig0<-matrix(0,m1,m)
x<-quantile(y,c(0.25,0.75,0.9));
alp0[1,]<-rep(1,m)/m; mu0[1,]<-x;sig0[1,]<-c(1,1,1)# starting values
eps<-0.0001;diff<-1;
t<-1
while(diff>eps)
  {
  for(i in 1:n)
    {
    tem<-0
    for(k in 1:m)
      tem<-tem+alp0[t,k]*dnorm(y[i],mu0[t,k],sqrt(sig0[t,k]))
    for(j in 1:m)
      p[i,j]<-alp0[t,j]*dnorm(y[i],mu0[t,j],sqrt(sig0[t,j]))/tem 
    }
  #conditional expectation of missing data
  for(j in 1:m)
    {
    tem1<-0;tem2<-0;tem3<-0;
    for(i in 1:n)
      {
      tem1<-tem1+p[i,j];    tem2<-tem2+p[i,j]*y[i]
      tem3<-tem3+p[i,j]*(y[i]*y[i])
      }
    mu0[t+1,j]<-tem2/tem1;
    alp0[t+1,j]<-tem1/n;
    sig0[t+1,j]<-tem3/tem1-mu0[t+1,j]^2
    }
  diff<-sum(mu0[t+1,]-mu0[t,])^2+sum(sig0[t+1,]-sig0[t,])^2
  diff<-diff+sum(alp0[t+1,]-alp0[t,])^2
  t<-t+1
}
par(mfrow=c(1,1))
matplot(mu0[1:t,],type='l')
abline(h=mu,lty=c(4:5))
matplot(sig0[1:t,],type='l')
abline(h=sig^2,lty=c(4:5))
matplot(alp0[1:t,],typ='l',ylim=c(0,0.8))
abline(h=alp,lty=c(4:5))
