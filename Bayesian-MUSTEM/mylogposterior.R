
mylogposterior<-function(beta,data){

#mu<-c(1.2,.12)
#sigma<-matrix(c(.01,.04,.04,0.2),2,2)
mu<-c(1.2,.12,.04,.02)
sigma<-matrix(c(.01,.04,.003,0.2,.04,.07,0.2,0.05,0.003,0.2,0.3,0.2,0.2,0.05,0.2,0.6),4,4)
#mu<-rep(2,4)
#cmu<-matrix(mu,ncol=4)
#nmu<-rnorm(4,25,4)
#V<-diag(1/nmu)
#sigma<-cmu %*% V %*% t(cmu)
#sigma<-matrix(c(.01,.04,.003,0.2,.04,.07,0.2,0.05,0.003,0.2,0.3,0.2,0.2,0.05,0.2,0.6),4,4)

#This is initial value of covariance matrix. We have to decide this values carefully. Set up this covariance matrix for jmax>3 and mmax>2 will be quite difficult and time consuming to check one by one the covariance matrix values. In order to do so, it is recommended to use test.R and running it first for initial covariance matrix. Because the dmnorm function might end up with singularity.
logobj<-0.0
for (k in 1:totaldata){
#cat("lan thu k=")
ik<-inspect1[k] #which equivalent to inspect1
jk<-inspect2[k] # which is second inspection
zk<-interval[k]# which is interval value
xk[1,]<-x0[k]
print(xk[1,])
#xk[2,]<-x1[k]
#xk[3,]<-x2[k]

######################################
#cat("Value of log-likelihood with respective to entire data and obtained unknown parameter beta")
print(k)
#Transforming beta as 
beta<-matrix(beta,nrow=mmax,ncol=(jmax-1))
source("thetavalue.R")
theta<-thetavalue(ik,jk,xk,beta)
print(theta)

thetasa<-thetasavalue(theta)
source("prob.R")
pi<-markovprob(ik,jk,zk,theta,thetasa)

logobj<-logobj+log(pi)

}
print(beta)
logg<-function(beta,mu,sigma){
dmnorm(beta,mu,sigma)
}
print("value cua val")
val<-sum(logobj) +sum(log(logg(beta,mu,sigma)))
print(val)
return(val)
}