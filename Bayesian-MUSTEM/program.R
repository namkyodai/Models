#This program belongs to Nam Lethanh. The author can be contacted via namkyodai@gmail.com
###some initial information for running Markov program#####
data<-read.csv("data.csv",header=TRUE)
attach(data)
totaldata<- 500#length(data[,2])
jmax=max(data[,2]) #Automatically find the maximum condition state in the ranges.
mmax=1 #Number of characteristic variables
#stop("debugging")
# It is noted that in order to run this program, users must select a good fitting covariance matrix   in function mylogposterior, which can be edited by opening file mylogposterior.R. If the value of jmax >3 and mmax>2, it will be difficult to set up covariance matrix because it is likely that the initial covariance matrix you choosed, might end up with singularity when inverting the matrix

#############################
xk<-matrix(double(1),nrow=mmax,ncol=(jmax-1))
beta<-matrix(double(1),nrow=mmax,ncol=(jmax-1))
betaa<-matrix(double(1),nrow=mmax,ncol=(jmax-1))
theta<-matrix(double(1),nrow=1, ncol=jmax)
thetasa<-matrix(double(1),nrow=jmax, ncol=jmax)
##########################################
# SETTING PARAMETER VALUES FOR BAYESIAN AND MCMC
maxinteration<-10000
burnin<-3000
##########################################
##First checking wherethere data of inspections are correct or not
####################################
###############################
######################
#calling function for multi-dimensional normal distribution
source("dmnorm.R")
#Calling Posterior of loglikelihood function.
source("mylogposterior.R")

#ex<-mylogposterior(beta)
# Using laplace function to summary the best starting values of MCMC
source("laplace.R")

#defining starting values for initial values of mean
##Reading initial value of betaa
B1<-0.00015
J1<-0.0105
for (j in 1:(jmax-1)){
  for (m in 1:mmax){
  betaa[m,j]<-j*J1+B1
}
}
#----------------------------------
start<-betaa
print(start)
#stop("debugging")
fit<-laplace(mylogposterior,start,data)
print(fit)
#This is to call log-likelihood values after obtaining values of beta
#calling Metro-Polis Hasting routine
#stop("debugging")
source("rwmetrop.R")
proposal<-list(var=fit$var,scale=2)
start<-betaa
s=rwmetrop(mylogposterior,proposal,start,maxinteration,data)
print(s$accept)
print(apply((s$par),2,quantile,c(.1,.5,.9)))

##THE END###








