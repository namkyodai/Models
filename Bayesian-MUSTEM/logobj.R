#this function is a part coding of Markov program, which can be referred to the doctoral thesis of Nam Lethanh. Users bear their own this if applying this part of coding. I highly recommend you to cite this part of my work if you wish to use in your writing and work. Thank you.
########## this is log likelihood function in Markov program ############################
logobj<-function(beta){
###Reading characteristic variables
logobj<-0.0
for (k in 1:totaldata){
ik<-inspect1[k] #which equivalent to inspect1
jk<-inspect2[k] # which is second inspection
zk<-interval[k]# which is interval value
xk[1,]<-x0[k]
#xk[2,]<-x1[k]
#xk[3,]<-x2[k]
######################################

#cat("Value of log-likelihood with respective to entire data and obtained unknown parameter beta")
#print(k)
source("thetavalue.R")
theta<-thetavalue(ik,jk,xk,beta)
thetasa<-thetasavalue(theta)
#print(theta)
#print(thetasa)
source("prob.R")
pi<-markovprob(ik,jk,zk,theta,thetasa)
logobj<-logobj+log(pi) #This is the log-likelihood of our objective function
#print(logobj)
}
return(logobj)
}

##################################################



