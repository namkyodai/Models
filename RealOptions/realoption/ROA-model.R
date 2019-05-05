#This program is developed by Nam Lethanh @IMG, IBI, ETHZ lethanh@ibi.baug.ethz.ch
#Please use this code at your own risk. 
TIME<-30 #total investigated number of years
#-------------------------------DATA---------------------------------
#source("I2-GCS1.R")
source("I3-GCS1.R")
#source("I2-GCS2.R")
#source("I3-GCS2.R")
#source("I2-GCS3.R")
#source("I3-GCS3.R")
#-------------------------------------------------------------------------

#------------------------------MODEL---------------------------------
#The model
N<-matrix(double(1),nrow=1,ncol=TIME)#number of annual visitor 

#For the reference strategy
B1<-matrix(double(1),nrow=1,ncol=TIME)#	revenue 
R1<-matrix(double(1),nrow=1,ncol=TIME)#	routine maintainance 
O1<-matrix(double(1),nrow=1,ncol=TIME)#	Annual Operational cost 
RE1<-matrix(double(1),nrow=1,ncol=TIME)#reliability in year t
FA1<-matrix(double(1),nrow=1,ncol=TIME)#failure probability in year t

#For the investigated strategy
B2<-matrix(double(1),nrow=1,ncol=TIME)#	revenue 
R2<-matrix(double(1),nrow=1,ncol=TIME)#	routine maintainance 
O2<-matrix(double(1),nrow=1,ncol=TIME)#	Annual Operational cost 
RE2<-matrix(double(1),nrow=1,ncol=TIME)#reliability in year t
FA2<-matrix(double(1),nrow=1,ncol=TIME)#failure probability in year t

#common term
S<-matrix(double(1),nrow=1,ncol=TIME)#benefit over time
K<-matrix(double(1),nrow=1,ncol=TIME)#strike price (strike benefit)
h<-matrix(double(1),nrow=1,ncol=TIME)# fare price
Benefit<-matrix(double(1),nrow=1,ncol=TIME)# fare price

#Define a function call call price (similar to European call option)
call.price <- function(x = 1, t = 0, T = 1, r = 1, sigma = 1, K=1){
d2 <- (log(x/K) + (r - 0.5 * sigma^2) * (T - t))/(sigma *sqrt(T - t))
d1 <- d2 + sigma * sqrt(T - t)
x * pnorm(d1) - K * exp(-r * (T - t)) * pnorm(d2)
}


for (T in 1:TIME){
  RE1[T]<-exp(-alpha1*(T**m1)) #reliability function
  FA1[T]<-alpha1*m1*(T**(m1-1))*exp(-alpha1*(T**m1)) #failure probability
  RE2[T]<-exp(-alpha2*(T**m2)) #reliability function
  FA2[T]<-alpha2*m2*(T**(m2-1))*exp(-alpha2*(T**m2)) #failure probability
  if (T==1){
    N[T]=N0
    R2[T]=R02
    O1[T]=O01
    O2[T]=O02
    h[T]=h0
    B1[T]=RE1[T]*(N[T]*h[T]-O1[T])-(1-RE1[T])*W1
    B2[T]=RE2[T]*(N[T]*h[T]-O2[T])-(1-RE2[T])*W2
    S[T]<-B2[T]-R2[T]-C2
    K[T]<-B1[T]-R1[T]
  }  else {
    N[T]=N[T-1]+N[T-1]*0.5/100
    R2[T]=R02
    O2[T]=O02
    h[T]=h[T-1]*exp(mu*T)
    B1[T]=RE1[T]*(N[T]*h[T]-O1[T])-(1-RE1[T])*W1
    B2[T]=RE2[T]*(N[T]*h[T]-O2[T])-(1-RE2[T])*W2
    S[T]<-B2[T]-R2[T]-C2
    K[T]<-B1[T]-R1[T]-C1
    }
  S0 <- S[T]
  K0 <- K[T]
  # t=T-1
  T <- T
  Benefit[T] <- call.price(x = S0, t = 0, T = T, r = r, K = K0,sigma = sigma)
    }
 #end of the model part
print(Benefit)
#print(S)
#print(K)
plot(c(1:TIME), Benefit/scalefactor, main="", pch=16,type="b",xlab="Years ", ylab="Expected benefits (10x6 mu)", lwd=2,col="darkorchid")
par(new=TRUE)
plot(c(1:TIME), S/scalefactor, axes=FALSE,main="", pch=15,type="b",xlab="",ylab="", lwd=2,col="red")
par(new=TRUE)
plot(c(1:TIME), K/scalefactor, axes=FALSE,main="", pch=14,type="b",xlab="",ylab="", lwd=2,col="blue1")
library(MASS)
write.matrix(Benefit, file="benefit.csv",sep=",")
tonghopK=0
for (t in 1:15){
  tonghopK=tonghopK+K[t]
}
tonghopS=0
for (t in 1:15){
  tonghopS=tonghopS+S[t]
}
tonghop=tonghopK+tonghopS
print(tonghop/1000000)

#print(sum(K))

#Maxbenefit
#
for (t in 1:TIME){
Maxbenefit=Benefit[1]
for (t in 1:YearMax){
if (Maxbenefit <= Benefit[t]){
Maxbenefit = Benefit[t]
Opttime=t
}
}
}
cat("optimal intervention time ",t,":",Opttime," years, cost:",Maxbenefit/1000000," year")
