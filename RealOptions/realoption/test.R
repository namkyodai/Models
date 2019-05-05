#This is real option program used for determining the optimal timing to carry out intervention
require(fOptions) # pls install package fOptions before running.
#function for estimate the call price
call.price <- function(x = 1, t = 0, T = 1, r = 1, sigma = 1, K=1){
d2 <- (log(x/K) + (r - 0.5 * sigma^2) * (T - t))/(sigma *sqrt(T - t))
d1 <- d2 + sigma * sqrt(T - t)
x * pnorm(d1) - K * exp(-r * (T - t)) * pnorm(d2)
}

W=350000
alpha1=0.002
m1=2.3
alpha2=0.002
m2=2.05
R01=206000 #annual routine maintenance
O01=727000 #annual operational cost
R02=237000 #annual routine maintenance
O02=745000 #annual operational cost
N0=200000 #annual visitor at time 0
C2<-2184000 #intervention cost
C1<-0 #660000 #intervention cost
h0=180 # CHF
mu=0.0003
#defining the total time to estimate
TIME=30 #years


B1<-matrix(double(1),nrow=1,ncol=TIME)#revenue
B2<-matrix(double(1),nrow=1,ncol=TIME)#revenue
N<-matrix(double(1),nrow=1,ncol=TIME)#number of annual visitor 
R1<-matrix(double(1),nrow=1,ncol=TIME)#routine maintainance
R2<-matrix(double(1),nrow=1,ncol=TIME)#routine maintainance
O1<-matrix(double(1),nrow=1,ncol=TIME)#Annual Operational cost
O2<-matrix(double(1),nrow=1,ncol=TIME)#Annual Operational cost
RE1<-matrix(double(1),nrow=1,ncol=TIME)#reliability in year t
RE2<-matrix(double(1),nrow=1,ncol=TIME)#reliability in year t
FA1<-matrix(double(1),nrow=1,ncol=TIME)#failure probability in year t
FA2<-matrix(double(1),nrow=1,ncol=TIME)#failure probability in year t
S<-matrix(double(1),nrow=1,ncol=TIME)#benefit over time
K<-matrix(double(1),nrow=1,ncol=TIME)#strike price (strike benefit)
h<-matrix(double(1),nrow=1,ncol=TIME)# fare price
Benefit<-matrix(double(1),nrow=1,ncol=TIME)# fare price

for (T in 1:TIME){
  RE1[T]<-exp(-alpha1*(T**m1)) #reliability function
  FA1[T]<-alpha1*m1*(T**(m1-1))*exp(-alpha1*(T**m1)) #failure probability
  RE2[T]<-exp(-alpha2*(T**m2)) #reliability function
  FA2[T]<-alpha2*m2*(T**(m2-1))*exp(-alpha2*(T**m2)) #failure probability
  if (T==1){
    N[T]=N0
    R2[T]=R02
    O2[T]=O02
    h[T]=h0
    B1[T]=RE1[T]*(N[T]*h[T]-O1[T])-(1-RE1[T])*W
    B2[T]=RE2[T]*(N[T]*h[T]-O2[T])-(1-RE2[T])*W
    S[T]<-B2[T]-R2[T]-C2
    K[T]<-B1[T]-R1[T]
  }  else {
    N[T]=N[T-1]+N[T-1]*0.5/100
    R2[T]=R02
    O2[T]=O02
    h[T]=h[T-1]*exp(mu*T)
    B1[T]=RE1[T]*(N[T]*h[T]-O1[T])-(1-RE1[T])*W
    B2[T]=RE2[T]*(N[T]*h[T]-O2[T])-(1-RE2[T])*W
    S[T]<-B2[T]-R2[T]-C2
    K[T]<-B1[T]-R1[T]-C1
    }
  S0 <- S[T]
  K0 <- K[T]
  r <- 0.02
 # t=T-1
  T <- T
  sigma <- 0.2
  Benefit[T] <- call.price(x = S0, t = 0, T = T, r = r, K = K0,sigma = sigma)
  }
print(Benefit)
plot(c(1:30), Benefit/1000, main="", pch=16,type="b",xlab="Years ", ylab="Expected benefits (1000 mu)", lwd=2,col="darkorchid")
library(MASS)
write.matrix(Benefit, file="benefit.csv",sep=",")

