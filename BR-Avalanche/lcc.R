#This program was coded by Nam Lethanh (lethanh@ibi.baug.ethz.ch)
#INPUT

#STEP 4 - calculating cost of corrective intervention for the bridge at each year if there is a failure

CItime<-matrix(double(1),nrow=length(T))  #expect corrective intervention cost at time t
for (t in 1:length(T)){
  CItime[t]<-pi[t,6]*CIimpact
}
cat("Expected IMPACT for CI of the bridge at time t \n")
#print(CItime)

#..............................................

#         MAIN PROGRAM

f1=function(x){
  (alpha1*m*x^(m-1)*exp(-alpha1*x^m))^2
}

f2=function(y){
  (exp(alpha1*y^m-rho*y))
}

q<-matrix(double(1),nrow=1,ncol=T)
C<-matrix(double(1),nrow=1,ncol=T)
canduoi=0
cantren=Inf

#scenario 1 - ***********************************************
ca=cpi*1
cm=cci*1
#library(pracma)
for (i in 1:length(T)){
 # timestep=i/d
  v1=integrate(f1,lower=canduoi,upper=Inf)
  v2=integrate(f2,lower=0,upper=T[i])
  C[i]=(exp(-rho*T[i])*(cs+ca)+(k0/rho)*(1-exp(-rho*T[i]))+(cm+CItime[T[i]])*v1$value*v2$value)/(1-exp(-rho*T[i]))
}




