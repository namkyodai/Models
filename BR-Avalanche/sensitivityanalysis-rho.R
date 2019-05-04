#Sensitivity analysis
yvalue = 30
alpha<-c(0.353,0.360,0.371,0.389,0.391)
gamma<-c(5.652,5.581,5.575,5.498,5.367)
d=10 # factor of time
YearMax<-50 #total no of year to be investigated.
#creating local variables
rate<-matrix(nrow=YearMax)     #hazard rate in respective year
colors=c("chocolate","coral3","coral4","cyan1","cyan4","darkblue","blue","chartreuse4","red","deeppink","deeppink4") #colors for plot
#SUBROUTINE FOr the bridge
#this routine is used to calculated the Intensity of rockfall, avalanche
#INPUT -----------------------------------
T=seq(0.1,YearMax, by =0.1)
#with example
#library(cubature)
factor=1
alpha1=0.0252
#rho=0.02
m=1.2665
cs=0
k0=0

Maxvolume <- 34 # e.g. m3 (volume of rockfall or avalanche)
Time<-YearMax
MaxCS<-5 #maximum no of CS for the bridge
Maxlatent<- 1 #number of latent condition state


#define dimension
mtp<- matrix(double(1),nrow=MaxCS,ncol=(MaxCS+Maxlatent))
mtplatent<- matrix(double(1),nrow=(MaxCS+Maxlatent),ncol=(MaxCS+Maxlatent))
mtp <- read.csv("mtp.csv",header=FALSE) 
#define the Markov transition probability for the bridge
#attach(mtp)
initialCS<-c(1,0,0,0,0,0) #initial state probability of the bridge
v3=10 #value of rock fall, or avalanche volume to be used

#STEP 1 --estimating the failure probability from hazard curve for each CS
#hazard curve
S<-matrix(double(1),nrow=Maxvolume,ncol=MaxCS)
#define a function for calculating intensity S
intensity<-function(alpha,gamma,volume){
  exp(-exp(-alpha*volume+gamma))
}
for (i in 1:MaxCS){
  for (v in 1:Maxvolume){
    S[v,i]<-intensity(alpha[i],gamma[i],v)/d
  }
}
#STEP 2 --Calculating the Markov transition probability - Latent
#choose the intensity (e.g. S=10 m3)
for (i in 1:(MaxCS+1)){
  if (i < (MaxCS+1)) {
    for (j in 1:(MaxCS+Maxlatent)){
      if (j < (MaxCS+Maxlatent)){
        mtplatent[i,j]<-mtp[i,j]*(1-S[v3,i])
      } else {
        mtplatent[i,j]<-S[v3,i]
      }
    }
  } else {
    mtplatent[i,j]<-1
  }
}

#----------------------------------------------------
#STEP 3. Define the initial state of the bridge at the time after the PI
pi<-matrix(double(1),nrow=length(T),ncol=(MaxCS+Maxlatent))
pioneyear<-matrix(double(1),nrow=Time,ncol=(MaxCS+Maxlatent))
for (t in 1:length(T)){
  if (t ==1){
    pi[t,]<-initialCS
  } else {
    pi[t,]<-pi[t-1,]%*%mtplatent
  }
}

#---------Start caculation SECTION

#Sensitivity analysis for ratio CI/PI
cpi<-1 #preventive intervention cost for the snow catcher system
CIimpact<-10 #corrective intervention cost for the bridge
cci<-2
maxrho<- 0.1
step=seq(0.01,maxrho, by =0.001)
rhosa<- matrix(double(1),nrow=length(step),ncol=1)
T_optimalsa<-matrix(double(1),nrow=length(step),ncol=1)
Min_Csa<-matrix(double(1),nrow=length(step),ncol=1)

for (u in 1:length(step)){
  rhosa[u]<-step[u]
#This program was coded by Nam Lethanh (lethanh@ibi.baug.ethz.ch)
#INPUT

#STEP 4 - calculating cost of corrective intervention for the bridge at each year if there is a failure

CItime<-matrix(double(1),nrow=length(T))  #expect corrective intervention cost at time t
for (t in 1:length(T)){
  CItime[t]<-pi[t,6]*CIimpact
}
#cat("Expected IMPACT for CI of the bridge at time t \n")
#print(CItime)

#..............................................

rho<-rhosa[u]
  
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
cm=cci
for (i in 1:length(T)){
  # timestep=i/d
  v1=integrate(f1,lower=canduoi,upper=Inf)
  v2=integrate(f2,lower=0,upper=T[i])
  C[i]=(exp(-rho*T[i])*(cs+ca)+(k0/rho)*(1-exp(-rho*T[i]))+(cm+CItime[T[i]])*v1$value*v2$value)/(1-exp(-rho*T[i]))
}

Min_Csa[u]<-min(C[c(10:length(T))])
T_optimalsa[u]<-which.min(C)/d
}

#drawing
plot.new()
par(mar=c(5, 4, 4, 6) + 0.1)
plot(rhosa, Min_Csa, pch=2, axes=FALSE, ylim=c(0,15), xlab="", ylab="", 
     type="l",col="blueviolet", lwd=5)
axis(2, ylim=c(0,20),col="black",las=1)  ## las=1 makes horizontal labels
mtext(expression(paste("Annual impacts  (",, "mus)")),side=2,col="black",line=2.5)
axis(1,pretty(range(rhosa),20))
box()

par(new=TRUE)

## Plot the second plot and put axis scale on right
plot(rhosa, T_optimalsa, pch=4,  xlab="", ylab="", ylim=c(15,25), 
     axes=FALSE, type="l", col="coral1",lwd=5,lty=2)
## a little farther out (line=4) to make room for labels
mtext("OTE (years)",side=4,col="black",line=2.3) 
axis(4, ylim=c(0,20), col="black",col.axis="black",las=1)
mtext(expression(paste("Discount factor ", rho)),side=1,col="black",line=2.5)  
grid(10, 10, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = TRUE)
legend("topright",inset=.08,legend=c("Annual impacts","OTE"),
       text.col=c("black"),lty=c(1,2),lwd=c(5,5),col=c("blueviolet","coral1"),horiz=F,cex=1.2,box.col = "white")

