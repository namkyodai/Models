#This program was coded by Nam Lethanh (lethanh@ibi.baug.ethz.ch)
#INPUT
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
rho=0.02
m=1.2665
cs=0
k0=0

Maxvolume <- 34 # e.g. m3 (volume of rockfall or avalanche)
Time<-YearMax
MaxCS<-5 #maximum no of CS for the bridge
Maxlatent<- 1 #number of latent condition state

#define dimension
mtp<- matrix(double(1),nrow=MaxCS,ncol=(MaxCS+Maxlatent))
mtponeyear<- matrix(double(1),nrow=MaxCS,ncol=(MaxCS+Maxlatent))
mtplatent<- matrix(double(1),nrow=(MaxCS+Maxlatent),ncol=(MaxCS+Maxlatent))
mtplatentoneyear<- matrix(double(1),nrow=(MaxCS+Maxlatent),ncol=(MaxCS+Maxlatent))
mtp <- read.csv("mtp.csv",header=FALSE) 
mtponeyear <- read.csv("mtponeyear.csv",header=FALSE) 
#define the Markov transition probability for the bridge
#attach(mtp)
initialCS<-c(1,0,0,0,0,0) #initial state probability of the bridge

v3=10 #value of rock fall, or avalanche volume to be used

#---------Start caculation SECTION
#STEP 1 --estimating the failure probability from hazard curve for each CS
#hazard curve
S<-matrix(double(1),nrow=Maxvolume,ncol=MaxCS)
Soneyear<-matrix(double(1),nrow=Maxvolume,ncol=MaxCS)
#define a function for calculating intensity S
intensity<-function(alpha,gamma,volume){
  exp(-exp(-alpha*volume+gamma))
}
for (i in 1:MaxCS){
  for (v in 1:Maxvolume){
    S[v,i]<-intensity(alpha[i],gamma[i],v)/d
    Soneyear[v,i]<-intensity(alpha[i],gamma[i],v)
  }
}

#STEP 2 --Calculating the Markov transition probability - Latent
#choose the intensity (e.g. S=10 m3)

for (i in 1:(MaxCS+1)){
  if (i < (MaxCS+1)) {
    for (j in 1:(MaxCS+Maxlatent)){
      if (j < (MaxCS+Maxlatent)){
        mtplatent[i,j]<-mtp[i,j]*(1-S[v3,i])
        mtplatentoneyear[i,j]<-mtponeyear[i,j]*(1-Soneyear[v3,i])
      } else {
        mtplatent[i,j]<-S[v3,i]
        mtplatentoneyear[i,j]<-Soneyear[v3,i]
      }
    }
  } else {
    mtplatent[i,j]<-1
    mtplatentoneyear[i,j]<-1
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

for (t in 1:Time){
  if (t ==1){
    pioneyear[t,]<-initialCS
  } else {
    pioneyear[t,]<-pioneyear[t-1,]%*%mtplatentoneyear
  }
}

source("stackedPlot.R")

pioneyear=pioneyear*100
colors=c("aliceblue","coral1","limegreen","orange","yellow2","black")
stackedPlot(data.frame(pioneyear[,]),col=colors,xlab="Time (years)",ylab="Percentage(%)")
legend("topright", inset=0.09, title="CSs",col=colors,lty=2,lwd=13,legend=c(1:(MaxCS+Maxlatent)),bg="azure2",cex=1.2)
title(main="", col.main="red", font.main=4)

#THE END



#---Case 1
cpi<-1 #preventive intervention cost for the snow catcher system
cci<-2 #corrective intervention cost for the snow catcher system
CIimpact<-10 #corrective intervention cost for the bridge

source("lcc.R")
plot.new()
par(mar=c(5, 4, 4, 6) + 0.1)
plot(T,C/factor, axes=FALSE, ylim=c(0,yvalue), xlab="", ylab="", type="l",col="red", lwd=4)

axis(2, ylim=c(0,max(C)/factor),col="black",las=1)  ## las=1 makes horizontal labels
mtext(expression(paste("Impact (mus)" )),side=2,col="black",line=2.3)
box()
axis(1,pretty(range(seq(T)/d),d))
mtext(expression(paste("Time (years)")),side=1,col="black",line=2.5)  
grid(10, 10, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = TRUE)

Min_C<-min(C[c(10:length(T))])
T_optimal<-which.min(C)/d

cat("lowest total impact is ")
print(Min_C)
cat("Optimal preventive intervention time")
print(T_optimal)

points(T_optimal,Min_C/factor,pch=23,bg="red",lwd=3)
text(T_optimal*0.9, Min_C/factor*1.2,"A(19, 6.59)",pos=4,cex = 1, srt = 0)

#---Case 2
cpi<-1 #preventive intervention cost for the snow catcher system
cci<-4 #corrective intervention cost for the snow catcher system
CIimpact<-10 #corrective intervention cost for the bridge
#-------case 2
source("lcc.R")
par(new=TRUE)
par(mar=c(5, 4, 4, 6) + 0.1)
plot(T,C/factor, axes=FALSE, ylim=c(0,yvalue), xlab="", ylab="", type="l",col="blueviolet", lwd=2)

Min_C<-min(C[c(10:length(T))])
T_optimal<-which.min(C)/d

cat("lowest total impact is ")
print(Min_C)
cat("Optimal preventive intervention time")
print(T_optimal)

points(T_optimal,Min_C/factor,pch=23,bg="red",lwd=3)
text(T_optimal*0.9, Min_C/factor*1.2,"B(15, 10.64)",pos=4,cex = 1, srt = 0)
#-----------case 3
#---Case 2
cpi<-1 #preventive intervention cost for the snow catcher system
cci<-1 #corrective intervention cost for the snow catcher system
CIimpact<-10 #corrective intervention cost for the bridge

source("lcc.R")
par(new=TRUE)
par(mar=c(5, 4, 4, 6) + 0.1)
plot(T,C/factor, axes=FALSE, ylim=c(0,yvalue), xlab="", ylab="", type="l",col="blue", lwd=2)

Min_C<-min(C[c(10:length(T))])
T_optimal<-which.min(C)/d
cat("lowest total impact is ")
print(Min_C)
cat("Optimal preventive intervention time")
print(T_optimal)

library("lattice")
points(T_optimal,Min_C/factor,pch=23,bg="red",lwd=3)
text(T_optimal*0.9, Min_C/factor*1.2,"C(23, 4.26)",pos=4,cex = 1, srt = 0)

legend("topright",inset=.05,legend=c("Scenario 1", "Scenario 2", "Scenario 3"),text.col=c("black"),lty=c(1,1,1),lwd=c(3,1,1),col=c("blueviolet","red","blue"),horiz=F,cex=1,box.col = "white")

library(MASS)
print(mtplatentoneyear)
write.matrix(mtplatentoneyear, file="mtplatentoneyear.csv",sep=",")

#Start the sensitivity analysis calculation
#stop("debugg")

#ratio CI/PI snow barrier system
source("sensitivityanalysis-ratioCIPI-snow.R")
source("sensitivityanalysis-ratioCIPI-bridge.R")
source("sensitivityanalysis-alpha.R")
source("sensitivityanalysis-m.R")
source("sensitivityanalysis-avalancheintensity.R") #age of the bridge
source("sensitivityanalysis-CS.R") #age of the bridge
source("sensitivityanalysis-rho.R") #discount factor
















