#this routine is used to calculated the Intensity of rockfall, avalanche using the concept of fragility curve

#Created by Nam Lethanh (ETH Zürich, namkyodai (at) gmail (dot) com )

#INPUT -----------------------------------
Maxvolume <- 34 # e.g. m3 (volume of rockfall or avalanche)
Time<-100
MaxCS<-5 #maximum no of CS for the bridge
Maxlatent<- 1 #number of latent condition state
alpha<-c(0.353,0.360,0.371,0.389,0.391)
gamma<-c(5.652,5.581,5.575,5.498,5.367)
#define dimension
mtp<- matrix(double(1),nrow=MaxCS,ncol=(MaxCS+Maxlatent))
mtplatent<- matrix(double(1),nrow=(MaxCS+Maxlatent),ncol=(MaxCS+Maxlatent))
mtp <- read.csv("mtp.csv",header=FALSE) 
#define the Markov transition probability for the bridge
attach(mtp)
initialCS<-c(1,0,0,0,0,0) #initial state probability of the bridge
v1=10 #value of rock fall, or avalanche volume to be used

#--------END OF INPUT--------------------------
#---------Start caculation SECTION
#STEP 1 --estimating the failure probability from hazard curve for each CS
#hazard curve
S<-matrix(double(1),nrow=Maxvolume,ncol=MaxCS)
#define a function for calculating intensity S
intensity<-function(alpha,gamma,volume){
  exp(-exp(-alpha*volume+gamma))
}
for (i in 1:MaxCS){
for (v in 1:Maxvolume){
  S[v,i]<-intensity(alpha[i],gamma[i],v)
}
}
cat("intensity - hazard curve \n")
print(S)

#STEP 2 --Calculating the Markov transition probability - Latent
#choose the intensity (e.g. S=10 m3)

for (i in 1:(MaxCS+1)){
  if (i < (MaxCS+1)) {
  for (j in 1:(MaxCS+Maxlatent)){
    if (j < (MaxCS+Maxlatent)){
      mtplatent[i,j]<-mtp[i,j]*(1-S[v1,i])
    } else {
    mtplatent[i,j]<-S[v1,i]
    }
  }
  } else {
    mtplatent[i,j]<-1
  }
}
cat("Manifest deterioration -Markov transition matrix \n")
print(mtp)
cat("Latent deterioration -Markov transition matrix \n")
print(mtplatent)

#----------------------------------------------------
#STEP 3. Define the initial state of the bridge at the time after the PI
pi<-matrix(double(1),nrow=Time,ncol=(MaxCS+Maxlatent))
for (t in 1:Time){
  if (t ==1){
    pi[t,]<-initialCS
  } else {
        pi[t,]<-pi[t-1,]%*%mtplatent
    }
}
cat("Latent deterioration (no intervention) -State probability \n")
print(pi)

scalefactor=1
trucy=1
colors=c("red", "darkorchid","green", "violet", "darkorange3","brown3")

plot(S[,1]/scalefactor,pch=16,axes=FALSE,ylim=c(0,trucy),ylab="", xlab="",col="red",type="o",lwd=4)
axis(2, ylim=c(0,trucy),col="darkblue",las=1)
mtext(expression(paste("Probability of failure CS l")),side=2,line=2.2, adj = 0.5 )
axis(1,pretty(range(c(1:Maxvolume)),10))
mtext(expression(paste("Rock fall volume" (m^3))),side=1,col="black",line=2.2)  
box()

par(new=TRUE)
plot(S[,2]/scalefactor,pch=15,axes=FALSE,ylim=c(0,trucy),ylab="", xlab="",col="darkorchid",type="o",lwd=4)
par(new=TRUE)
plot(S[,3]/scalefactor,pch=14,axes=FALSE,ylim=c(0,trucy),ylab="", xlab="",col="green",type="o",lwd=4)
grid(10, 10, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)

legend("topleft",inset=.02,legend=c("i=1","i=2","i=3"),text.col=colors,pch=c(16,15,14),col=colors,horiz=F,bty != "n",bg="white")
text(15, 18.5,expression(paste("I3-GCS1 (14.73x",10^6, "mu-year 15)")),pos=1,cex = 0.8, srt = 0)


