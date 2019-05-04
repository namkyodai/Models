jmax=5
theta<-c(0.1726,0.0644,0.0744,0.0931,0) # obj1-8
z=1  # Please select the interval or elapsed time in Markov processs
totalyear=50 # Total number of years used to compute the state probabilities over time.
########################################
source("mtp.R") #The source code of function to estimate Markov transition probability
mtp=mtp(jmax,z,theta) #Markov transition probability.

RMD<-matrix(double(1),nrow=1,ncol=(jmax-1))# duration of each condition state
CRMD<-matrix(double(1),nrow=1,ncol=(jmax)) #cummulative life time of condition state
#subroutine to compute the life expectancy based on hazard rate.
for (i in 1:(jmax-1)){
  RMD[i]=1/theta[i] #in years
}
#subroutine to compute the cummulative life expectancy for plotting
for (i in 1:(jmax)){
  if (i==1){
    CRMD[i]<-0
  } else {
    CRMD[i]=CRMD[i-1]+RMD[i-1]
  }
}
#par(mfrow=c(1,2),mar=rep(1, 4))
plot.new() #plot new graph if needed
plot(CRMD[,],c(1:jmax),xlim=c(0,100),ylim=rev(range(c(1,jmax))),xlab="",ylab="",type="p",axes=FALSE,)
predict(smooth.spline(CRMD[,],c(1:jmax),df=4.5)) #function to smooth the curve, df is the smoothness level.
lines(predict(smooth.spline(CRMD[,],c(1:jmax),df=4.5),x=seq(0,50,length=100)),lwd=2)
axis(2, ylim=c(0,5),col="darkblue",las=1)  ## las=1 makes horizontal labels
#title(ylab="Condition states")
mtext(expression(paste("Condition states")),side=2,line=2, adj = 0.5 )

#stop("debugg")


par(new=TRUE)
# Evolution
data=read.csv("failureprobability.csv", header=T)
data=data.frame(data)
YearMax<-100
interventionno<-length(data[,1])
survivalweib<-function(x,a,b){exp(-(a*x)^b)}
failureweib<-function(x,a,b){a*b*x^(b-1)*exp(-(a*x)^b)}
cdfweib<-function(x,a,b){1-exp(-(a*x)^b)}
x=c(1:YearMax)#"I1-GCS1"

#Bridge
par(mar=c(5, 4, 4, 6) + 0.1)

plot(survivalweib(x,data$alpha[1], data$m[1]),pch=14,axes=FALSE,ylim=c(0,1),ylab="", xlab="",col="red",type="b",lwd=1,cex=0.5)
axis(4, ylim=c(0,1),col="darkblue",las=1)  ## las=1 makes horizontal labels
mtext(expression(paste("Survival probability")),side=4,line=2.8, adj = 0.5 )
axis(1,pretty(range(c(1:100)),10))
mtext(expression(paste("Time (years)")),side=1,col="black",line=2)  
colors=c("black","red")
grid(10, 10, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = TRUE)
legend("topright",inset=.08,legend=c("Manifest deterioration", "Latent deterioration"),text.col="black",lty=c(1,1),pch=c(NA,14),lwd=c(2,2),col=colors, horiz=F,cex=1,box.col = "white")

box()





#THE END