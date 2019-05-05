# Evolution
data=read.csv("survival-failureprobability.csv", header=T)
data=data.frame(data)
YearMax<-30
interventionno<-length(data[,1])
survivalweib<-function(x,a,b){exp(-a*x^b)}
failureweib<-function(x,a,b){a*b*x^(b-1)*exp(-a*x^b)}
x=c(1:YearMax)

#"I2-GCS1"
plot(survivalweib(x,data$alpha[4], data$m[4]),pch=13,axes=FALSE,ylim=c(0,1),ylab="", xlab="",col="red",type="o",lwd=2)
axis(2, ylim=c(0,1),col="darkblue",las=1)  ## las=1 makes horizontal labels
mtext(expression(paste("Probability")),side=3,line=0.2, adj = -0.06 )
axis(1,pretty(range(x),10))
mtext(expression(paste("Time (years)")),side=1,col="black",line=2)  
box()

text(10, 0.1,"p.d.f")
text(15, 0.7,"r.f")


par(new=TRUE)
#"I2-GCS2"
plot(survivalweib(x,data$alpha[5], data$m[5]),pch=14,axes=FALSE,ylim=c(0,1),ylab="", xlab="",col="red",type="o",lwd=2)
par(new=TRUE)
#"I3-GCS2"
plot(survivalweib(x,data$alpha[6], data$m[6]),pch=11,axes=FALSE,ylim=c(0,1),ylab="", xlab="",col="red",type="o",lwd=2)

#failure probability
par(new=TRUE)
#"I1-GCS2"
plot(failureweib(x,data$alpha[4], data$m[4]),pch=13,axes=FALSE,ylim=c(0,1),ylab="", xlab="",col="red",type="o",lwd=0.2)
par(new=TRUE)
#"I2-GCS2"
plot(failureweib(x,data$alpha[5], data$m[5]),pch=14,axes=FALSE,ylim=c(0,1),ylab="", xlab="",col="red",type="o",lwd=0.2)
par(new=TRUE)
#"I3-GCS2"
plot(failureweib(x,data$alpha[6], data$m[6]),pch=11,axes=FALSE,ylim=c(0,1),ylab="", xlab="",col="red",type="o",lwd=0.2)



## Add Legend
legend("topright",inset=.05,legend=c("I1-GCS2","I2-GCS2","I3-GCS2"),text.col=c("red","red","red"),pch=c(13,12,11),col=c("red","red","red"),horiz=F)
#legend("topright",inset=.02,legend=c("I1-GCS1","I2-GCS1","I3-GCS1","I1-GCS2","I2-GCS2","I3-GCS2","I1-GCS3","I2-GCS3","I3-GCS3"),text.col=c("blue1","blue1","blue1","red","red","red","darkmagenta","darkmagenta","darkmagenta"),pch=c(16,15,14,13,12,11,10,9,8),col=c("blue1","blue1","blue1","red","red","red","darkmagenta","darkmagenta","darkmagenta"),horiz=F)

#plot the cdf
