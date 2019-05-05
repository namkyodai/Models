# Evolution
data=read.csv("survival-failureprobability.csv", header=T)
data=data.frame(data)
YearMax<-30
interventionno<-length(data[,1])
survivalweib<-function(x,a,b){exp(-a*x^b)}
failureweib<-function(x,a,b){a*b*x^(b-1)*exp(-a*x^b)}
x=c(1:YearMax)

#"I1-GCS3"
plot(survivalweib(x,data$alpha[7], data$m[7]),pch=10,axes=FALSE,ylim=c(0,1),ylab="", xlab="",col="darkmagenta",type="b",lwd=2)
axis(2, ylim=c(0,1),col="darkblue",las=1)  ## las=1 makes horizontal labels
mtext(expression(paste("Probability")),side=3,line=0.2, adj = -0.06 )
axis(1,pretty(range(x),10))
mtext(expression(paste("Time (years)")),side=1,col="black",line=2)  
box()

text(10, 0.1,"p.d.f")
text(15, 0.7,"r.f")
#surivival probability

par(new=TRUE)
#"I2-GCS3"
plot(survivalweib(x,data$alpha[8], data$m[8]),pch=9,axes=FALSE,ylim=c(0,1),ylab="", xlab="",col="darkmagenta",type="b",lwd=2)
par(new=TRUE)
#"I3-GCS3"
plot(survivalweib(x,data$alpha[9], data$m[9]),pch=8,axes=FALSE,ylim=c(0,1),ylab="", xlab="",col="darkmagenta",type="b",lwd=2)


#failure probability
par(new=TRUE)
#"I1-GCS3"
plot(failureweib(x,data$alpha[7], data$m[7]),pch=10,axes=FALSE,ylim=c(0,1),ylab="", xlab="",col="darkmagenta",type="b",lwd=1)
par(new=TRUE)

#"I2-GCS3"
plot(failureweib(x,data$alpha[8], data$m[8]),pch=9,axes=FALSE,ylim=c(0,1),ylab="", xlab="",col="darkmagenta",type="b",lwd=1)
par(new=TRUE)

#"I3-GCS3"
plot(failureweib(x,data$alpha[9], data$m[9]),pch=8,axes=FALSE,ylim=c(0,1),ylab="", xlab="",col="darkmagenta",type="b",lwd=1)


## Add Legend
legend("topright",inset=.05,legend=c("I1-GCS3","I2-GCS3","I3-GCS3"),text.col=c("darkmagenta","darkmagenta","darkmagenta"),pch=c(10,9,8),col=c("darkmagenta","darkmagenta","darkmagenta"),horiz=F)
#legend("topright",inset=.02,legend=c("I1-GCS1","I2-GCS1","I3-GCS1","I1-GCS2","I2-GCS2","I3-GCS2","I1-GCS3","I2-GCS3","I3-GCS3"),text.col=c("blue1","blue1","blue1","red","red","red","darkmagenta","darkmagenta","darkmagenta"),pch=c(16,15,14,13,12,11,10,9,8),col=c("blue1","blue1","blue1","red","red","red","darkmagenta","darkmagenta","darkmagenta"),horiz=F)

#plot the cdf




