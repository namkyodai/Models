# Evolution
data=read.csv("survival-failureprobability.csv", header=T)
data=data.frame(data)
YearMax<-30
interventionno<-length(data[,1])
survivalweib<-function(x,a,b){exp(-a*x^b)}
failureweib<-function(x,a,b){a*b*x^(b-1)*exp(-a*x^b)}
x=c(1:YearMax)
#"I1-GCS1"
plot(survivalweib(x,data$alpha[1], data$m[1]),pch=16,axes=FALSE,ylim=c(0,1),ylab="", xlab="",col="blue1",type="b",lwd=2)
axis(2, ylim=c(0,1),col="darkblue",las=1)  ## las=1 makes horizontal labels
mtext(expression(paste("Survival probability")),side=3,line=0.2, adj = -0.06 )
axis(1,pretty(range(c(1:100)),10))
mtext(expression(paste("Time (years)")),side=1,col="black",line=2)  
box()

#text(15, 0.7,"r.f")
#par(new=TRUE)
#"I2-GCS1"
#plot(survivalweib(x,data$alpha[2], data$m[2]),pch=15,axes=FALSE,ylim=c(0,1),ylab="", xlab="",col="blue1",type="b",lwd=2)
#par(new=TRUE)
#"I3-GCS1"
#plot(survivalweib(x,data$alpha[3], data$m[3]),pch=14,axes=FALSE,ylim=c(0,1),ylab="", xlab="",col="blue1",type="b",lwd=2)
#par(new=TRUE)
#plot(survivalweib(x,data$alpha[4], data$m[4]),pch=13,axes=FALSE,ylim=c(0,1),ylab="", xlab="",col="red",type="o",lwd=2)
#par(new=TRUE)
#plot(survivalweib(x,data$alpha[7], data$m[7]),pch=10,axes=FALSE,ylim=c(0,1),ylab="", xlab="",col="darkmagenta",type="b",lwd=2)

#par(new=TRUE)
#failure probability
#plot(failureweib(x,data$alpha[1], data$m[1]),pch=16,axes=FALSE,ylim=c(0,1),ylab="", xlab="",col="blue1",type="b",lwd=1)
#par(new=TRUE)
#plot(failureweib(x,data$alpha[2], data$m[2]),pch=15,axes=FALSE,ylim=c(0,1),ylab="", xlab="",col="blue1",type="b",lwd=1)
#par(new=TRUE)
#"I3-GCS1"
#plot(failureweib(x,data$alpha[3], data$m[3]),pch=14,axes=FALSE,ylim=c(0,1),ylab="", xlab="",col="blue1",type="b",lwd=1)
#par(new=TRUE)

## Add Legend
#legend("bottomleft",inset=.05,legend=c("I1-GCS1","I2-GCS1","I3-GCS1"),text.col=c("blue1","blue1","blue1"),pch=c(16,15,14),col=c("blue1","blue1","blue1"),horiz=F)
#
#legend("bottomleft",inset=.05,legend=c("I1-GCS1","I2-GCS1","I3-GCS1","I1-GCS2","I2-GCS2","I3-GCS2","I1-GCS3","I2-GCS3","I3-GCS3"),text.col=c("blue1","blue1","blue1","red","red","red","darkmagenta","darkmagenta","darkmagenta"),pch=c(16,15,14,13,12,11,10,9,8),col=c("blue1","blue1","blue1","red","red","red","darkmagenta","darkmagenta","darkmagenta"), horiz=F)
#segments(13, 0.55, 13, 1, col= 'blue1')
grid(10, 10, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = TRUE)
legend("topright",inset=.05,legend=c("Do nothing","DS1","DS2"),text.col=c("blue1","red","darkmagenta"),pch=c(16,15,14),col=c("blue1","red","darkmagenta"), horiz=FALSE)
#segments(13, 0.55, 13, 1, col= 'blue1')

text(12.8, 0.95,"Execution of interventions",pos=2,cex = 0.9, srt = 90)
arrows(13, 0.3, 13, 0.98,length = 0.2, angle = 20,code = 2,col= 'burlywood4')
#par(new=TRUE)
#"I2-GCS1"
y=c(13:30)
survivalweib2<-function(y,a,b){exp(-a*(y-12)^b)}
#plot(survivalweib2(y,data$alpha[2], data$m[2]),pch=15,axes=FALSE,ylim=c(0,1),ylab="", xlab="",col="red",type="b",lwd=2)
#GCS1
lines(y,survivalweib2(y,data$alpha[2], data$m[2]),pch=15,col="red",type="b",lwd=1)
lines(y,survivalweib2(y,data$alpha[3], data$m[3]),pch=14,col="darkmagenta",type="b",lwd=1)
#GCS2
#segments(15, 0.44, 15, 1, col= 'red')
#arrows(15, 0.3, 15, 0.88,length = 0.2, angle = 20,code = 2,col= 'red')
#text(14.8, 0.88,"Intervention on GCS2",pos=2,cex = 0.9, srt = 90)
#y=c(15:30)
#survivalweib2<-function(y,a,b){exp(-a*(y-14)^b)}
#lines(y,survivalweib2(y,data$alpha[5], data$m[5]),pch=12,col="red",type="b",lwd=1)
#lines(y,survivalweib2(y,data$alpha[6], data$m[6]),pch=11,col="red",type="b",lwd=1)
#GCS3
#segments(18, 0.3, 18, 1, col= 'darkmagenta')
#arrows(18, 0.15, 18, 0.88,length = 0.2, angle = 20,code = 2,col= 'darkmagenta')
#text(17.8, 0.82,"Intervention on GCS3",pos=2,cex = 0.9, srt = 90)
#y=c(18:30)
#survivalweib2<-function(y,a,b){exp(-a*(y-17)^b)}
#lines(y,survivalweib2(y,data$alpha[8], data$m[8]),pch=9,col="darkmagenta",type="b",lwd=1)
#lines(y,survivalweib2(y,data$alpha[9], data$m[9]),pch=8,col="darkmagenta",type="b",lwd=1)