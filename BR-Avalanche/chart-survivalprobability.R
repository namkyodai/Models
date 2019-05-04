# Evolution
YearMax<-100
alpha=0.0252
m=1.2665


survivalweib<-function(x,a,b){exp(-(a*x)^b)}
failureweib<-function(x,a,b){a*b*x^(b-1)*exp(-(a*x)^b)}
cdfweib<-function(x,a,b){1-exp(-(a*x)^b)}
x=c(1:YearMax)#"I1-GCS1"

#Bridge
par(mar=c(5, 4, 4, 6) + 0.1)

plot(survivalweib(x,alpha,m),pch=14,axes=FALSE,ylim=c(0,1),ylab="", xlab="",col="red",type="b",lwd=1,cex=0.5)
axis(2, ylim=c(0,1),col="darkblue",las=1)  ## las=1 makes horizontal labels
#mtext(expression(paste("Survival probability")),side=3,line=0.2, adj = -0.06 )
axis(1,pretty(range(c(1:100)),10))
mtext(expression(paste("Time (years)")),side=1,col="black",line=2.5)  
mtext(expression(paste("Reliability")),side=2,col="black",line=2.5)  
box()

#par(new=TRUE)
#plot(cdfweib(x,data$alpha[2], data$m[2]),pch=13,axes=FALSE,ylim=c(0,1),ylab="", xlab="",col="blue",type="o",lwd=1,cex=0.5)

par(new=TRUE)
plot(failureweib(x,alpha,m),pch=13,axes=FALSE,ylim=c(0,0.05),ylab="", xlab="",col="violet",type="o",lwd=1,cex=0.5)
axis(4, ylim=c(0,0.05), col="black",col.axis="black",las=1)
mtext(expression(paste("Density")),side=4,col="black",line=3.5)  
grid(10, 10, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = TRUE)

colors=c("red","violet")

legend("topright",inset=.05,legend=c("Reliability", "Failure density"),text.col="black",lty=c(1,3),pch=c(14,12),col=colors, horiz=F,cex=0.8)
