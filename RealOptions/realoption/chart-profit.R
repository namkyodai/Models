# Evolution
data=read.csv("profit.csv", header=T)
attach(data)
#object 1
scalefactor=1000000
trucy=18
colors=c("red", "darkorchid","green", "violet", "darkorange3","brown3")

plot(data$DS1/scalefactor,pch=16,axes=FALSE,ylim=c(0,trucy),ylab="", xlab="",col="red",type="o",lwd=1)
axis(2, ylim=c(0,trucy),col="darkblue",las=1)
mtext(expression(paste("Option value  (", 10^6," mu)")),side=3,line=0.2, adj = -0.06 )
axis(1,pretty(range(data$Time),10))
mtext(expression(paste("Time (years)")),side=1,col="black",line=2)  
box()

par(new=TRUE)
plot(data$DS2/scalefactor,pch=15,axes=FALSE,ylim=c(0,trucy),ylab="", xlab="",col="darkorchid",type="o",lwd=1)
#par(new=TRUE)
#plot(data$I2.GCS2/scalefactor,pch=14,axes=FALSE,ylim=c(0,trucy),ylab="", xlab="",col="green",type="o",lwd=1)
#par(new=TRUE)
#plot(data$I3.GCS2/scalefactor,pch=13,axes=FALSE,ylim=c(0,trucy),ylab="", xlab="",col="violet",type="o",lwd=1)
#par(new=TRUE)
#plot(data$I2.GCS3/scalefactor,pch=12,axes=FALSE,ylim=c(0,trucy),ylab="", xlab="",col="darkorange3",type="o",lwd=1)
#par(new=TRUE)
#plot(data$I3.GCS3/scalefactor,pch=11,axes=FALSE,ylim=c(0,trucy),ylab="", xlab="",col="brown3",type="o",lwd=1)
grid(10, 10, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = TRUE)
legend("topright",inset=.02,legend=c("DS1","DS2"),text.col=colors,pch=c(16,15),col=colors,horiz=F)

arrows(15, 15.2, 15, 17,length = 0.15, angle = 20,code = 1,col= 'blue3',lty = 5,lwd=1)
text(15, 18.5,expression(paste("DS2 (14.73x",10^6, "mu-year 15)")),pos=1,cex = 0.8, srt = 0)
text(14.5, 16.2,"A",pos=1,cex = 0.8, srt = 0)
points(15,14.72,pch=23,bg="black",lwd=3)
segments(15, 0, 15, 14.72, col= 'cadetblue',lty=2,lwd=0.5)
segments(0, 14.72, 15, 14.72, col= 'cadetblue',lty=2,lwd=0.5)

text(19.5, 11.5,"B",pos=1,cex = 0.8, srt = 0)
points(20,11.24,pch=23,bg="black",lwd=3)
segments(20, 0, 20, 11.24, col= 'cadetblue',lty=2,lwd=0.5)
segments(0, 11.24, 20, 11.24, col= 'cadetblue',lty=2,lwd=0.5)

text(6.5, 10,"C",pos=1,cex = 0.8, srt = 0)
points(7,8.9,pch=23,bg="black",lwd=3)
segments(7, 0, 7, 8.9, col= 'cadetblue',lty=2,lwd=0.5)
segments(0, 8.9, 7, 8.9, col= 'cadetblue',lty=2,lwd=0.5)


#arrows(16, 13.9, 18, 16,length = 0.15, angle = 20,code = 1,col= 'blue3',lty = 5,lwd=1)
#text(21, 17.4,expression(paste("I3-GCS2 (13.93x",10^6, "mu-year 16)")),pos=1,cex = 0.8, srt = 0)
#points(16, 13.9,pch=23,bg="chartreuse4",lwd=1)
#arrows(17, 12.7, 20, 14.5,length = 0.15, angle = 20,code = 1,col= 'blue3',lty = 5,lwd=1)
#text(23, 15.9,expression(paste("I3-GCS3 (12.77x",10^6, "mu-year 17)")),pos=1,cex = 0.8, srt = 0)
#points(17, 12.7,pch=23,bg="chartreuse4",lwd=1)

