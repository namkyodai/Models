# Evolution
data=read.csv("profit_final.csv", header=T)
attach(data)
#object 1
scalefactor=1000000
trucy=100
colors=c("red", "darkorchid","green", "violet", "darkorange3","brown3")

plot(data$DS2/scalefactor,pch=16,axes=FALSE,ylim=c(0,trucy),ylab="", xlab="",col="red",type="o",lwd=1)

axis(2, ylim=c(0,trucy),col="darkblue",las=1)
mtext(expression(paste("intervention window value (", 10^6," mu)")),side=3,line=0.2, adj = -0.06 )

axis(1,pretty(range(data$Time),15))

mtext(expression(paste("Time units")),side=1,col="black",line=2)  
box()

grid(10, 10, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = TRUE)
#legend("topright",inset=.02,legend=c("Renewal"),text.col=colors,pch=c(16,15),col=colors,horiz=F)

#arrows(which.max(data$DS2), max(data$DS2)/scalefactor*1.05, which.max(data$DS2), max(data$DS2)*1.2/scalefactor,length = 0.15, angle = 20,code = 1,col= 'blue3',lty = 1,lwd=2)
text(which.max(data$DS2), max(data$DS2)/scalefactor*1.2,expression(paste("(69.42x",10^6, "mus-at 46 tus)")),pos=1,cex = 1, srt = 0)
text(which.max(data$DS2)*1.061, max(data$DS2)/scalefactor*1.1,"A",pos=1,cex = 1, srt = 0)

points(which.max(data$DS2),max(data$DS2)/scalefactor,pch=23,bg="black",lwd=3)

segments(which.max(data$DS2), 0, which.max(data$DS2), max(data$DS2)/scalefactor, col= 'cadetblue',lty=2,lwd=0.5)
segments(0, max(data$DS2)/scalefactor, which.max(data$DS2), max(data$DS2)/scalefactor, col= 'cadetblue',lty=2,lwd=0.5)

X=65
text(X*1.051, data$DS2[X]/scalefactor*1.12,"B",pos=1,cex = 1, srt = 0)
points(X,data$DS2[X]/scalefactor,pch=23,bg="black",lwd=3)
segments(X, 0, X, data$DS2[X]/scalefactor, col= 'cadetblue',lty=2,lwd=0.5)
segments(0, data$DS2[X]/scalefactor, X, data$DS2[X]/scalefactor, col= 'cadetblue',lty=2,lwd=0.5)

Y=20
text(Y, data$DS2[Y]/scalefactor*1.4,"C",pos=1,cex = 1, srt = 0)
points(Y,data$DS2[Y]/scalefactor,pch=23,bg="black",lwd=3)
segments(Y, 0, Y, data$DS2[Y]/scalefactor, col= 'cadetblue',lty=2,lwd=0.5)
segments(0, data$DS2[Y]/scalefactor, Y, data$DS2[Y]/scalefactor, col= 'cadetblue',lty=2,lwd=0.5)


