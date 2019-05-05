require(graphics)
totalyear=100
jmax=5
csstate1<-matrix(double(1),nrow=totalyear,ncol=jmax) #state probability in
csstate2<-matrix(double(1),nrow=totalyear,ncol=jmax) #state probability in
dif1<-matrix(double(1),nrow=totalyear,ncol=jmax) #s
dif2<-matrix(double(1),nrow=totalyear,ncol=1) #s

source("MarkovPo.R")
csstate1=csstate/100
source("Hazard2Markov.R")
csstate2=csstate/100

dif1=csstate2-csstate1
dif2=dif1[,1]+dif1[,2]+dif1[,3]+dif1[,4]+dif1[,5]


#print Figure pmcsdide
postscript(file="evaluation.eps", onefile=FALSE, horizontal=FALSE,
           width=7, height = 5, paper="a4", family="Times")
# Trim off excess margin space (bottom, left, top, right)
par(mar=c(3, 3, 0.2, 0.7))
# Trim off excess outer margin space (bottom, left, top, right)
par(oma=c(0,0,0,0))
# Trim off excess space for label and ticks (label, ticks, line)
par(mgp=c(1.9,0.6,0))
# lty: line styles (1=solid, 2=dash, 3=dot, 4=dash-dot)
# lab: (# of x-ticks, # of y-ticks, len of ticks), approximately
# lwd: line-width
# cex.lab: fontsize scaling-factor for labels

colors=c("black","darkorchid","red")
# 4 figures arranged in 2 rows and 2 columns
par(mfrow=c(3,3))
plot(dif1[,1],pch=1, ylab="error term",xlab="Time (years)",col="darkorchid",type="b",lwd=0.1)
abline(h = 0, col = "red",lwd=0.8)
grid(10, 10, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)
legend("bottomright",inset=.02,legend=c("CS1", "error term","0 line"),text.col=colors,pch=c(NA,1,NA),lty=c(0,0,1),lwd=c(0,0.1,0.1),col=colors,horiz=F,bty != "n",bg="white",cex=0.8)
#text(80, 0.1,expression(paste(sum(epsilon[i],i=1,5)^{CS1},"=")),pos=1,cex = 1, srt = 0)


plot(dif1[,2],pch=1, ylab="error term",xlab="Time (years)",col="darkorchid",type="b",lwd=0.1)
abline(h = 0, col = "red",lwd=0.8)
grid(10, 10, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)
legend("bottomright",inset=.02,legend=c("CS2", "error term","0 line"),text.col=colors,pch=c(NA,1,NA),lty=c(0,0,1),lwd=c(0,0.1,0.1),col=colors,horiz=F,bty != "n",bg="white",cex=0.8)
#text(0, 0,expression(paste("CS1")),pos=1,cex = 0.8, srt = 0)

plot(dif1[,3],pch=1, ylab="error term",xlab="Time (years)",col="darkorchid",type="b",lwd=0.1)
abline(h = 0, col = "red",lwd=0.8)
grid(10, 10, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = TRUE)
legend("bottomright",inset=.02,legend=c("CS3", "error term","0 line"),text.col=colors,pch=c(NA,1,NA),lty=c(0,0,1),lwd=c(0,0.1,0.1),col=colors,horiz=F,bty != "n",bg="white",cex=0.8)
#text(0, 0,expression(paste("CS1")),pos=1,cex = 0.8, srt = 0)

plot(dif1[,4],pch=1, ylab="error term",xlab="Time (years)",col="darkorchid",type="b",lwd=0.1)
abline(h = 0, col = "red",lwd=0.8)
grid(10, 10, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = TRUE)
legend("bottomright",inset=.02,legend=c("CS4", "error term","0 line"),text.col=colors,pch=c(NA,1,NA),lty=c(0,0,1),lwd=c(0,0.1,0.1),col=colors,horiz=F,bty != "n",bg="white",cex=0.8)
#text(0, 0,expression(paste("CS1")),pos=1,cex = 0.8, srt = 0)

plot(dif1[,5],pch=1, ylab="error term",xlab="Time (years)",col="darkorchid",type="b",lwd=0.1)
abline(h = 0, col = "red",lwd=0.8)

grid(10, 10, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)
legend("bottomright",inset=.02,legend=c("CS5", "error term","0 line"),text.col=colors,pch=c(NA,1,NA),lty=c(0,0,1),lwd=c(0,0.1,0.1),col=colors,horiz=F,bty != "n",bg="white",cex=0.8)
#text(0, 0,expression(paste("CS1")),pos=1,cex = 0.8, srt = 0)

plot(dif2,pch=1, ylab="error term",xlab="Time (years)",col="darkorchid",type="b",lwd=0.1)
abline(h = 0, col = "red",lwd=0.8)
grid(10, 10, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)
legend("bottomright",inset=.02,legend=c("Total", "error term","0 line"),text.col=colors,pch=c(NA,1,NA),lty=c(0,0,1),lwd=c(0,0.1,0.1),col=colors,horiz=F,bty != "n",bg="white",cex=0.8)
#text(0, 0,expression(paste("CS1")),pos=1,cex = 0.8, srt = 0)

dev.off()





