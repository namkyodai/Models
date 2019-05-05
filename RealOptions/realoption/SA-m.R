data=read.csv("SA-m.csv", header=T)
scalefactor=1000000
## add extra space to right margin of plot within frame
par(mar=c(5, 4, 4, 6) + 0.1)

## Plot first set of data and draw its axis
plot(data$m, data$Extra, pch=2, axes=FALSE, ylim=c(5,30), xlab="", ylab="", 
     type="b",col="blueviolet", lwd=2)
axis(2, ylim=c(0,2),col="black",las=1)  ## las=1 makes horizontal labels
mtext(expression(paste("Option value  (", 10^6, "mu)")),side=2,col="black",line=2.5)
box()

## Allow a second plot on the same graph
par(new=TRUE)

## Plot the second plot and put axis scale on right
plot(data$m, data$Time, pch=4,  xlab="", ylab="", ylim=c(7,21), 
     axes=FALSE, type="b", col="coral1",lwd=2,lty=2)
## a little farther out (line=4) to make room for labels
mtext("OTD (years)",side=4,col="black",line=2) 
axis(4, ylim=c(0,20), col="black",col.axis="black",las=1)

## Draw the time axis
axis(1,pretty(range(data$m),10))
mtext(expression(paste("Deterioration parameter of DS2  ", m)),side=1,col="black",line=2.5)  
grid(10, 10, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = TRUE)
## Add Legend
legend("topright",inset=.04,legend=c("Option value (mu)","Time (Year)"),
       text.col=c("blueviolet","coral1"),pch=c(2,4),col=c("blueviolet","coral1"),horiz=F,cex=0.9)