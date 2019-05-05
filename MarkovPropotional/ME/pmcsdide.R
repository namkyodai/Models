require(graphics)

#print Figure pmcsdide
postscript(file="pmcsdide.eps", onefile=FALSE, horizontal=FALSE,
           width=7, height = 4, paper="a4", family="Times")
# Trim off excess margin space (bottom, left, top, right)
par(mar=c(3, 5, 0.2, 0.7))
# Trim off excess outer margin space (bottom, left, top, right)
par(oma=c(0,0,0,0))
# Trim off excess space for label and ticks (label, ticks, line)
par(mgp=c(1.9,0.6,0))
# lty: line styles (1=solid, 2=dash, 3=dot, 4=dash-dot)
# lab: (# of x-ticks, # of y-ticks, len of ticks), approximately
# lwd: line-width
# cex.lab: fontsize scaling-factor for labels

source("Hazard2Markov.R")
dev.off()





