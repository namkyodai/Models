
#Creates binomial tree with possible future values of x0
require(fOptions)
years=9
S=0.9466
K=0.95
r=0.05
sigma=0.2349
Time=45

#years=2
#S=50
#K=52
#r=0.05
#sigma=0.2231
#Time=2


source("up-down.R") #this is compute the value of price go up and down along the tree.
source("bioption.R")
nodevalue=upanddown(TypeFlag = "ce", S = S, X = K, Time = Time, r = r, b =r, sigma = sigma, n = years)

option = bioption(TypeFlag = "ce", S = S, X = K, Time = Time, r = r, b =r, sigma = sigma, n = years)


#plot the evolution of price
BinomialTreePlot(nodevalue, dy = 1, cex = 0.8, ylim = c(-15, 15),xlab = "n", ylab = "Option Value")
title(main = "Option Tree")

# plot the value of option at each node
BinomialTreePlot(option, dy = 1, cex = 0.8, ylim = c(-11, 11),xlab = "Time", ylab = "Option Value")
title(main = "Option Tree")
