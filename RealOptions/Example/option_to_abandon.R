#This example replicates the calculation in Chapter 7 of the book "Project Valuation using Real Option"
#example in page 102
library(fOptions)

S=100
X=95
Y=65 #value abandon
Time=5 #year to maturity

sigma=0.35
r=0.05
b=0.05
t=1 #interval of analysis
n=5 #number of steps

TypeFlag="pa"

optionvalue1=CRRBinomialTreeOption(TypeFlag = TypeFlag, S = S, X = Y, 
                      Time = Time, r = r, b = b, sigma = sigma, n = n)

print(optionvalue1)

optionvalue2=BinomialTreeOption(TypeFlag = TypeFlag, S=S, X=Y, 
                   Time=Time, r=r, b=b, sigma=sigma, n=n, title = NULL, description = NULL)

print(optionvalue2)

BinomialTreePlot(optionvalue2, dy = 1, cex = 0.8, ylim = c(-6, 7),
                 xlab = "n", ylab = "Option Value")
title(main = "Option Tree")

