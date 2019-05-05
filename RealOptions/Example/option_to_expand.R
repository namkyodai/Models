#This example replicates the calculation in Chapter 7 of the book "Project Valuation using Real Option"
#example in page 111
library(fOptions)

S=80
X=200

Time=4 #year to maturity

sigma=0.3
r=0.05
b=0.05
t=1 #interval of analysis
n=4 #number of steps
E=3 #expansion factor




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

