MCPrice <- function(x = 1, t = 0, T = 1, r = 1, sigma = 1,M = 10000,f){
h <- function(m) {
u <- rnorm(m) #random generating value of u following normal distribution with mean 0 and standard deviation 1
tmp <- c(x * exp((r - 0.5 * sigma^2) * (T - t) + sigma *sqrt(T - t) * u), x * exp((r - 0.5 * sigma^2) * (T -t) + sigma * sqrt(T - t) * (-u)))
mean(sapply(tmp, function(xx) f(xx)))
}
p <- h(M)
p * exp(-r * (T - t))
}
S0 <- 100
K <- 110
r <- 0.02
T <- 1/4
sigma <- 0.25
library(fOptions)
price=GBSOption(TypeFlag = "c", S = S0, X = K, Time = T, r = r,b = r,sigma = sigma)
print(price)
#start the MCMC 

f <- function(x) max(0, x - K)
set.seed(123)
M <- 50000 #number of simulation.
price=MCPrice(x = S0, t = 0, T = T, r = r, sigma, M = M, f = f)
cat("with MCMC, the price of the call option is ")
print(price)
