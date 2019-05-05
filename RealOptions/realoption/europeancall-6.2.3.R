call.price <- function(x = 1, t = 0, T = 1, r = 1, sigma = 1, K=1){
d2 <- (log(x/K) + (r - 0.5 * sigma^2) * (T - t))/(sigma *sqrt(T - t))
d1 <- d2 + sigma * sqrt(T - t)
x * pnorm(d1) - K * exp(-r * (T - t)) * pnorm(d2)
}
