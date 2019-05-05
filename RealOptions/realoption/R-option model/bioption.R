bioption=function (TypeFlag = c("ce", "pe", "ca", "pa"), S, X, Time, r, 
    b, sigma, f, n, title = NULL, description = NULL) 
{
    TypeFlag = TypeFlag[1]
    if (TypeFlag == "ce" || TypeFlag == "ca") 
        z = +1
    if (TypeFlag == "pe" || TypeFlag == "pa") 
        z = -1
    dt = Time/n
    u = exp(sigma * sqrt(dt))
    d = 1/u
    p = (exp(b * dt) - d)/(u - d)
    Df = exp(-r * dt)
    OptionValue = z * (S * u^(0:n) * d^(n:0) - X)
    offset = 1
    Tree = OptionValue = (abs(OptionValue) + OptionValue)/2
    if (TypeFlag == "ce" || TypeFlag == "pe") {
        for (j in (n - 1):0) {
            Tree <- c(Tree, rep(0, times = n - j))
            for (i in 0:j) {
                OptionValue[i + offset] = (p * OptionValue[i + 
                  1 + offset] + (1 - p) * OptionValue[i + offset]) * 
                  Df
                Tree = c(Tree, OptionValue[i + offset])
            }
        }
    }
    if (TypeFlag == "ca" || TypeFlag == "pa") {
        for (j in (n - 1):0) {
            Tree <- c(Tree, rep(0, times = n - j))
            for (i in 0:j) {
                OptionValue[i + offset] = max((z * (S * u^i * d^(abs(i - j)) - X)), (p * OptionValue[i + 1 + offset] + (1 - p) * OptionValue[i + offset]) * 
                  Df)
                Tree = c(Tree, OptionValue[i + offset])
            }
        }
    }
    Tree = matrix(rev(Tree), byrow = FALSE, ncol = n + 1)
    invisible(Tree)
}