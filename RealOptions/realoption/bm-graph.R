require("ape")
time_steps = 30
sigma=0.2
mu=0.0003
year=30
h0=180

h<-matrix(double(1),nrow=1,ncol=year)
h=h0
for (t in 1:year){
if (t==1){
h[t]=h0
}
else {
h[t]=h[t-1]*exp(mu*t) #+rnorm()
}
}
print(h)


sim.ou = function() { 
x.ou = numeric(time_steps)
for (i in 2:time_steps)
x.ou[i] = h + rnorm(h,0,sigma)
x.ou   # returns the value of x
}



X.ou = replicate(10, sim.ou())

layout(matrix(1:2, 2, 1))
yl = range(h)
y2 = range(h)
par(new=TRUE)
matplot(h, ylim = yl, type = "l", col = 1, main = "Brownian")
