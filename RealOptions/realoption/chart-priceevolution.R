require(sde) #reading stochastic differential equation package
gbm=GBM(180,0.0003,0.2,30,1000)
#plot(gbm)
sigma=0.2
mu=0.0003
year=30
h0=180 #CHF the price at the begining

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
plot(h,pch=13,axes=FALSE,ylim=c(170,210),ylab="", xlab="",col="red",type="o",lwd=2)
h1=h+rnorm(h,0,sigma)
par(new=TRUE)
print(h)
plot(h1,pch=13,axes=FALSE,ylim=c(170,210),ylab="", xlab="",col="blue1",type="o",lwd=1)
axis(2, ylim=c(0,1),col="darkblue",las=1)  ## las=1 makes horizontal labels
mtext(expression(paste("Price (mu)")),side=3,line=0.2, adj = -0.06 )
axis(1,pretty(range(x),10))
mtext(expression(paste("Time (years)")),side=1,col="black",line=2)  
box()
par(new=TRUE)
h2=h+rnorm(h,0,sigma)
plot(h2,pch=13,axes=FALSE,ylim=c(170,210),ylab="", xlab="",col="cyan3",type="o",lwd=1)



legend("topleft",inset=.05,legend=c("Fixed increase of price","Stochastic variation 1","Stochastic variation 2"),text.col=c("red","blue1","cyan3"),pch=c(16,15,14),col=c("red","blue1","cyan3"),horiz=F)


