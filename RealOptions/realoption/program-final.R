#This is real option program used for determining the optimal timing to carry out intervention
require(fOptions) # pls install package fOptions before running.
#function for estimate the call price
call.price <- function(x = 1, t = 0, T = 1, r = 1, sigma = 1, K=1){
d2 <- (log(x/K) + (r - 0.5 * sigma^2) * (T - t))/(sigma *sqrt(T - t))
d1 <- d2 + sigma * sqrt(T - t)
x * pnorm(d1) - K * exp(-r * (T - t)) * pnorm(d2)
}

u1=125
u2=120
v1=3.2
v2=2.6
coeff01_u=1.02
coeff02_u=1.0
coeff03_u=0.92
coeff04_u=1.05
coeff05_u=1.03
coeff06_u=0.99
coeff07_u=1.1
coeff08_u=0.95
coeff09_u=0.97
coeff10_u=1.07


coeff01_v=1.02
coeff02_v=1.0
coeff03_v=0.92
coeff04_v=1.05
coeff05_v=1.03
coeff06_v=0.99
coeff07_v=1.1
coeff08_v=0.95
coeff09_v=0.97
coeff10_v=1.07



#define the reliabiliy parameter of each section
lambda=c(u1*coeff01_u,u1*coeff02_u,u1*coeff03_u,u1*coeff04_u,u1*coeff05_u,u1*coeff06_u,u1*coeff07_u,u1*coeff08_u,u1*coeff09_u,u1*coeff10_u)
k=c(v1*coeff01_v,v1*coeff02_v,v1*coeff03_v,v1*coeff04_v,v1*coeff05_v,v1*coeff06_v,v1*coeff07_v,v1*coeff08_v,v1*coeff09_v,v1*coeff10_v)

lambda2=c(u2*coeff01_u,u2*coeff02_u,u2*coeff03_u,u2*coeff04_u,u2*coeff05_u,u2*coeff06_u,u2*coeff07_u,u2*coeff08_u,u2*coeff09_u,u2*coeff10_u)
k2=c(v2*coeff01_v,v2*coeff02_v,v2*coeff03_v,v2*coeff04_v,v2*coeff05_v,v2*coeff06_v,v2*coeff07_v,v2*coeff08_v,v2*coeff09_v,v2*coeff10_v)


O01=500000 #annual operational cost

O02=500000 #annual operational cost
N0=300000 #annual goods
C2<-4000000#intervention cost
C1<-0 #660000 #intervention cost
h0=15 # CHF
mu=0.0003
#defining the total time to estimate
TIME=80 #years
W=500000

B1<-matrix(double(1),nrow=1,ncol=TIME)#revenue
B11<-matrix(double(1),nrow=1,ncol=TIME)#revenue
B2<-matrix(double(1),nrow=1,ncol=TIME)#revenue
B22<-matrix(double(1),nrow=1,ncol=TIME+54)#revenue
B222<-matrix(double(1),nrow=1,ncol=TIME)#revenue
C<-matrix(double(1),nrow=1,ncol=TIME)#revenue
CC<-matrix(double(1),nrow=1,ncol=TIME)#revenue
N<-matrix(double(1),nrow=1,ncol=TIME)#number of annual visitor 
R1<-matrix(double(1),nrow=1,ncol=TIME)#routine maintainance
R2<-matrix(double(1),nrow=1,ncol=TIME)#routine maintainance
O1<-matrix(double(1),nrow=1,ncol=TIME)#Annual Operational cost
O2<-matrix(double(1),nrow=1,ncol=TIME)#Annual Operational cost


A_rel01<-matrix(double(1),nrow=1,ncol=TIME)#reliability of section 1 in year t - 
A_rel02<-matrix(double(1),nrow=1,ncol=TIME)#reliability of section 2 in year t - 
A_rel03<-matrix(double(1),nrow=1,ncol=TIME)#reliability of section 3 in year t - 
A_rel04<-matrix(double(1),nrow=1,ncol=TIME)#reliability of section 4 in year t - 
A_rel05<-matrix(double(1),nrow=1,ncol=TIME)#reliability of section 5 in year t - 
A_rel06<-matrix(double(1),nrow=1,ncol=TIME)#reliability of section 6 in year t - 
A_rel07<-matrix(double(1),nrow=1,ncol=TIME)#reliability of section 7 in year t - 
A_rel08<-matrix(double(1),nrow=1,ncol=TIME)#reliability of section 8 in year t - 
A_rel09<-matrix(double(1),nrow=1,ncol=TIME)#reliability of section 9 in year t - 
A_rel10<-matrix(double(1),nrow=1,ncol=TIME)#reliability of section 10 in year t - 

B_rel01<-matrix(double(1),nrow=1,ncol=TIME)#reliability of section 1 in year t - 
B_rel02<-matrix(double(1),nrow=1,ncol=TIME)#reliability of section 2 in year t - 
B_rel03<-matrix(double(1),nrow=1,ncol=TIME)#reliability of section 3 in year t - 
B_rel04<-matrix(double(1),nrow=1,ncol=TIME)#reliability of section 4 in year t - 
B_rel05<-matrix(double(1),nrow=1,ncol=TIME)#reliability of section 5 in year t - 
B_rel06<-matrix(double(1),nrow=1,ncol=TIME)#reliability of section 6 in year t - 
B_rel07<-matrix(double(1),nrow=1,ncol=TIME)#reliability of section 7 in year t - 
B_rel08<-matrix(double(1),nrow=1,ncol=TIME)#reliability of section 8 in year t - 
B_rel09<-matrix(double(1),nrow=1,ncol=TIME)#reliability of section 9 in year t - 
B_rel10<-matrix(double(1),nrow=1,ncol=TIME)#reliability of section 10 in year t - 




RE1<-matrix(double(1),nrow=1,ncol=TIME)#reliability in year t
RE2<-matrix(double(1),nrow=1,ncol=TIME)#reliability in year t
FA1<-matrix(double(1),nrow=1,ncol=TIME)#failure probability in year t
FA2<-matrix(double(1),nrow=1,ncol=TIME)#failure probability in year t
S<-matrix(double(1),nrow=1,ncol=TIME)#benefit over time
K<-matrix(double(1),nrow=1,ncol=TIME)#strike price (strike benefit)
h<-matrix(double(1),nrow=1,ncol=TIME)# fare price
Benefit<-matrix(double(1),nrow=1,ncol=TIME)# fare price



for (T in 1:TIME){
  
  
  A_rel01[T]<-exp(-(T/lambda[1])^k[1])
  A_rel02[T]<-exp(-(T/lambda[2])^k[2])
  A_rel03[T]<-exp(-(T/lambda[3])^k[3])
  A_rel04[T]<-exp(-(T/lambda[4])^k[4])
  A_rel05[T]<-exp(-(T/lambda[5])^k[5])
  A_rel06[T]<-exp(-(T/lambda[6])^k[6])
  A_rel07[T]<-exp(-(T/lambda[7])^k[7])
  A_rel08[T]<-exp(-(T/lambda[8])^k[8])
  A_rel09[T]<-exp(-(T/lambda[9])^k[9])
  A_rel10[T]<-exp(-(T/lambda[10])^k[10])
  
  
  
  B_rel01[T]<-exp(-(T/lambda2[1])^k2[1])
  B_rel02[T]<-exp(-(T/lambda2[2])^k2[2])
  B_rel03[T]<-exp(-(T/lambda2[3])^k2[3])
  B_rel04[T]<-exp(-(T/lambda2[4])^k2[4])
  B_rel05[T]<-exp(-(T/lambda2[5])^k2[5])
  B_rel06[T]<-exp(-(T/lambda2[6])^k2[6])
  B_rel07[T]<-exp(-(T/lambda2[7])^k2[7])
  B_rel08[T]<-exp(-(T/lambda2[8])^k2[8])
  B_rel09[T]<-exp(-(T/lambda2[9])^k2[9])
  B_rel10[T]<-exp(-(T/lambda2[10])^k2[10])
  
  
  
    RE1[T]<-A_rel01[T]*A_rel02[T]*A_rel03[T]*A_rel04[T]*A_rel05[T]*A_rel06[T]*A_rel07[T]*A_rel08[T]*A_rel09[T]*A_rel10[T]
    
  
  RE2[T]<-B_rel01[T]*B_rel02[T]*B_rel03[T]*B_rel04[T]*B_rel05[T]*B_rel06[T]*B_rel07[T]*B_rel08[T]*B_rel09[T]*B_rel10[T]

  
    
  if (T==1){
    N[T]=N0
    
    O1[T]=O01
    O2[T]=O02
    h[T]=h0
    B1[T]=(RE1[T]*(N[T]*h[T]-O1[T])-(1-RE1[T])*W)#*exp(-r*(T-1))
    B2[T]=(RE2[T]*(N[T]*h[T]-O2[T])-(1-RE2[T])*W)#*exp(-r*(T-1))
    S[T]<-B2[T]-C2
    K[T]<-B1[T]
  }  else {
    N[T]=N[T-1]+N[T-1]*0.5/100
    
    O1[T]=O01
    O2[T]=O02
    h[T]=h[T-1]*exp(mu*T)
    B1[T]=(RE1[T]*(N[T]*h[T]-O1[T])-(1-RE1[T])*W)*exp(-r*(T-1))
    B2[T]=(RE2[T]*(N[T]*h[T]-O2[T])-(1-RE2[T])*W)*exp(-r*(T-1))
    S[T]<-B2[T]-C2
    K[T]<-B1[T]-C1
    }
  S0 <- S[T]
  K0 <- K[T]
  r <- 0.02
 # t=T-1
  T <- T
  sigma <- 0.25
  Benefit[T] <- call.price(x = S0, t = 0, T = T, r = r, K = K0,sigma = sigma)
  }

color1=c("gray87","gray80","gray74","gray68","gray64","gray54","gray49","gray38","gray23","gray8")

color2=c("skyblue","skyblue2","skyblue3","skyblue4","slateblue","springgreen","springgreen3","steelblue","steelblue2","tan")


plot.new()
par(mar=c(5, 4, 4, 6) + 0.1)
plot(c(1:TIME), A_rel01, main="", pch=1,type="l",xlab=" ", ylab="", lwd=1,col=color1[1],axes=FALSE,lty=11)
axis(2, ylim=c(0,1),col="black",las=1)
axis(1, xlim=c(0,TIME),col="black",las=1)


mtext(expression(paste("Time units")),side=1,col="black",line=2.2) 
mtext(expression(paste("Reliability")),side=2,col="black",line=3) 


grid(10, 10, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = TRUE)
box()

#stop("debugg")


par(new=TRUE)
plot(c(1:TIME), A_rel02, main="", pch=2,type="l",xlab=" ", ylab="", lwd=1,col=color1[2],axes=FALSE,lty=2)

par(new=TRUE)
plot(c(1:TIME), A_rel03, main="", pch=3,type="l",xlab=" ", ylab="", lwd=1,col=color1[3],axes=FALSE,lty=3)

par(new=TRUE)
plot(c(1:TIME), A_rel04, main="", pch=4,type="l",xlab=" ", ylab="", lwd=1,col=color1[4],axes=FALSE,lty=4)

par(new=TRUE)
plot(c(1:TIME), A_rel05, main="", pch=5,type="l",xlab=" ", ylab="", lwd=1,col=color1[5],axes=FALSE,lty=5)

par(new=TRUE)
plot(c(1:TIME), A_rel06, main="", pch=6,type="l",xlab=" ", ylab="", lwd=1,col=color1[6],axes=FALSE,lty=6)

par(new=TRUE)
plot(c(1:TIME), A_rel07, main="", pch=7,type="l",xlab=" ", ylab="", lwd=1,col=color1[7],axes=FALSE,lty=7)

par(new=TRUE)
plot(c(1:TIME), A_rel08, main="", pch=8,type="l",xlab=" ", ylab="", lwd=1,col=color1[8],axes=FALSE,lty=8)

par(new=TRUE)
plot(c(1:TIME), A_rel09, main="", pch=9,type="l",xlab=" ", ylab="", lwd=1,col=color1[9],axes=FALSE,lty=9)

par(new=TRUE)
plot(c(1:TIME), A_rel10, main="", pch=10,type="l",xlab=" ", ylab="", lwd=1,col=color1[10],axes=FALSE,lty=10)


par(new=TRUE)
plot(c(1:TIME), B_rel01, main="", pch=1,type="l",xlab=" ", ylab="", lwd=1,col=color2[1],axes=FALSE)

par(new=TRUE)
plot(c(1:TIME), B_rel02, main="", pch=2,type="l",xlab=" ", ylab="", lwd=1,col=color2[2],axes=FALSE)

par(new=TRUE)
plot(c(1:TIME), B_rel03, main="", pch=3,type="l",xlab=" ", ylab="", lwd=1,col=color2[3],axes=FALSE)

par(new=TRUE)
plot(c(1:TIME), B_rel04, main="", pch=4,type="l",xlab=" ", ylab="", lwd=1,col=color2[4],axes=FALSE)

par(new=TRUE)
plot(c(1:TIME), B_rel05, main="", pch=5,type="l",xlab=" ", ylab="", lwd=1,col=color2[5],axes=FALSE)

par(new=TRUE)
plot(c(1:TIME), B_rel06, main="", pch=6,type="l",xlab=" ", ylab="", lwd=1,col=color2[6],axes=FALSE)

par(new=TRUE)
plot(c(1:TIME), B_rel07, main="", pch=7,type="l",xlab=" ", ylab="", lwd=1,col=color2[7],axes=FALSE)

par(new=TRUE)
plot(c(1:TIME), B_rel08, main="", pch=8,type="l",xlab=" ", ylab="", lwd=1,col=color2[8],axes=FALSE)

par(new=TRUE)
plot(c(1:TIME), B_rel09, main="", pch=9,type="l",xlab=" ", ylab="", lwd=1,col=color2[9],axes=FALSE)

par(new=TRUE)
plot(c(1:TIME), B_rel10, main="", pch=10,type="l",xlab=" ", ylab="", lwd=1,col=color2[10],axes=FALSE)

par(new=TRUE)
plot(c(1:TIME), RE1, main="", pch=16,type="l",xlab=" ", ylab="", lwd=2,col="black",axes=FALSE)

par(new=TRUE)
plot(c(1:TIME), RE2, main="", pch=16,type="l",xlab=" ", ylab="", lwd=2,col="red",axes=FALSE)

legend("topright",inset=.02,legend=c("Before intervention","After intervention"),col=c("black","red"),lty=c(1,1),lwd=c(2,2),horiz=F)


cat("value of option")
print(Benefit)

tuu=46
for (T in 1:TIME){
  
  if (T==1){
    C[T]=sum(B2)
  }  else  {
    C[T]=sum(B1[c(1:T)])+sum(B2[c(1:(TIME-T))])-C2*exp(-r*(T-1))
    CC[T]=sum(B1[c(1:T)])+sum(B1[c(1:(TIME-T))])-C1*exp(-r*(T-1))
     B22[T+tuu]=(RE2[T]*(N[T+tuu]*h[T+tuu]-O2[T])-(1-RE2[T])*W)*exp(-r*(T+0))
    B222[T]=B22[T]
  }
  
}

plot.new()
trucy=1000000
plot(c(1:TIME), (C)/trucy/100, main="", pch=6,cex=0.6,type="b",xlab="", ylab="", lwd=2,col="blue",axes=FALSE,ylim=c(0,9))
library(MASS)
axis(2, ylim=c(0,max(B222)),col="darkblue",las=1)
mtext(expression(paste("Annual profit  (", 10^8," mu)")),side=2,line=2.2, adj = 0.5 )
axis(1,pretty(c(0:80),10))

abline(a = sum(B1)/trucy/100, b = 0,lty=4,lwd=2, pch=13,cols="black")

#lines(0,sum(B1),80,sum(B1))
#par(new=TRUE)
#plot(c(1:TIME), B2/trucy, main="", pch=10,cex=0.5,type="b",xlab="", ylab="", lwd=1,col="coral4",axes=FALSE,ylim=c(0,1*max(B222)/trucy))


#par(new=TRUE)
#plot(c(1:TIME), C/trucy/100, main="", pch=1,type="b",xlab="", ylab="", lwd=1,col="blue",axes=FALSE,ylim=c(0,1.2*max(B1)/trucy))

par(new=TRUE)
plot(c(1:TIME), (C-sum(B1))/trucy/100, main="", pch=16,cex=0.8,type="b",xlab="", ylab="", lwd=1,col="red",axes=FALSE,ylim=c(0.5*min((C-sum(B1)))/trucy/100,1.2*max((C-sum(B1)))/trucy/100))
axis(4, ylim=c(min((C-sum(B1)))/trucy/100,max((C-sum(B1)))/trucy/100),col="darkblue",las=1)
#axis(4,c(seq(min((C-sum(B1))),max((C-sum(B1))),by=1000)),c(seq(min((C-sum(B1)))/trucy/100,max((C-sum(B1)))/trucy/100,by=1)))

mtext(expression(paste("Difference  (", 10^8," mu)")),side=4,line=2.5, adj = 0.5 )
mtext(expression(paste("Time units")),side=1,col="black",line=2)  
box()


#par(new=TRUE)
#plot(c(1:TIME), C/trucy/100, main="", pch=6,cex=0.6,type="b",xlab="", ylab="", lwd=1,col="blue",axes=FALSE,ylim=c(0,5))


#grid(10, 10, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = TRUE)
legend("topleft",inset=.02,legend=c("Do-nothing", "Execute intervention", "Difference"),lty=c(2,2,3),pch=c(-1,6,1),col=c("black","blue","red"),horiz=F)

print(sum(S))

write.matrix((C-sum(B1)), file="benefit.csv",sep=",")


print(which.max((C-sum(B1))))
print(max((C-sum(B1))))


