#Reading data for running the model
#Define the parameters of the model with dimensions specified.

N0=200000 #annual visitor at time 0
h0=180 # CHF
mu=0.0003 #the growth rate of Geometric brownian motion process (Winear proces)
r <- 0.02 #interest rate
sigma <- 0.2 #standard deviation of normal distribution
scalefactor=1000000

#for reference strategy
alpha1=0.002
m1=2.5
C1<-0*scalefactor 		# intervention cost
R01=0.2*scalefactor	#annual routine maintenance
O01=0.72*scalefactor 	#annual operational cost
W1=0.4*scalefactor

#for investigated strategy
alpha2=0.002
m2=2.5
C2<-1.82*scalefactor #intervention cost
R02=0.213*scalefactor 	#annual routine maintenance
O02=0.72*scalefactor 	#annual operational cost
W2=0.4*scalefactor