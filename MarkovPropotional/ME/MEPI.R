#Purpose: generate proportional data from a ME model for a bridge element
# by Nam Lethanh (namkyodai@gmail.com)
#date: 30 November, 2015
#------------------------------------------------------
rm(list=ls()) #: clear the memory and objects in R console
#------------------------------------------------------
#.....read inputparameters from file
input <- read.csv("input.csv", header =TRUE, sep = ",") #: import data from excel
data <- data.frame(input) #: making data to be in frame
Omax <- length(data[,1])/2 #: number of objects
#----------------------------------------------------
#.....STEP 1: Define model parameter
T <- 100 #: years of investigation
limit_chrolide <- 0.48 #: limit concentration of chloride [wt.-%Cl/binder]
limit_crack <- 0.50 #: limit crack width [mm]
N <- 20000 #: maximum number used for random sampling with Monte Carlo or sampling techniques [-]
ftspl <- 2.6 #: tensil splitting strength [MPa]
epsilon <- 10e-16 #: this is used to prevent value of simulated standard deviation to be less than 0.005, which gives a distribution with value less than 10^(-16), and in R, if value is less than 10^(-16), the PC cannot understand, basically when that situation happens the result returns to 0, which is not true.
#.....Define the number of condition state
I <- array(dim=c(Omax,1)) #: number of CS related to chloride-induced onset of corrosion (inition phase), the last condition state related to failure state
for (o in 1:Omax){
  I[o] = data[(o-1)*2+1,1]
}
J <- array(dim=c(Omax,1)) #: number of CS related to propagation of corrosion (propagation phase "cracking")
for (o in 1:Omax){
  J[o] = data[(o-1)*2+2,1]
}
IM <- max(I)
a <- array(dim=c(Omax,IM)) #: range value of condition state for chloride [wt.-%Cl/binder]
for (o in 1:Omax){
  for (i in 1:IM){
    a[o,i] = 0
  }
}
for (o in 1:Omax){
  a[o,(1:I[o])] = seq(from=0, to=limit_chrolide,length.out=I[o])
}
JM <- max(J)
b <- array(dim=c(Omax,JM)) #: range value of condition state for crack [mm]
for (o in 1:Omax){
  for (i in 1:JM){
    b[o,i] = 0
  }
}
for (o in 1:Omax){
  b[o,(1:J[o])] = seq(from=0, to=limit_crack,length.out=J[o])
}
#.....Define parameter for initiation phase
d <- array(dim=c(Omax,1)) #: distance from the concrete surface to rebar [mm] (if the onset of corrosion is considered, d is equal to the concrete cover depth dc)
for (o in 1:Omax){
  d[o] = data[(o-1)*2+1,2]
}
D0 <- array(dim=c(Omax,2)) #: mean value and standard deviation of chloride migration coefficient at defined compaction, measured at time t0 [mm2/yr]
for (o in 1:Omax){
  for (v in 1:2){
    D0[o,v] = data[(o-1)*2+v,3]
  }
}
Cs <- array(dim=c(Omax,2)) #: mean value and standard deviation of the surface chloride level [wt.-%Cl/binder]
for (o in 1:Omax){
  for (v in 1:2){
    Cs[o,v] = data[(o-1)*2+v,4]
  }
}
ke <- array(dim=c(Omax,2)) #: shape and scale value of the parameter which considers the influence of environment on D0 [-]
for (o in 1:Omax){
  for (v in 1:2){
    ke[o,v] = data[(o-1)*2+v,5]
  }
}
kt <- array(dim=c(Omax,2)) #: mean value and standard deviation of the parameter which considers the influence of test method on D0 [-]
for (o in 1:Omax){
  for (v in 1:2){
    kt[o,v] = data[(o-1)*2+v,6]
  }
}
kc <- array(dim=c(Omax,1)) #: parameter which considers the influence of curing on D0 [-]
for (o in 1:Omax){
  kc[o] = data[(o-1)*2+1,7]
}
t0 <- array(dim=c(Omax,1)) #: reference periode [yr]
for (o in 1:Omax){
  t0[o] = data[(o-1)*2+1,8]
}
nn <- array(dim=c(Omax,2)) #: shape value 1 and 2 of the age factor [-]
for (o in 1:Omax){
  for (v in 1:2){
    nn[o,v] = data[(o-1)*2+v,9]
  }
}
#.....Define parameter for propagation phase
w0 <- array(dim=c(Omax,2)) #: mean value and standard deviation of the crack width when it is visible
for (o in 1:Omax){
  for (v in 1:2){
    w0[o,v] = data[(o-1)*2+v,10]
  }
}
beta <- array(dim=c(Omax,2)) #: mean value and standard deviation of the propagation controlling parameter [-]
for (o in 1:Omax){
  for (v in 1:2){
    beta[o,v] = data[(o-1)*2+v,11]
  }
}
V0 <- array(dim=c(Omax,2)) #: mean value and standard deviation of the corrosion rate when corrosion is active [mm/yr]
for (o in 1:Omax){
  for (v in 1:2){
    V0[o,v] = data[(o-1)*2+v,12]
  }
}
wet <- array(dim=c(Omax,2)) #: mean value and standard deviation of the wetness period (sum of all raining days in a year) [-]
for (o in 1:Omax){
  for (v in 1:2){
    wet[o,v] = data[(o-1)*2+v,13]
  }
}
alpha <- array(dim=c(Omax,2)) #: mean value and standard deviation of the pitting factor [-]
for (o in 1:Omax){
  for (v in 1:2){
    alpha[o,v] = data[(o-1)*2+v,14]
  }
}
a1 <- array(dim=c(Omax,2)) #: mean value and standard deviation of the regression parameter 1 [mm]
for (o in 1:Omax){
  for (v in 1:2){
    a1[o,v] = data[(o-1)*2+v,15]
  }
}
a2 <- array(dim=c(Omax,2)) #: mean value and standard deviation of the regression parameter 2 [mm]
for (o in 1:Omax){
  for (v in 1:2){
    a2[o,v] = data[(o-1)*2+v,16]
  }
}
a3 <- array(dim=c(Omax,2)) #: mean value and standard deviation of the regression parameter 3 [mm/MPa]
for (o in 1:Omax){
  for (v in 1:2){
    a3[o,v] = data[(o-1)*2+v,17]
  }
}
phi <- array(dim=c(Omax,1)) #: rebar diameter [mm]
for (o in 1:Omax){
  phi[o] = data[(o-1)*2+1,18]
}
#----------------------------------------------------
#.....STEP 2: Define functions to be used
#.....Define error function
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
#.....Model parameters 
#......Define the function of chloride concentration in initiation phase
Ccl <- function(Cs,x,ke,kt,kc,D0,t0,n,t){
  var=ke*kt*kc*D0*(t0/t)^n*t
  if (var >  epsilon) {
    Cs*(1-erf(x/(2*sqrt(var))))
  }
  else {
    Cs*(1-erf(x/(2*sqrt(epsilon))))
  }
}
#.....Cs: surface chloride level [wt.-%Cl/binder] (probabilistic, ND (LN))
#.....x: cover depth [mm]
#.....ke: influence of environment on D0 [-] (probabilistic, Gamma (ND))
#.....kt: influence of test method on D0 [-] (probabilistic, ND)
#.....kc: influence of curing on D0 [-] (D (probabilistic, Beta))
#.....D0: chloride migration coefficient [mm2/yr] (probabilistic, ND)
#.....t0: reference periode [yr]
#.....t: exposure periode [yr]
#.....n: age factor [-] (probabilistic, Beta (ND))
#----------------------------------------------------
#.....Define the dimension
pi_cl <- array(dim=c(T,IM,Omax)) #: state probability for initiation stage
for (o in 1:Omax){
  for (t in 1:T){
    for (i in 1:IM){
      pi_cl[t,i,o] = 0
    }
  }
}
pi_cr <- array(dim=c(T,JM,Omax)) #: state probability for propagation stage
for (o in 1:Omax){
  for (t in 1:T){
    for (j in 1:JM){
      pi_cr[t,j,o] = 0
    }
  }
}
pi <- array(dim=c(T,(IM+JM-1),Omax)) #: state probability of the entire process -1 mean we consider the last state I as failure state, therefore, the total state of the system is I+J-1
for (o in 1:Omax){
  for (t in 1:T){
    for (i in 1:(IM+JM-1)){
      pi[t,i,o] = 0
    }
  }
}
Cclt <- array(dim=c(T,N,Omax))
Ccrt <- array(dim=c(T,N,Omax))
Ccltd <- array(dim=c(T,Omax))
Ccrtd <- array(dim=c(T,Omax))
#----------------------------------------------------
#......Define the function for the crack width (propagation phase)
crackwidth <- function(w0,beta,V0,alpha,wet,a1,a2,a3,c,phi,ftspl,t){
  w0+beta*(V0*wet*alpha*t-(a1+a2*c/phi+a3*ftspl))
}
# .....w(x): crack width at time x [mm]
# .....w(0): crack width when it is visible 0.05 [mm] (probabilistic, ND)
# .....beta: propagation controlling parameter [-] (probabilistic, ND)
# .....V0: mean corrosion rate when corrosion is active [mm/yr] (probabilistic, ND (Weibull))
# .....V: corrosion rate [mm/yr]
# .....wet: wetness period [-] (sum of all raining days in a year) (probabilistic, ND)
# .....alpha: pitting factor taking non-uniform corrosion of the rebars into account [-] (probabilistic, ND)
# .....a1: regression parameter 1 [mm] (probabilistic, ND)
# .....a2: regression parameter 2 [mm] (probabilistic, ND)
# .....a3: regression parameter 3 [mm/MPa] (probabilistic, ND)
# .....c: cover depth [mm] (D)
# .....phi: rebar diameter [mm] (D)
# .....ftspl: tensile splitting strength [MPa] (D)
#----------------------------------------------------
#.....STEP 3: ESTIMATION
  #.....for parameter concerning the initiation phase
  i1 <- array(dim=c(N,Omax)) #: i1 represent for surface Cl concentration Cs
  for (o in 1:Omax){
    i1[,o] = rnorm(N,Cs[o,1],Cs[o,2]) 
  }
  i2 <- array(dim=c(N,Omax)) #: i2 represent for diffusion coefficient D0
  for (o in 1:Omax){
    i2[,o] = rnorm(N,D0[o,1],D0[o,2]) 
  }
  i3 <- array(dim=c(N,Omax)) #: i3 represent the environmental factor ke
  for (o in 1:Omax){
    i3[,o] = rgamma(N,ke[o,1],scale = ke[o,2])
  }
  i4 <- array(dim=c(N,Omax)) #: i4 represent the test method factor kt
  for (o in 1:Omax){
    i4[,o] = rnorm(N,kt[o,1],kt[o,2])
  }
  i5 <- array(dim=c(N,Omax)) #: i5 represent the age factor n
  for (o in 1:Omax){
    i5[,o] = rbeta(N,nn[o,1],nn[o,2],ncp = 0)
  }
  #.....for parameter concerning the propagation phase
  p1 <- array(dim=c(N,Omax)) #: p1 represent the crack width when it is visible w0
  for (o in 1:Omax){
    p1[,o] = rnorm(N,w0[o,1],w0[o,2])
  }
  p2 <- array(dim=c(N,Omax)) #: p2 represent the propagation controlling parameter beta
  for (o in 1:Omax){
    p2[,o] = rnorm(N,beta[o,1],beta[o,2])
  }
  p3 <- array(dim=c(N,Omax)) #: p3 represent the mean corrosion rate V0
  for (o in 1:Omax){
    p3[,o] = rnorm(N,V0[o,1],V0[o,2])
  }
  p4 <- array(dim=c(N,Omax)) #: p4 represent the pitting factor alpha
  for (o in 1:Omax){
    p4[,o] = rnorm(N,alpha[o,1],alpha[o,2])
  }
  p5 <- array(dim=c(N,Omax)) #: p5 represent the wetness period wet
  for (o in 1:Omax){
    p5[,o] = rnorm(N,wet[o,1],wet[o,2])
  }
  p6 <- array(dim=c(N,Omax)) #: p6 represent the regression parameter 1 a1
  for (o in 1:Omax){
    p6[,o] = rnorm(N,a1[o,1],a1[o,2])
  }
  p7 <- array(dim=c(N,Omax)) #: p7 represent the regression parameter 2 a2
  for (o in 1:Omax){
    p7[,o] = rnorm(N,a2[o,1],a2[o,2])
  }
  p8 <- array(dim=c(N,Omax)) #: r8 represent the regression parameter 3 a3
  for (o in 1:Omax){
    p8[,o] = rnorm(N,a3[o,1],a3[o,2])
  }
for (o in 1:Omax){
  for (t in 1:T){
    #.....STEP 3.1: estimate the statistical properties of chloride Ccl at any time t using random generation (e.g. Monte Carlo simulation or random sampling)
    for (n in 1:N){
      Cclt[t,n,o] = Ccl(i1[n,o],d[o],i3[n,o],i4[n,o],kc[o],i2[n,o],t0[o],i5[n,o],t)
      Ccrt[t,n,o] = crackwidth(p1[n,o],p2[n,o],p3[n,o],p4[n,o],p5[n,o],p6[n,o],p7[n,o],p8[n,o],d[o],phi[o],ftspl,t)
    }
    Ccltd[t,o] = Ccl(mean(i1[,o]),d[o],mean(i3[,o]),mean(i4[,o]),kc[o],mean(i2[,o]),t0[o],mean(i5[,o]),t)
    Ccrtd[t,o] = crackwidth(mean(p1[,o]),mean(p2[,o]),mean(p3[,o]),mean(p4[,o]),mean(p5[,o]),mean(p6[,o]),mean(p7[,o]),mean(p8[,o]),d[o],phi[o],ftspl,t)    
    #.....STEP 3.2: Calculating the state probability for each phase
    #.....for chloride
      #..................................................
    for (k in 1:I[o]){
      if (k==1){
        pi_cl[t,k,o] = (length(which(Cclt[t,,o]<a[o,(k+1)])))/N
      } else if (k<I[o]) {
          pi_cl[t,k,o] = (length(which(Cclt[t,,o]<a[o,(k+1)])))/N - sum(pi_cl[t,(1:(k-1)),o])
        } else if (k==I[o]) {
            pi_cl[t,k,o] = 1 - sum(pi_cl[t,(1:(k-1)),o])
          }    
    }
    #for propagation
    #....................................................
    for (k in 1:J[o]){
      if (k==1){
        pi_cr[t,k,o] = (length(which(Ccrt[t,,o]<b[o,(k+1)])))/N
      } else if (k<J[o]) {
          pi_cr[t,k,o] = (length(which(Ccrt[t,,o]<b[o,k+1])))/N - sum(pi_cr[t,(1:(k-1)),o])
        } else if (k==J[o]) {
            pi_cr[t,k,o] = 1 - sum(pi_cr[t,(1:(k-1)),o])
          }  
    }
    #.....STEP 3.3: calculating the state probability for the entire system
    #.....combining two probability
    for (k in 1:(I[o]+J[o]-1)){
      if (k<I[o]){
        pi[t,k,o] = pi_cl[t,k,o]
      } else {
          pi[t,k,o] = pi_cl[t,I[o],o]*pi_cr[t,(k-I[o]+1),o]
        }
    }
  }
}
#####################################################
detl <- array(dim=c(T,Omax)) #: deterioration process for the inition phase
for (o in 1:Omax){
  for (t in 1:T){
    detl[t,o] = 0
  }
}
for (o in 1:Omax){
  for (t in 1:T){
    if (I[o]>2){
      for (i in 1:(I[o]-1)){
        if (i==1){
          if (Ccltd[t,o]<a[o,2]){
            detl[t,o] = 1
          }
        } else {
          if (detl[t,o]==0){
            if (Ccltd[t,o]<a[o,(i+1)]){
              detl[t,o] = i
            } else {
              detl[t,o] = I[o]
            }
          } else {  
          }
        }
      }
    } else {
      for (i in 1:(I[o]-1)){
        if (Ccltd[t,o]<a[o,(i+1)]){
          detl[t,o] = i
        } else {
          detl[t,o] = I[o]
        }
      }
    }
  }
}    
detr <- array(dim=c(T,Omax)) #: deterioration process for the propagation phase
for (o in 1:Omax){
  for (t in 1:T){
    detr[t,o] = 0
  }
}
for (o in 1:Omax){
  for (t in 1:T){
    if (J[o]>2){
      for (j in 1:(J[o]-1)){
        if (j==1){
          if (Ccrtd[t,o]<b[o,2]){
            detr[t,o] = 1
          }
        } else {
          if (detr[t,o]==0){
            if (Ccrtd[t,o]<b[o,(j+1)]){
              detr[t,o] = j
            } else {
              detr[t,o] = J[o]
            }
          } else {  
          }
        }
      }
    } else {
      for (j in 1:(J[o]-1)){
        if (Ccrtd[t,o]<b[o,(j+1)]){
          detr[t,o] = j
        } else {
          detr[t,o] = J[o]
        }
      }
    }
  }
}
det <- array(dim=c(T,Omax)) #: combined deterioration process
for (o in 1:Omax){
  for (t in 1:T){
    if (detl[t,o]==I[o]){
      det[t,o] = detr[t,o] + (I[o]-1)
    } else {
      det[t,o] = detl[t,o]
    }
  }
}
#####################################################
cat("state probability of the initiation phase \n")
print(pi_cl)
cat("state probability of the propagation phase \n")
print(pi_cr)
cat("state probability of the entire process \n")
print(pi)

#------------------------------SAVE THE RESULTS------------------------------
file.remove("pi1.csv")
file.create("pi1.csv")
write.table(pi[,,1], file = "pi1.csv", sep = ",", append = TRUE,col.names = FALSE) 
#------------------------------THE END------------------------------