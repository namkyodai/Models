# File      : lst.R 
# Creation  : 09 Jul 2015
# Time-stamp: <Fre 2015-07-09 14:26 Nam>
#
# Copyright (c) 2015 Nam Lethanh <namkyodai@gmail.com>
#               http://www.ibi.ethz.ch
# $Id$ 
#
# Description : Unrestricted conventional Least Square estimator for Markovian transition probability
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>. 


#Least Square fit to the data to compute the transition matrix from proportional data

#Example taken from Miller (1952) - Finite Markov Process in psychology, Vol. 17 issue 2, p149-167

#function to compute the transition matrix
mtp<-function(a,b){
  c<-a%*%t(b)%*%solve(b%*%t(b))
  return(c)
}

#data
#data <- read.csv("example.csv",header=FALSE,sep=",") 
data <- read.csv("inp_cigarette.csv",header=FALSE,sep=",") 
#data <- read.csv("inp_sim.csv",header=FALSE,sep=",") 

attach(data)
Nmax<-ncol(data) #maximum number of condition state
Kmax<-length(data[,1]) # (K-1) number of trials, K is the total number of trial

M<-matrix(nrow=Nmax,ncol=(Kmax-1))
N<-matrix(nrow=Nmax,ncol=(Kmax-1))
CN<-matrix(nrow=Nmax,ncol=(Kmax-1)) #correction state probability
u<-matrix(nrow=Nmax,ncol=(Kmax-1)) #error term from 2 to T

P<-matrix(nrow=Nmax,ncol=Nmax)

solieu<-t(data)
#load data into function
for (k in 1:Kmax){
  if (k < Kmax){
  M[,k]<-solieu[,k]
  }
  if (k>1){
  N[,k-1]<-solieu[,k]
  }
}
p<-mtp(N,M)

P<-t(p)

cat("Transition probability matrix P \n")
print(P)

#calculating the corrected state probability and error term

for (k in 1:(Kmax-1)){
  if (k==1){
         CN[,k]<-M[,1]%*%P
      u[,k]<-CN[,k]-N[,k]
  } else {
        CN[,k]<-N[,k-1]%*%P
        u[,k]<-CN[,k]-N[,k]
}
}

cat("expected value of state probability \n")
print(CN)
cat("error term u \n")
print(u)

#calculate the squared deviations (covariance matrix)

cat("covariance matrix \n")
uu<-u%*%t(u)
print(uu)

cat("the best estimate dispersion of the calculated from the observed value is \n")

delta<-sqrt(uu[1,1]/(Kmax-Nmax-1))
print(delta)

cat("Variance-Covariance matrix V is  \n")

V<-delta^2*solve(M%*%t(M))
print(V)

#question????: Why there is negative value in the transition probability matrix
