##setwd("/home/namazu/Documents/Markovmodelwith-R/")
##source("markovgiventheta.R")
#This program was coded by Nam Lethanh from Osaka University
############INPUT PART############################
jmax=5  # Please select the maximum numbers of condition states in your system
theta<-c(0.158405202,0.061781277,0.026189521,0.009189522,0)
#theta<-c(0.205217083,0.092037386,0.045462598,0.019722855,0) # Please select values of hazard rate corresponding to each condition states
#jmax=3
#theta<-c(0.1,0.08,0)
z=2  # Please select the interval or elapsed time in Markov processs
########################################
thetasa<-matrix(double(1),nrow=jmax,ncol=jmax)
probability<-matrix(double(1),jmax,jmax)
##############defining thetasa value#################

##############################################
markovprob<-function(jmax,z,theta,probb){
probb<-matrix(double(1),jmax,jmax)
##theta<-matrix(double(1),nrow=1,ncol=jmax)
thetasa<-matrix(double(1),nrow=jmax,ncol=jmax)
#################################################
for (i in 1:jmax){
for (j in 1:jmax)
thetasa[i,j]=theta[i]-theta[j]
}
print(thetasa)
###############################################
for (i in 1:jmax){
for (j in 1: jmax){
prob1=0.0
reserve<-1.0
for (k in i:(j-1)){
if (j<=i) {
reserve=1
} else {
reserve=reserve*theta[k]
}
}
print(reserve)
#################################################
if (i>j){
probb[i,j]=0.0
} else {
for (k in i:j){
prod11=1.0
######################
for (e in i:j){
if(e !=k) {
prod11=thetasa[e,k]*prod11
}
}
#####################
prob1=prob1+exp(-theta[k]*z)/prod11
}
#############
prob1<-prob1*reserve
probb[i,j]=prob1
}
}
}
print(probb)
}
#########################
pro<-markovprob(jmax,z,theta,probb)
# After running this code in R, from R console, you just type down pro, results will appear.