markovprob<-function(i,j,z,theta,thetasa){
probb<-matrix(double(1),jmax,jmax)
#theta<-matrix(double(1),nrow=1,ncol=jmax)
#thetasa<-matrix(double(1),nrow=jmax,ncol=jmax)
#################################################

prob1=0.0
reserve<-1.0
for (k in i:(j-1)){
if (j<=i) {
reserve=1
} else {
reserve=reserve*theta[k]
}
}
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
probb<-prob1
}
return(probb)
}
#########################



