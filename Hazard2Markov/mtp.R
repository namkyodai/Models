mtp<-function(jmax,z,theta){
thetasa<-matrix(double(1),nrow=jmax,ncol=jmax)
probb<-matrix(double(1),jmax,jmax)
thetasa<-matrix(double(1),nrow=jmax,ncol=jmax)
#################################################
for (i in 1:jmax){
for (j in 1:jmax){
thetasa[i,j]=theta[i]-theta[j]
}
}
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
return(probb)
}
#########################
