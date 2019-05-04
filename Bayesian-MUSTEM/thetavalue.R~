##################################################
##BEGIN OF FUNCTION
thetavalue<-function(ik,jk,xk,beta){
valtheta<-matrix(double(1),nrow=1, ncol=jmax)
for (j in 1:jmax){
   if (j>=ik && j<=jk && j!=jmax){
  for (m in 1:mmax){ 
valtheta[j]<-valtheta[j]+beta[m,j]*xk[m,j]
}
}
 else {
valtheta[j]<-0

}
}
return(valtheta)
}
###############################################
thetasavalue<-function(theta){
for (i in 1:jmax){
  for (j in 1: jmax){
thetasa[i,j]<-theta[i]-theta[j]
}
}
return(thetasa)
}
#############################################


##END OF FUNCTION THETAVALUE
