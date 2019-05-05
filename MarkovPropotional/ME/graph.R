# This program was coded by  Nam Lethanh (lethanh@ibi.baug.ethz.ch)
#Steady state Markov model
#Intervention strategy for bridge with Markov model
#------------------------------------------------------
rm(list=ls()) #clear the memory and objects in R console


data <- read.csv("data.csv",header=FALSE,sep=",") #for example 2
stp<-c(1,0,0,0,0) #state probability #for example 3

attach(data)
Nmax<-length(data[,1]) #maximum number of condition states
P<-array(dim=c(Nmax,Nmax)) #define the dimension for transition matrix P
YearMax=100 #no of year
cat('# -------------------INPUT ------------------------------- \n')
#read value of transition matrix from csv file
for (i in 1:Nmax){
  for (j in 1:Nmax){
    P[i,j]=data[i,j]
  }
}
cat('value of transition matrix P \n')
print(P)
cat('# -------ESTIMATION ------------ \n')
#STEP 2: Define the state probabilitpi
pi<-array(dim=c(YearMax,Nmax)) #this is the state probability
for (t in 1:YearMax){
  if (t==1){
    pi[t,]=stp
    cat('Value of state probabiliy at time t \n')
    print(pi[t,])
    
  } else {
  pi[t,]= pi[t-1,]%*%P[,]   
  }
}
plot.new()
stackedPlot <- function(data, time=NULL, col=1:length(data),...) {
  if (is.null(time)) {
    time <- 1:length(data[[1]])
    plot(0,0, xlim = range(time), ylim = c(0,max(rowSums(data))), t="n" , ...)
    for (i in length(data):1){
      #the summup to the current collumn
      prep.data <- rowSums(data[1:i]);     
      # The polygon must have his first and last point on the zero line
      prep.y <- c(0, prep.data,0)
      prep.x <- c(time[1], time, time[length(time)] )
      polygon(prep.x, prep.y, col=col[i], border = NA);
    }
  }
}
#..............................................................
plot.new()
#define evolution of conditison states over time into data frame
csstate<-matrix(double(1),nrow=T,ncol=(Nmax)) #state probability in year
csstate=pi #turning probability distribution of condition states to percentage
csstate <- data.frame(csstate)
colors=c("white","coral1","limegreen","orange","yellow2")
stackedPlot(csstate,col = colors,xlab="",ylab="") #this is used for defaul colour
mtext(expression(paste("Probability")),side=2,line=2.2, adj = 0.5)
mtext(expression(paste("Time (years)")),side=1,line=2.2, adj = 0.5)
legend("topright", inset=0.09, title="CSs",col=colors,lty=2,lwd=13,legend=c(1:Nmax),bg="azure2",cex=0.8)
file.remove("pi.csv")
file.create("pi.csv")
write.table(csstate, file="pi.csv", sep = ",", append = TRUE,col.names = FALSE) 
#THE END




cat("THE END")






