# This R code is the Financial Model used to evaluate the feasibility of Investment to build a new Thermal Power Plant.
# Developed by Nam Lethanh 
#All righr reserved
#Distributed as an open source, please use it at your own risk

# published in Github for general use

# for specific use of any project, users need to customize and use the project specific input files.

#Verion 3.0

#Loading dependent package
library(date)
library(FinCal) #calling the financial package
#*********************************************************
#**********************************************************
#**********************************************************
# --------------INPUT-------------------------------
#**********************************************************
marketinput <- read.csv(file="marketinput.csv", header=TRUE, sep=",")
marketinput<-data.frame(marketinput)
primaryinput<-read.csv(file="primaryinput.csv", header=TRUE, sep=",")
primaryinput<-data.frame(primaryinput)
datetobe<-read.csv(file="date.csv", header=TRUE, sep=",",stringsAsFactors=F)
datetobe<-data.frame(datetobe)

# subroutine - generation unit 1
source("financial-model-unit1.R") # this is the subroutine for unit 1
source("financial-model-unit2.R") # this is the subroutine for unit 1
source("financial-model-unit3.R") # this is the subroutine for unit 1

# Voila --> the End :)



