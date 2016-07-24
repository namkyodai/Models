# This R code is the Financial Model used to evaluate the feasibility of Investment to build a new Thermal Power Plant.
# Developed by Nam Lethanh (24 May, 2016)
#All righr reserved
#Distributed as an open source, please use it at your own risk

#Verion 2.0

#Loading dependent package
library(date)
library(FinCal) #calling the financial package
#*********************************************************
#**********************************************************
#**********************************************************
# --------------INPUT-------------------------------
#**********************************************************
#**********************************************************
#**********************************************************
#**********************************************************

# 1) General information
T<-39 #Total time to considered
TT<-31 # Total numbers of year of the project life cycle
Tc <- 3 # construction to be at year 3 (start to construction on 09 June 2018)
Time<-4 # numbers of year require for construction work
hours<-24 #Number of hours to operate the plant in a day

Startdate<-as.date("1Jan2016") # date to consider in the project

Opedate<-as.date("08August2021") # date of operation of the plant
LastdayofOpeYear<-as.date("31Dec2021")
StartYear<-2016
OpeYear<-2021

#FuelC<-2.05 #($/GJ)
#Epr<- 85 # ($/Mwh)

marketinput <- read.csv(file="marketinput.csv", header=TRUE, sep=",")
marketinput<-data.frame(marketinput)

#------------------------
DayRemain<-LastdayofOpeYear-Opedate # the remaining 

LeapYear<-matrix(double(1),nrow=1,ncol=T) #Leap Year
# leap year occurs every 4 year. Thus for the first 4 year, binary value must be defined.
for (t in 1:T){
  if (t==1) {
    LeapYear[t] = 1 } # this has to be defined.
  else if (t==2){
    LeapYear[t] = 0 } # this has to be defined.
  else if (t==3){
    LeapYear[t] = 0 } # this has to be defined.
  else if (t==4){
    LeapYear[t] = 0 } # this has to be defined.
  else {
    LeapYear[t] = LeapYear[t-4]*1 
  } # this has to be defined.
}


YearPlan<-matrix(double(1),nrow=3,ncol=T) #Leap Year
for (t in 1:T){
  if (t==1){
    YearPlan[1,t] <- StartYear
    YearPlan[2,t] <- LeapYear[t]
    YearPlan[3,t] <- 365+LeapYear[t]
  } else {
    YearPlan[1,t]<-YearPlan[1,t-1]+1
    YearPlan[2,t] <- LeapYear[t]
    YearPlan[3,t] <- 365+LeapYear[t]
  }
}
# 2)  Plant output and efficientcy
Ci <- 952.247   # Installed capacity (in MW)
Nf <- 0.944892   # Netoutput factor (%)
Eff <- 0.3679    # Efficiency % (HHV basic)
Ca <- 0.85        # Capacity factor (%)
Cafa <- 0.996        # change in yearly capacity factor
Ne <-900        # Net operating power output (MW)
Chr <- 3600/Eff # cycle heat rate (Kj/Kwh)
Coeff<-3.6 #Coefficient if WHAT?

# 3) Coal Heat values
FHHV1 <- 4221  # Fuel HHV (Kcal/kg)
FLHV1 <- 3889  # Fuel HHV (Kcal/kg)
FHHV2 <- FHHV1*4.185265/1000  # Fuel HHV (Mj/kg)
FLHV2 <- FLHV1*4.185265/1000  # Fuel HHV (Mj/kg)

limestoneP<-24.128

# 3) Allocation of investment
AoI<-matrix(double(1),nrow=1,ncol=Time) # In percentage of total investment
AoI[1]<- 0.19
AoI[2]<- 0.31
AoI[3]<- 0.31
AoI[4]<- 0.19


# 4) Allocation of investment
EPC0 <- 1399659092 # Power Plant EPC cost $
CsC <- 218064753 # Contractor soft cost $
OsC <- 125000000 # Owner soft cost $
FIDC <- 100000 # Financial and IDC cost $
WCA <- 500000# Working capital allowance $
ITPC <- EPC0+CsC+OsC+FIDC+WCA # Initial total project cost $
Occ <- EPC0/Ci/1000 # Overnight cost of construction ($/Kw)

# 5) Construction and Financing

Mlf <- 1.005 #Marginal lost factor 
Eqt <- 30/100 # equity as % of total project cost
Cir <- 6.5/100 # Construction interest rate (%)
FirSe <- 6/100 # Financing interest rate (%) senior debt
DbFeeSe <- 0.0225 # Debt fee (% of debt loan)
DbTermSe <- 15 # Debt term in number of years
FirSu <- 5.5/100 # Financing interest rate (%) subordinate debt
DbFeeSu <- 0.027 # Debt fee (% of debt loan)
DbTermSu <- 20 # Debt term in number of years
BuffCost <- 0 # Buffer for cost overruns, % of total amount to be financed

#  6) O&M
OMfix <- 12500000 # Fixed O&M cost($/year)
OMvar <- 4 # Variable O&M cost ($/Mwh)

#  7) Additional operating expenses
CoC <- 0.01 # Cost of consumables ($/Mwh)
Pex <- 500000 # Personnel expenses, ($/year)
MnR <-0.02 # equipment maintnenace and report cost ($/Mwh)
Ovh <- 50000 # General and admin expenses ($/year)
IIP <- 0.05/100 # Initial insurance premium (% of toal CAPEX)


#  8) Taxation
Ctr <- 30/100 # Corporate taxation rate
Dth <- 0 # Detail of tax holiday
Dwt <- 0 # Devidends withholding tax

# 9) Depreciation and residual value
DepTime <- 30 # Deperciation period, years (for tax purposes)
DepRate <- 1/DepTime # Depreciation rates
Res <- 20/100 # Plant residual value
 
# 10) Discount rate and return on equity

normalDisrate<-6/100

#ROE <- 0 # Return on Equity
TargetROE <-12/100 # Targeted Return on Equity

# 11) Escalation and exchange rate

Iflat <- 2.5/100 # Inflation rate (%) #not this this value can be 2.5

Disrate <- (normalDisrate+1)/(1+Iflat)-1 # discount rate (real term)


Iflatusd <- 2/100 # Inflation rate (%) in USA
EsOM <- 2/100 # O&M escalation rate (%)
Exc <- 0.02196 # Exchange rate 
ExcVar <- 0 # Exchange rate variation (%/year)


#*********************************************************
#**********************************************************
#**********************************************************
# --------------MODEL-------------------------------
#**********************************************************
#**********************************************************
#**********************************************************
#**********************************************************


# -----------INPUT GROWTH PATH-----------------------
# This section is about to calculate the annual OPEX of the plant base on unit price

PriceSet <-2 # this is a set of price (e.g. Normal price or Constant price in )


# -------Macroeconomic indicators

# Philippines CPI index (2016=1)
CPI <-matrix(double(1),nrow=1,ncol=T)
for (t in 1:T){
  if (t==1) {
    CPI[t] = 1 } # this has to be defined.
  else {CPI[t]=CPI[t-1]*(1+Iflat)
    }
  }

# USA CPI index (2016=1)
CPIusd <-matrix(double(1),nrow=1,ncol=T)
for (t in 1:T){
  if (t==1) {
    CPIusd[t] = 1 } # this has to be defined.
  else {CPIusd[t]=CPIusd[t-1]*(1+Iflatusd)
  }
}

# O&M cost index (2016=1)
OMIndix <-matrix(double(1),nrow=1,ncol=T)
for (t in 1:T){
  if (t==1) {
    OMIndix[t] = 1 } # this has to be defined.
  else {OMIndix[t]=OMIndix[t-1]*(1+EsOM)
  }
}

# Exchange rate forecasting
ExcAn <-matrix(double(1),nrow=1,ncol=T) #Annual exchange rate
for (t in 1:T){
  if (t==1) {
    ExcAn[t] = Exc } # this has to be defined.
  else {ExcAn[t]=ExcAn[t-1]*(1+ExcVar)
  }
}

# Fuel cost ($/GJ)
Fuel <-matrix(double(1),nrow=1,ncol=T) #Annual exchange rate

# if the fuel follow a linear model
#for (t in 1:T){
#  if (t==1){
#    Fuel[t]<-FuelC
#  } else {
#    Fuel[t]<-Fuel[t-1]*1 #1 is the coefficient
#  }
#}

# if the fuel is taken from the market study
for (t in 1:T){
Fuel[t]<-marketinput[t,2]
}


cat("Yearly fuel price \n")
print(Fuel)

# electricprice ($/GJ)
Energyprice <-matrix(double(1),nrow=1,ncol=T) #Annual exchange rate

#if energy follows a linear function
#for (t in 1:T){
#  if (t==1){
#    Energyprice[t]<-Epr
#  } else {
#    Energyprice[t]<-Energyprice[t-1]*1 #1 is the coefficient
#  }
#}

#if energy price is taken from the market study

for (t in 1:T){
  Energyprice[t]<-marketinput[t,4]
}



cat("Yearly fuel price \n")
print(Fuel)


# -----------PRODUCTIOn-----------------------

# operating hour per year

Production<-matrix(double(1),nrow=6,ncol=T) #Leap Year
NPVProduction<-matrix(double(1),nrow=2,ncol=T) #Leap Year
for (t in 1:T){
      if (YearPlan[1,t]-OpeYear < 0){
      Production[1,t]=0 #hour per year (hours)
      Production[2,t]=0 #capacity (%)
    } else if ((YearPlan[1,t]-OpeYear == 0)) {
      Production[1,t]=hours*DayRemain
      #Production[2,t]=Ca # if this is to follow a linear model
      Production[2,t]<-marketinput[t,3]/100} # if this is from the market outlook
      else {
      Production[1,t]=hours*YearPlan[3,t]
      #Production[2,t]=Production[2,t-1]*Cafa #if this is to follow a linear model
      Production[2,t]<-marketinput[t,3]/100 # if this is from the Market outlook
      }
      Production[3,t]=Production[1,t]*Production[2,t]*Ci*Nf/1000 #Energy sold (GWh)
      Production[4,t]= Energyprice[t]*CPIusd[t]# Energy price (Norminal usd/Mwh taking into consideration of CPI index)
      Production[5,t]<- Production[3,t]*Production[4,t]/Mlf/1000# Energy Market revenue (normial)
      Production[6,t]<- Energyprice[t]*Production[3,t]/Mlf/1000# Energy Market revenue (constant)
      Production[,T]<-Production[,T-1]
      NPVProduction[1,t]<-Production[5,t]/((1+Disrate)^(t-1)) #Yearly NPV of Energy Market Revenue (normial)
     NPVProduction[2,t]<-Production[6,t]/((1+Disrate)^(t-1)) #Yearly NPV of Energy Market Revenue (constant)
  }


cat("Summary of yearly production factors \n")
print(Production)

#stop("debugg")

# -----------OPERATION COST-----------------------

#fixed annual O&M cost

OM<-matrix(double(1),nrow=7,ncol=T) #Leap Year
NPVOM<-matrix(double(1),nrow=4,ncol=T) #Leap Year
for (t in 1:T){
  if (YearPlan[1,t]-OpeYear < 0){
    OM[1,t]=0 #hour per year (hours)
  } else {
    OM[1,t]=Production[3,t]/(Eff*FHHV2)*Coeff #Fuel consumption (kt)
    OM[2,t]=Production[3,t]*Coeff*marketinput[t,2]*CPIusd[t]/(1000*Eff*Nf) #Fuel cost (nominal USD millions normial)
    OM[3,t]=Production[3,t]*Coeff*marketinput[t,2]/(1000*Eff*Nf) #Fuel cost (nominal USD millions constant)
    OM[4,t]=OMfix/1000000*OMIndix[t]+Production[3,t]*(OMvar/1000*OMIndix[t]) #O&M normial
    OM[5,t]= OM[4,t]/CPIusd[t] #O&M constant
    OM[6,t]= 100*limestoneP*Production[2,t]*CPIusd[t]/90 #limestone cost (normial)
    OM[7,t]= OM[6,t]/CPIusd[t] #limestone cost (costatnt)
    
    OM[,T]=OM[,T-1]
    NPVOM[1,t]<-OM[2,t]/((1+Disrate)^(t-1)) #
    NPVOM[2,t]<-OM[3,t]/((1+Disrate)^(t-1)) #
    NPVOM[3,t]<-OM[4,t]/((1+Disrate)^(t-1)) #
    NPVOM[4,t]<-OM[5,t]/((1+Disrate)^(t-1)) #
  } 
}
cat("Summary of O&M factors \n")
print(OM)
#stop("aaaa")

OP<-matrix(double(1),nrow=12,ncol=T) #Other Operation Cost need to be considered or aggregated

for (t in 1:T){
  if (t<T){
  for (i in 1:12){
    if (Production[3,t] > 0 ){
      OP[1,t]<-CoC # yearly cost of consumables (usd/MWh) - constant price
      OP[2,t]<-Pex/1000000 # yearly personall expenses (million usd/year) - constant price
      OP[3,t]<-MnR # Equipment maintenance (usd/MWh) - constant price
      OP[4,t]<-Ovh/1000000 # yearly genera expenses (million usd/year) - constant price
      OP[5,t]<-OMfix/1000000*OMIndix[t] # yearly fixed annual O&M cost (million $) - normial price
      OP[6,t]<-OMvar/1000*OMIndix[t]  # yearly variable annual O&M cost (million $) - normial price
    } else {
      OP[1,t]<-0
      OP[2,t]<-0
      OP[3,t]<-0
      OP[4,t]<-0
      OP[5,t]<-0
     # OP[6,t]<-0
    }
      OP[7,t]<-OP[5,t]/CPIusd[t]  # yearly fixed annual O&M cost (million $) - constant price
      OP[8,t]<-OP[6,t]/CPIusd[t]  # yearly variable annual O&M cost (million $)- constant price
      OP[9,t]<-OP[1,t]*CPIusd[t] #normial price
      OP[10,t]<-OP[2,t]*CPIusd[t]#normial price
      OP[11,t]<-OP[3,t]*CPIusd[t]#normial price
      OP[12,t]<-OP[4,t]*CPIusd[t]#normial price
  }
  } else {
    for (i in 1:12){
      OP[i,t]  <-OP[i,t-1]
    }
    
  }
}
cat("Other OM cost \n")
print(OP)


#stop("deubb")

# -----------CAPITAL COST-----------------------

CP<-matrix(double(1),nrow=17,ncol=T) #Capital Cost
for (t in 1:T){
  if (t < Tc){
    CP[1,t]<- 0 # Annual ratio of allocation of investment  
  } else if (t < Tc+Time) {
    CP[1,t]<- AoI[t-Tc+1] # Annual ratio of allocation of investment  
  } else {
    CP[1,t]<- 0 # Annual ratio of allocation of investment  
  }
  
  #------------
  CP[2,t]<- CP[1,t]*EPC0*(1+BuffCost)/1000000 # Power Plant EPC cost million $
  CP[3,t]<- CP[1,t]*CsC*(1+BuffCost)/1000000 # Contractor soft cost million $
  CP[4,t]<- CP[1,t]*OsC*(1+BuffCost)/1000000 # Owner soft cost million $
  CP[5,t]<- CP[1,t]*FIDC*(1+BuffCost)/1000000 # Financial and IDC cost million $
  CP[6,t]<- CP[1,t]*WCA*(1+BuffCost)/1000000 # Working capital allowance million $
  CP[7,t]<- CP[2,t]+CP[3,t]+CP[4,t]+CP[5,t]+CP[6,t] # Total investment million $
  CP[8,t]<- CP[7,t]*(1-Eqt)*Cir # interest during construction million $
  CP[9,t]<- CP[7,t]*(1-Eqt) # Financing amount million $
  CP[10,t]<- CP[8,t]*(1-Eqt) #Financing amount million $
  CP[11,t]<- CP[7,t]-(CP[9,t]+CP[10,t]) # Initial equity million $
  } 

# ------
for (t in 1:T){
  if (t < Tc+Time){
    CP[12,t]=0 # Loand Repayment - senior loan (milliion usd)
  } else if (t<Tc+Time+DbTermSe){
    CP[12,t]= sum(CP[9,])/DbTermSe#
  } else {
    CP[12,t]= 0
  }
}

# -----
for (t in 1:T){
  if (t < OpeYear+DbTermSe-StartYear+1){
    CP[13,t]=sum(CP[9,seq(1:t)])-max(CP[12,t]*(YearPlan[1,t]-OpeYear),0) # Loan balance end of year - senior loan (usd million)
  } else  {
    CP[13,t]= 0
  }
  if (t < Tc+Time){
    CP[14,t]=0 # Interest - senior loand (usd million)
  } else if (t<Tc+Time+DbTermSe){
    CP[14,t]= CP[13,t-1]*FirSe# 
  } else {
    CP[14,t]= 0
  }
}
# ----

for (t in 1:T){
  if (t < Tc+Time){
    CP[15,t]=0 # Loand Repayment - senior loan (milliion usd)
  } else if (t<Tc+Time+DbTermSu){
    CP[15,t]= sum(CP[10,])/DbTermSu#
  } else {
    CP[15,t]= 0
  }
}


# -----
# Interest - senior loan (usd million)
for (t in 1:T){
  if (t < OpeYear+DbTermSu-StartYear+1){
    CP[16,t]=sum(CP[10,seq(1:t)])-max(CP[15,t]*(YearPlan[1,t]-OpeYear),0) # Loan balance end of year - senior loan (usd million)
  } else  {
    CP[16,t]= 0
  }
  if (t < Tc+Time){
    CP[17,t]=0 # Interest - senior loand (usd million)
  } else if (t<Tc+Time+DbTermSu){
    CP[17,t]= CP[16,t-1]*FirSu# 
  } else {
    CP[17,t]= 0
  }
}


CP[,T]<-CP[,T-1] #residual of the last year equal to its year before.

#now time to calculate the NPV for each of the capital factors

NPVCP<-matrix(double(1),nrow=16,ncol=T) #
for (t in 1:T){
  for (i in 1:16){
    NPVCP[i,t]<-CP[i+1,t]/((1+Disrate)^(t-1)) #
  }
}


cat("Summary of Capital Cost factors \n")
print(CP)

#stop("debugg")

# ------------NPV values
NPV<-matrix(double(1),nrow=22,ncol=1) # NPV value of various factor
for (i in 1:22){
  if (i < 3){
    NPV[i]<-sum(NPVProduction[i,])
  } else if (i < 7){
    NPV[i]<-sum(NPVOM[i-2,])
  } else {
    NPV[i]<-sum(NPVCP[i-6,])
  }
}  
cat("NPV of key factors \n")
print(NPV,digits =16)

#Reframe NPV

#stop("Debugg")

#*********************************************************
#**********************************************************
#**********************************************************
# --------------OUTPUT-------------------------------
#**********************************************************
#**********************************************************
#**********************************************************
#**********************************************************


# -----INCOME STATEMENT


InNo<-21 # total number of items in the Income statement sheet

INS<-matrix(double(1),nrow=InNo,ncol=T) #
INSNPVyear<-matrix(double(1),nrow=InNo,ncol=T) # NPV of INS

for (t in 1:T){
  if (t<T){
  for (i in 1:InNo){
    #Sales Revenue
    INS[1,t]<-Production[5,t] #energy market 
    INS[2,t]<-INS[1,t] # Total sales revenue, if exists only one market (A)
    #Cost of Sale
    INS[3,t]<-OM[2,t] # Fuel
    INS[21,t]<-OM[6,t] # limestone
    INS[4,t]<-OM[4,t] # O&M
    INS[5,t]<-INS[3,t]+INS[4,t]+INS[21,t] # Total cost of Sale (B)
    #Gross Profit
    INS[6,t]<-INS[2,t]-INS[5,t] # Gross Profit = Total savel revenue - Total cost of sale (C=A-B)
    # General and Administrative
    INS[7,t]<-OP[9,t]*Production[3,t]/1000 #Cost of consumables
    INS[8,t]<-OP[10,t] # personal expense
    INS[9,t]<-OP[11,t]*Production[3,t]/1000 # equipment M&R
    INS[10,t]<-OP[12,t] # other general and admin expenses
    INS[11,t]<-CP[9,t]*DbFeeSe+CP[10,t]*DbFeeSu # Debt Fee
    INS[12,t]<-IIP*sum(CP[2,seq(1:t)])*CPIusd[t] # Insurance
    INS[13,t]<- INS[7,t]+INS[8,t]+INS[9,t]+INS[10,t]+INS[11,t]+INS[12,t]# Total General and Admin Expenses ---CHECK THIS (D)
    #Taxable Income (E=C-D)
    INS[14,t]<-INS[6,t]-INS[13,t] 
    INS[15,t]<-CP[14,t]+CP[17,t] #Interest expense (F)
    if (t >= (OpeYear-StartYear+1) & t<=(OpeYear-StartYear+1)+DepTime){
    INS[16,t]<-(sum(CP[7,seq(1:t)])+sum(CP[11,seq(1:t)]))/DepTime*(1-Res) #Depreciation (G)
    } else {
    INS[16,t]  <-0
    }
    #Tax
    INS[17,t]  <-max(INS[14,t]*Ctr,0) # income taxes
    INS[18,t]  <-INS[17,t] # Total taxes as the sume of all taxes. here theer is only one so there is no need to sum (H)
    INS[19,t] <- INS[14,t]- INS[15,t]-INS[16,t]-INS[18,t]# Net After Tax Income (I=E-F-G-H) 
    INS[20,t]   <-FirSe*sum(CP[7,t]) # Interest during construction
  }
  } else {
    for (i in 1:InNo){
      INS[i,t]<-INS[i,t-1] # this is the value of residual
    }
  } 
}

cat("----------------------------------------------------------------- \n")
cat ("INCOME STATEMENT \n")
cat("----------------------------------------------------------------- \n")
print(INS)

#stop("debugg")
#  to calculate the NPV of income statement at each year

for (t in 1:T){
  for (i in 1:InNo){
    INSNPVyear[i,t]<-INS[i,t]/((1+Disrate)^(t-1)) ##
  }
}
# ------------NPV values
INSNPV<-matrix(double(1),nrow=InNo,ncol=1) # NPV value of various factor
for (i in 1:InNo)    {
  INSNPV[i]<-sum(INSNPVyear[i,seq(1:TT)])
}
cat("NPV of Income statement \n")
print(INSNPV)


# -----CASH FLOW


cat("---------------------------------------------------------- \n")
cat ("Cash Flow Factors \n")
cat("---------------------------------------------------------- \n")
CFno<-9
CF<-matrix(double(1),nrow=CFno,ncol=T) #Casd flow
CFNPVyear<-matrix(double(1),nrow=CFno,ncol=T) #NPV of CF

for (t in 1:T){
  if (t <T){
  for (i in 1:CFno){
    if (i==1){
      CF[i,t]<--(CP[7,t]+CP[8,t]) # Investment cost
    } else if (i==2){
      CF[i,t]<--(OM[3,t]) # Fuel
    } else if (i==3){
      CF[i,t]<--(OM[5,t]) # O&M
    } else if (i==4){
      CF[i,t]<-(sum(CF[1,seq(1:t)]))*IIP # Insurance
    } else if (i==5){
      CF[i,t]<--((OP[1,t]+OP[3,t])*Production[3,t]/1000 +(OP[2,t]+OP[4,t])) #overhead
    } else if (i==6){
      CF[i,t]<-Production[6,t] #value of energy sales
    } else if (i==7){
      CF[i,t]<--INS[18,t] #tax
    } else if (i==8){
      CF[i,t]<--(OM[6,t])  # limestone 
    } else {
      CF[i,t]<-CF[1,t]+CF[2,t]+CF[3,t]+CF[4,t]+CF[5,t]+CF[6,t]+CF[7,t]+CF[8,t]# Net after tax cashflow
    }
  }
  } else {
    for (i in 1:CFno){
    CF[i,t]<-CF[i,t-1]
    }
  }
}
print(CF)
#  to calculate the NPV of cash flow at each year

for (t in 1:T){
  for (i in 1:CFno){
    CFNPVyear[i,t]<-CF[i,t]/((1+Disrate)^(t-1)) ##
  }
}
# ------------NPV values
CFNPV<-matrix(double(1),nrow=CFno,ncol=1) # NPV value of various factor
for (i in 1:CFno)    {
  CFNPV[i]<-sum(CFNPVyear[i,seq(1:T)])
}
cat("NPV of cash flow \n")
print(CFNPV)


cat("---------------------------------------------------------- \n")
cat ("BALANCE SHEET \n")
cat("----------------------Asset----------------------------- \n")

Assetno<-15
Asset<-matrix(double(1),nrow=Assetno,ncol=T) #
AssetNPVyear<-matrix(double(1),nrow=Assetno,ncol=T) #
for (t in 1:T){
  for (i in 1:Assetno){
    #Current Asset
    if (i==1){
      Asset[i,t]<-WCA/1000000*CPIusd[t] # Yearly Cash
    } else if (i==2){
      Asset[i,t]<-0 #Investments
    } else if (i==3){
      Asset[i,t]<-0 # Inventories
    } else if (i==4){
      Asset[i,t]<-0 # Accounts receivable
    } else if (i==5){
      Asset[i,t]<-0 # Loans payable and current portion long-term debt
    } else if (i==6){
      Asset[i,t]<-0 # Pre-paid expenses
    } else if (i==7){
      Asset[i,t]<-sum(Asset[seq(1:i-1),t]) #Total assets
    } else if (i==8){
    #  Fixed Asset
      if (t==1){
      Asset[i,t]<-(CP[7,t]+CP[8,t]) # Property and equipment
      } else {
      Asset[i,t]<-(CP[7,t]+CP[8,t])+Asset[i,t-1] # Property and equipment  
      }
    } else if (i==9){
      Asset[i,t]<-0 # Leasehold improvements
      } else if (i==10){
      Asset[i,t]<-0 # Equity and other investments
      } else if (i==11){
      if (t==1){
      Asset[i,t]<-INS[16,t] # Less accumulated depreciation
      } else {
      Asset[i,t]<-min(INS[16,t]+Asset[i,t-1],Asset[8,t]) # Less accumulated depreciation
      }
      } else if (i==12){ # Net fixed assets
      Asset[i,t]<-max(Asset[8,t]+Asset[9,t]+Asset[10,t]-Asset[11,t],0)
      # Other Assets
      } else if (i==13) { #goodwill
      Asset[i,t] <-0 
      } else if (i==14){
        Asset[i,t] <-Asset[i-1,t]
      } else if (i==15){
        Asset[i,t] <-Asset[7,t]+Asset[12,t]+Asset[14,t]# Total assets 
    }
  }
}
Asset[,T]<-Asset[,T-1]

print(Asset)
#  to calculate the NPV of asset at each year

for (t in 1:T){
  for (i in 1:15){
    AssetNPVyear[i,t]<-Asset[i,t]/((1+Disrate)^(t-1)) ##
  }
}
# ------------NPV values
AssetNPV<-matrix(double(1),nrow=15,ncol=1) # NPV value of various factor
for (i in 1:15)    {
  AssetNPV[i]<-sum(AssetNPVyear[i,seq(1:T)])
}
cat("NPV of Asset \n")
print(AssetNPV)



cat("----------------------reliabilities and Onwer's equity----------------------------- \n")

Relino<-15
Reliability<-matrix(double(1),nrow=Relino,ncol=T) #
ReliabilityNPVyear<-matrix(double(1),nrow=Relino,ncol=T) #
for (t in 1:T){
  for (i in 1:15){
    if (i<8){
      Reliability[i,t]<-0 # this assumption holds true unless there is a change.
    } else if (i==8){
      Reliability[i,t]<-sum(Reliability[seq(1:7),t]) # this is the sum of reliabilities and owner's equity
    } else if (i==9){ #Long-term reliabilities
      Reliability[i,t]<-CP[13,t] #Long term debt
    } else if (i==10){
      Reliability[i,t]<-0 # Deferred income tax
    } else if (i==11){
      Reliability[i,t]<-Reliability[9,t]+Reliability[10,t] # Total long-term reliabilities
    } else if (i==12){
      Reliability[i,t]<-CP[11,t] # Investment capital
    } else if (i==13){
      Reliability[i,t]<-Asset[15,t]-Reliability[11,t]-Reliability[8,t] #total owner's equity
    } else if (i==14){
      Reliability[i,t]<-Reliability[13,t]-Reliability[12,t] #accmulated retained earnings
    } else if (i==15){
      Reliability[i,t]<-Reliability[13,t]+Reliability[11,t]+Reliability[8,t] # Total reliabiliteis and owners equity
    }
  }
}

Reliability[,T]<-Reliability[,T-1] # residual of reliability in the last year
print(Reliability)


#  to calculate the NPV of reliability at each year

for (t in 1:T){
  for (i in 1:15){
    ReliabilityNPVyear[i,t]<-Reliability[i,t]/((1+Disrate)^(t-1)) ##
  }
}
# ------------NPV values
ReliabilityNPV<-matrix(double(1),nrow=Relino,ncol=1) # NPV value of various factor
for (i in 1:Relino)    {
  ReliabilityNPV[i]<-sum(ReliabilityNPVyear[i,seq(1:T)])
}
cat("NPV of Reliability \n")
print(ReliabilityNPV)


cat("---------------------------------------------------------- \n")
cat("----------------------COVERAGE RATIOS---------------------- \n")
cat("---------------------------------------------------------- \n")
 
COVno<-4 # total numbers of coverage factors
COCRatio<-matrix(double(1),nrow=COVno,ncol=T) 
for (t in 1:T){
  if (Reliability[15,t] != 0){
  COCRatio[1,t]<-INS[19,t]/Reliability[15,t]
  } else {
    COCRatio[1,t]<-0
  }
  COCRatio[2,t]<-Asset[7,t]-Reliability[8,t] # Working capital (millions)
  COCRatio[3,t]<-(Reliability[9,t]+Asset[5,t])/Reliability[13,t] # = (long term debt+loans payable and current portion long-term debt)/Total onwer equity
  COCRatio[4,t]<-(Reliability[9,t]+Asset[5,t])/Asset[15,t] # Debt ratio
}
COCRatio[4,T]<-COCRatio[4,T-1] #residual values
print(COCRatio)

# calculating present value of coverage ratios


COCRatioNPV<-matrix(double(1),nrow=COVno,ncol=1) 
cat("------------------Present value of Coverate ratios------------------ \n")


COCRatioNPV[1]<-INSNPV[19]/ReliabilityNPV[15]*100 # Return on equity
COCRatioNPV[2]<-AssetNPV[7]-ReliabilityNPV[8] # Working capital (millions)
COCRatioNPV[3]<-(AssetNPV[5]+ReliabilityNPV[9])/ReliabilityNPV[13] # Debt-to-equity ratio
COCRatioNPV[4]<-(AssetNPV[5]+ReliabilityNPV[9])/AssetNPV[15] # Debt ratio

print(COCRatioNPV)

cat("---------------------------------------------------------- \n")

cat("----------------CALCULATION THE INTERNAL RATE OF RETURN (IRR)-------------- \n")

IRR<-irr(CF[9,]*CPI)*100
cat("value of IRR is \n")
print(IRR)
cat("---------------------------------------------------------- \n")
cat("---------------------------------------------------------- \n")


# Stacked Bar Plot with Colors and Legend

#barplot(CF[seq(1:7),], main="Cash Flow",xlab="Years", col=c("gold2","lightgreen","deepskyblue1","darkorchid1","firebrick2","gray61","indianred1"),legend = c("Investment","Fuel","O&M","Insurance","Overheads","Energy sales", "Tax"))

#box()
