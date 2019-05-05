require(fOptions)
source("bioption.R")
y1=-10
y2=10
#input value
AH=53000 # heated area in m2
AF=26000 #m2 facade surface area in m2
Fuel=10 #kW fuel demand in I per KW
Rent=1084 #CHF rent revenue
COST=7000 # CHF/m2
heatdemand=40 # Kwh/m2
sigma=0.3 #votality
P=0.95 #CHF/l initial fuel price
r=0.03 #discount factor
inflation=0.02 #inflation rate
years=10


Value1=AH*Rent-Fuel*AH*heatdemand*P


CRRTree = bioption(TypeFlag = "ce", S = P, X = P, 
     Time = years, r = r, b = 0.1, sigma = sigma, f=inflation,n = years)
   BinomialTreePlot(CRRTree, dy = 1, cex = 0.8, ylim = c(y1, y2),
     xlab = "n", ylab = "Option Value")
   title(main = "Option Tree")


  