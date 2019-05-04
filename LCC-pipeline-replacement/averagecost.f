*******************************************************
! This program is for estimating the parameter of Weibull hazard function
! and estimating the life cycle cost, optimal year of renewal for
! pipelines system
! Paper name: Optimal Renewal of Water Supply Pipeline
! Written and programed by Nam, Le Thanh on April, 30th 2008
! Re-programed and mofified on Nov, 23th, 2008
*****************************************************************
	Module variable
	real*8:: lcc
!	the input file must be checked carefully for the exact No. of data
	character(len=19), parameter:: file100="input\rho.csv" ! 
	character(len=19), parameter:: file110="input\ccost.csv" !
	character(len=19), parameter:: file120="input\icost.csv" !	
	character(len=19), parameter:: file130="input\solieu.csv" !

	character(len=19), parameter:: file200="output\arho.csv" !
	character(len=19), parameter:: file210="output\aopyearc.csv" !
	character(len=19), parameter:: file220="output\aopyeari.csv" !
	character(len=19), parameter:: file230="output\ketqua.csv" !

!	real*8:: lcc1(n)

	End module
***********************************************************************
! Main program starts here
***********************************************************************	
	Program optimalrenewaltime
	use variable; implicit none
	real*8:: m, eps, obj, r, alpha1, m1, lcc2
	real*8:: c1A, c2A, c1C, c2C, c1F, c2F, c1FL, c2FL
	real*8:: r1(100), c1(100), c2(100)
	real*8:: gamma(100), lamda(100)

	real*8:: alpha, hess(2,2), ts(2), optts(2), df(2), grad(2), tv(2)
!	real*8, dimension(2):: ts
	integer:: k, time, i, j, eval, dim, ii, iii

!	Here new definition for LCC estimation
	real*8:: x1r, x2r, tpr, minur, ddr, ttr, JJJr
	real*8:: rr(100), cir, AC(100), MinAC, ACC
!	real*8, dimension(n):: ur
	integer:: kr

	! first to read the data from input file
!	alpha=3.44E-05
!	m=2.144000053
	
	alpha=1.24E-05
	m=2.484087137
	
	open(110, file="input\ccost.csv")
	open(130, file="input\solieu.csv")
	open(230, file="output\ketqua.csv")
	
!	stop
	
	do ii=1, 100
	read(110,*) r1(ii), c1(ii), c2(ii)
	end do
	
	do iii=1, 100
		read(130,*) gamma(iii), lamda(iii)
	

	end do
	
	do ii=1, 100
		do iii=1, 100
	AC(iii)= ((c1(ii)+c2(ii))*alpha*m*gamma(iii)
     #	+lamda(iii))/iii
	end do
	
	MinAC=AC(1)
	Do iii=1, 100
	if (MinAC.GT.AC(iii)) MinAC=AC(iii)
	end do

	ACC=MinAC
	
	do iii=1, 100
	if (AC(iii).EQ.MinAC) k=iii
	Print*, k
	end do
	end do
	
	
!	Print*, AC


	
	
!	Print*, gamma





	End program optimalrenewaltime
*******************************************************************	
!	Pattern method to estimate the value of m and alpha (or ts(1) and ts(2)


