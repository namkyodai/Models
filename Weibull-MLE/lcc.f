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

	character(len=19), parameter:: file200="output\rho.csv" !
	character(len=19), parameter:: file210="output\opyearc.csv" !
	character(len=19), parameter:: file220="output\opyeari.csv" !

	integer, parameter:: n=3*10**4
	real*8:: lcc1(n)



	End module
		
***********************************************************************
! Main program starts here
***********************************************************************	
	Program optimalrenewaltime
	use variable; implicit none
	real*8:: m, eps, obj, r, alpha1, m1, lcc2
	real*8:: c1A, c2A, c1C, c2C, c1F, c2F, c1FL, c2FL
	real*8:: r1(100), c1(100), c2(100)

	real*8:: alpha, hess(2,2), ts(2), optts(2), df(2), grad(2), tv(2)
!	real*8, dimension(2):: ts
	integer:: k, time, i, j, eval, dim, ii

!	Here new definition for LCC estimation
	real*8:: x1r, x2r, tpr, minur, ddr, ttr, JJJr
	real*8:: rr(100), cir
	real*8, dimension(n):: ur
	integer:: kr

	! This part will alternatively estimate the alpha and m value of each type of pipelines
	! Thus, it is require to seperate the section of caculation in following steps

	! first to read the data from input file
	alpha=3.44E-05
	m=2.144000053
	
!	alpha= 7.46E-06
!	m= 2.335

	!---------------------------------------------
	! From here, we use subroutine to estimate the life cycle cost and optimal renewal time

!-----------------------------------------------------------------------------------
!	Defining the social cost, direct cost and the interest rate

!	Range for discount rate
	r=0.04 ! (4%) this is interest rate (Benchmark)
!	r=0.03 !(3%)
!	r=0.05 ! 5%
 !	r=0.0000000000000000000000000000000005
	! this step is for type A
	c1A=5000000 ! this is social cost (Benchmark)
!	c1A=6000000 ! this is social cost (upper bound)
!	c1A=4000000 ! this is social cost (lower bound)
	
	
!	c2A=500000 ! this is direct repairing cost (Benchmark)
!	c2A=1500000 ! this is direct repairing cost (upper)
	c2A=1000000 ! this is direct repairing cost (lower)

	call step1(alpha,m,c1A,c2A,r)
	! for type C
	alpha1=1.24E-05
	m1=2.484087137

	
	! for type F
!	alpha1=2.72E-05
!	m1=2.288024872

!	! for type FL
!	alpha1=1.93E-05
!	m1=2.390954544


	open(23, file='output\output.csv')


	
	c1C=c1A ! this is social cost
	c2C=c2A ! this is direct repairing cost

	! if we keep the value of LCC change according to years
	call step2(alpha1,m1,c1C,c2A,r)
	
!	lcc2=lcc
!	print*, lcc2
	

!	If we keep same value of LCC for estimation
!	call step3(alpha1,m1,c1C,c2A,r,lcc2)

	call step1(alpha1,m1,c1C,c2C,r)



	End program optimalrenewaltime
*******************************************************************	
!	Pattern method to estimate the value of m and alpha (or ts(1) and ts(2)



*****************************************************************************
!	STEP 1
!	Subroutine for estimating the life cycle cost and optimal value of renewal time
! The assmuption is to use the same type of pipe for renewal the system

	subroutine step1(alpha,m,c1,c2,r)
	use variable; implicit none
!	integer, parameter:: n=3*10**4
	real*8:: alpha, m, r
	real*8:: c1, c2 ! c1 is social cost, c2 is direct repair cost
	real*8:: x1, x2, tp, minu, ci, dd, tt, JJ
	real*8, dimension(n):: u
	integer:: i, k
	
	open(23, file='output\output.csv')

	! Defining the ininitual value

	x1=0.0
	x2=1.0
	tp=0.0
	k=0.0
	tt=0.01
	dd=0.01

	do i=1, n
			x1=x2
			x2=exp(-alpha*tt**m-r*tt)
			tp=tp+(x1+x2)*dd/2
			JJ=(c1+c2-c1*x2)/(r*tp) - (c1+c2)
			u(i) =JJ
			tt=tt+dd
!		Print*, u(i)

!		write(1,*) u(i)
	write(23,*) u(i)
		

	end do
	! Finding Minimin of LCC and the respective number of year
	
	Minu=u(1)
	do i=1, n
	if (Minu.GT.u(i)) Minu=U(i)
	end do
	
	lcc=Minu

	do i=1, n
	if (U(i).EQ.Minu) k=i
	

!	Print*, k
	end do

	Print*, "Mininimum LCC will be", Minu
	Print*, "Z vaule", k/100
	
	close(23)

	end subroutine

*************************************************************************************************
!	STEP 2
!	Subroutine for estimating the life cycle cost and optimal value of renewal time
!	The assmuption is to use the new type of pipe for renewal the old pipe 

	subroutine step2(alpha,m,c1,c2,r)
	use variable; implicit none

	real*8:: alpha, m, r
	real*8:: c1, c2 ! c1 is social cost, c2 is direct repair cost
	real*8:: x1, x2, tp, minu, ci, dd, tt, JJJ
	real*8, dimension(n):: u
	integer:: i, k
	

	open(24, file='output\newlcc.csv')
	! Defining the ininitual value

	x1=0.0
	x2=1.0
	tp=0.0
	k=0.0
	tt=0.01
	dd=0.01

	do i=1, n

			read(23,*) lcc1(i)
	
				x1=x2
			x2=exp(-alpha*tt**m-r*tt)
			tp=tp+(x1+x2)*dd/2
			JJJ=(c1+c2+lcc1(i))*(1-x2-r*tp)+(c2+lcc1(i))*x2
			u(i) =JJJ
			tt=tt+dd
!		Print*, u(i)

!		write(1,*) u(i)
	write(24,*) u(i)
		

	end do
	! Finding Minimin of LCC and the respective number of year
	
		Print*, lcc1(1)

!		pRINT*, alpha, m
	
	Minu=u(1)
	do i=1, n
	if (Minu.GT.u(i)) Minu=U(i)
	end do
	
	do i=1, n
	if (U(i).EQ.Minu) k=i
	

!	Print*, k
	end do

	Print*, "NEW LCC will be", Minu
	Print*, "Z vaule", k/100
	
	close (23)
	close(24)

	end subroutine

*************************************************************************

!Note: If we use the subroutine Newton for estimating the value of alpha and m, it is possible but it has limitation that 
! the inverse maxtrix might not be estimated due to the unavailability of matrix rank after some iterations. Thus, it is 
! better to combine the pattern method and Newton Rapson as supplementaries for each other. Pattern method is used to 
! estimate alpha and m, then we use this result as the input for Newton Method just to estimate the t-value. 

*************************************************************************************************
!
