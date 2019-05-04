*******************************************************
! This program is for estimating the parameter of Weibull hazard function
! and estimating the life cycle cost, optimal year of renewal for
! pipelines system
! Paper name: Optimal Renewal of Water Supply Pipeline
! Written and programed by Nam, Le Thanh on April, 30th 2008
! Re-programed and mofified on Nov, 23th, 2008
*****************************************************************
	Module variable
	integer, parameter:: totaldata=548
	real(8), dimension(totaldata):: d, t, le, x0, x1, x2
!	real(8), dimension(totaldata):: d, t, lcc
	integer, parameter:: timemax=1000000 !000000
	
	character(len=19), parameter:: file1="input\FL.csv" !
!	the input file must be checked carefully for the exact No. of data


	character(len=19), parameter:: file2="output\result.text"

	End module
		
***********************************************************************
! Main program starts here
***********************************************************************	
	Program optimalrenewaltime
	use variable; implicit none
	real*8:: m, eps, obj, c1, c2, r, AIC, RSS, alpha0, alpha1
	real*8:: alpha, hess(2,2), ts(2), optts(2), df(2), grad(2), tv(2)
!	real*8, dimension(2):: ts
	integer:: k, time, i, j, eval, dim


	open(10, file=file1)
	do k=1, totaldata
		read(10,*) t(k), d(k), x0(k), x1(k), x2(k)
	end do
	

	! This part will alternatively estimate the alpha and m value of each type of pipelines
	! Thus, it is require to seperate the section of caculation in following steps

	
	! first to read the data from input file
	!-----------------------------
	!TYPE C
	! for length
!	alpha0=2.51E-06
!	alpha1=1.94E-04
	
!	For rut
	!alpha0=1.38E-05
	!alpha1=-3.92E-06
!	m=2.484
	!-------------------------------

	!-----------------------------
	!TYPE F
	! for length
!	alpha0=4.92E-06
!	alpha1=3.25E-04
	
!	For rut
!	alpha0=2.56E-05
!	alpha1=6.24E-07
!	m=2.288
	!-------------------------------
	!TYPE FL
	! for length
	alpha0=6.73E-06
	alpha1=1.22E-04
	
!	For rut
!	alpha0=4.80E-06
!	alpha1=2.16E-05
	m=2.391
	!-------------------------------

	!TYPE A
	! for length
!	alpha0=8.27E-06
!	alpha1=4.18E-04
	
!	For rut
!	alpha0=2.56E-05
!	alpha1=6.24E-07
!	m=2.144
	!-------------------------------


	!---------------------------------------------

!	AKaike information cretiria (AIC), estimating Residual Sum of Square

	rss=0.0
		do k=1, totaldata
		rss=rss+ (1-d(k))*((alpha0*x0(k)+alpha1*x1(k))*t(k)**m) 
     #		+ d(k)*(log(m) + (m-1.0)*log(t(k))+
     #	log(alpha0*x0(k)+alpha1*x1(k)))
     #	- (alpha0*x0(k)+alpha1*x1(k))
     #      *t(k)**m
		end do
	Print*, "Value of RSS", Rss
	close(10)

	! From here, we use subroutine to estimate the life cycle cost and optimal renewal time


	End program optimalrenewaltime
*******************************************************************	
!	Pattern method to estimate the value of m and alpha (or ts(1) and ts(2)

	subroutine pattern(ts,dim)
	integer,intent(in)::dim
	real(8),dimension(dim+1)::MAXOBJ
	real(8),dimension(dim)::ts,b,dels
	real(8),dimension(dim+1,dim)::tt,del

	real(8)::OBJ,eps,op1,op2,op3 
	integer::i,j

	!****Setting of width of search****
		MAXOBJ=-10
		MAXOBJ(1)=10000
		dels(1)=2.5D-7
		dels(2)=2.5D-7


	!*******************
	eps=0.00000001
	b=ts

	do j=1,1000000000

	!	print*,j
	!	print*,dels
		if(dels(1)>eps)then

			tt(1,:)=ts
			del=0
			do i=1,dim
				del(i+1,i)=dels(i)
			end do


			do i=1,dim
				op1=OBJ(tt(i,:))
				op2=OBJ(tt(i,:)+del(i+1,:))
				op3=OBJ(tt(i,:)-del(i+1,:))


				if(op2>op1.and.op2>op3)then
					tt(i+1,:)=tt(i,:)+del(i+1,:)
				else if(op3>op1.and.op3>op2)then
					tt(i+1,:)=tt(i,:)-del(i+1,:)
				else if(op1>=op3.and.op1>=op2)then
					tt(i+1,:)=tt(i,:)
				end if
			end do


				if(OBJ(tt(dim+1,:))>OBJ(b))then
					ts=2.0*tt(dim+1,:)-b
					b=tt(dim+1,:)
				else
					ts=b
					dels=1.0/2.0*dels
				end if

		else 
			exit
		end if
		if(OBJ(b)>MAXOBJ(1))then
			MAXOBJ(1)=OBJ(b)
			do i=1,dim
				MAXOBJ(i+1)=b(i)
			end do
		end if
	print*, "No. of iteration", j
	print*,b,OBJ(b)
	end do

	end subroutine pattern

************************************************************************************
!	Subroutine Newton Raphson to estimate the unknown value of likelihood funtion
!	This subroutine only takes alpha and m into account (do not consider yet the length variable)

	Subroutine Newton(ts,hess,df)
	use variable; implicit none

	real*8, dimension(2):: ts, df, grad, hs, tv, ms, ss, mms, c, optts
	real*8, dimension(2,2):: hes, hess 

	real*8:: alpha, m, eps, m0, alpha0, maxS, ts10, ts20
	integer:: i,j,k, time

	eps=0.0000001

	ts10=ts(1)
	ts20=ts(2)


***********************************************************	
	do 400 time=1, 1

	print*, "No. of iteration", time

!	Asumming the initial value 
	df=0.0
	hes=0.0
	hess=0.0
	grad=0.0

	! here ts(1) = alpha and ts(2)=m

	do  i=1,2
		do k=1, totaldata
	! Formula for the first derivative of Objective likihood function
!------------------------------------------------
		df(1)=df(1) + d(k)/ts(1) - t(k)**ts(2)
		!with respect to alpha
		df(2)=df(2)+ d(k)*(1.0/ts(2)+
     #		log(t(k)))-ts(1)*ts(2)*(t(k)**(ts(2)-1))
	    ! with respect to m
	! Formula for the second  derivative of Objective likihood function
	
	hes(1,1)=hes(1,1)-d(k)/(ts(1)**2)
	hes(1,2)=hes(1,2)-ts(2)*(t(k)**(ts(2)-1))

	hes(2,1)=hes(2,1)-ts(2)*(t(k)**(ts(2)-1))
	hes(2,2)=hes(2,2)-d(k)/(ts(2)**2)-ts(1)*(t(k)**(ts(2)-1)
     #	+ts(2)*(ts(2)-1)* (t(k)**(ts(2)-2)))
!-------------------------------------------------		
		end do
	end do

	grad=df
	hess=hes

!	do i=1, 2
!		print*, (hess(i,j), j=1,2)
!	end do
!	call gaussj(hess,2,grad)
!	do i=1, 2
!	c(i)=grad(i)
!	end do
!	do i=1, 2
!		print*, c(i)
!	end do	
!	stop

! This is another way to estimate the value df(i)/hess(i,j) = c(i)
! in this case,if we do not use gaussj, and use matrix_inverse instead
! the c(i) = hs(i)


	call matrix_inverse(hess,2)

	hs=0.0
	
	do i=1, 2
		do j=1, 2
			hs(i)=hs(i)+hess(i,j)*df(j)
		end do
	end do
!	do i=1, 2
!		print*, hs(i)
!	end do

	ss(1)=ts(1)
	ss(2)=ts(2)
	ms=ss-hs
	mms=ms-ss
	
	do i=1,2
		mms(i)=mms(i)/ms(i)
	end do	

!	maxs=maxval(abs(mms))

!	if (maxs.lt.eps) then
!	exit
!	else
!	ts(1)=ms(1)
!	ts(2)=ms(2)
!	end if
!	if(ts(2)<0.or.ts(2)>20.or.ts(1)>700) then
	ts(1)=ts10
	ts(2)=ts20
!	end if



	Print*, "Value of alpha and m"
	Print*, "alpha=", ts(1)
	Print*, "m=", ts(2)

!	stop

400	continue

	Print*, "hessian matrix"
	
	do i=1, 2
		Print*, (hess(i,j), j=1,2)
	end do

	

	tv(1)=ts(1)/(sqrt(-hess(1,1)))
	tv(2)=ts(2)/(sqrt(-hess(2,2)))

	Print*, "T-value"

	Print*, "T-value for alpha", tv(1)
	Print*, "T-value for m", tv(2)
	open(12, file='output\t-value.text')
	Write(12,*) "=====t-value========="
	write(12,*) "t- for alpha=", tv(1)
	write(12,*) "===t- for  m=", tv(2)




	End subroutine Newton

*******************************************************************
! Subroutine for invert matrix and matrix multiply based on Gaussian Estimation
! Chuong trinh tinh ma tran nghich dao cua A
c
	SUBROUTINE gaussj(a,n,b)
	INTEGER n
	integer::NMAX
	REAL*8 a(n,n),b(n)
	PARAMETER (NMAX=50)	
	INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
	REAL big,dum,pivinv
	do j=1,n
	  ipiv(j)=0
	end do
	do i=1,n
	 big=0.
	 do j=1,n
	  if(ipiv(j)/=1) then
      do k=1,n
        if (ipiv(k)==0) then
          if (abs(a(j,k))>=big) then
            big=abs(a(j,k))
            irow=j
            icol=k
          endif
        else if (ipiv(k)>1) then
          pause 'singular matrix in gaussj'
        endif
      end do
	  endif
	end do
	ipiv(icol)=ipiv(icol)+1
	if (irow/=icol) then
	  do l=1,n
      dum=a(irow,l)
      a(irow,l)=a(icol,l)
      a(icol,l)=dum
	   end do
	   dum=b(irow)
	b(irow)=b(icol)
	b(icol)=dum
	 endif
	indxr(i)=irow
	indxc(i)=icol
	 if (a(icol,icol)==0.) pause 'singular matrix in gaussj'
	 pivinv=1./a(icol,icol)
	 a(icol,icol)=1.
	 do l=1,n
	  a(icol,l)=a(icol,l)*pivinv
	 end do
	 b(icol)=b(icol)*pivinv
	do ll=1,n
	   if(ll/=icol) then
      dum=a(ll,icol)
      a(ll,icol)=0.
      do l=1,n
        a(ll,l)=a(ll,l)-a(icol,l)*dum
      end do
      b(ll)=b(ll)-b(icol)*dum
	  endif
	end do
	end do
	do l=n,1,-1
	if(indxr(l)/=indxc(l)) then
	do k=1,n
		dum=a(k,indxr(l))
		a(k,indxr(l))=a(k,indxc(l))
		a(k,indxc(l))=dum
	end do
	endif
	end do
	END SUBROUTINE gaussj
****************************************************************************

*******************************************************************************
!	Function for Likelihood function OBJ without considering the length of pipelines
	Real(8) function obj(x)
	use variable
	real(8), dimension(2):: x
	real(8):: m, alpha, ts(2)

!	Here, we change all alpha and m value into ts


	alpha=x(1)
	m=x(2)
!	ts(1)=x(1)
!	ts(2)=x(2)
	obj=0

	do k=1, totaldata
		obj=obj+(1-d(k))*(-alpha)*t(k)**m  
     #	+ d(k)*(log(alpha)+log(m)+(m-1.0)
     #    *log(t(k))-alpha*t(k)**m)
	end do
	end function obj

*****************************************************************************
!!!!!!! Subroutine tinh Matran nghich dao cua A voi N la kich thuoc cua matran
!!!!!!! This subroutine estimates the inversve matrix of A with N is the dimension of the Matrix
	SUBROUTINE MATRIX_INVERSE(A,N)
	IMPLICIT NONE
	INTEGER,INTENT(IN):: N
	REAL*8:: A(N,N),AW,MM(N)
	INTEGER:: K,I,J,P,Q
	REAL*8 EPS
	INTEGER LW,MW,L(N),M(N)

	EPS=1.0D-10

	DO J=1,N
		L(J)=J;M(J)=J
	END DO

	DO 100 K=1,N
		AW=DABS(A(K,K))
		P=K;Q=K
		DO J=K,N; DO I=K,N
			IF(AW < DABS(A(I,J))) THEN
				AW=DABS(A(I,J)); P=I; Q=J
			END IF
		END DO; END DO

		IF(AW < EPS) THEN
			PRINT *, "This is abnormal" ;STOP
		END IF

		IF(K /= P) THEN
			DO J=1,N
				AW=A(K,J)
				A(K,J)=A(P,J)
				A(P,J)=AW
			END DO
			MW=M(K);M(K)=M(P);M(P)=MW
		END IF

		IF(K /= Q) THEN
			DO I=1,N
				AW=A(I,K)
				A(I,K)=A(I,Q)
				A(I,Q)=AW
			END DO
			LW=L(K);L(K)=L(Q);L(Q)=LW
		END IF

		AW=A(K,K);	A(K,K)=1.0
		DO J=1,N
			A(K,J)=A(K,J)/AW
		END DO

		DO I=1,N
			IF(I /= K) THEN
				AW=A(I,K);A(I,K)=0.0
				DO J=1,N
					A(I,J)=A(I,J)-AW*A(K,J)
				END DO
			END IF
		END DO
 100	CONTINUE

	DO J=1,N
		DO I=1,N
			P=L(I)
			MM(P)=A(I,J)
		END DO
		DO I=1,N
			A(I,J)=MM(I)
		END DO
	END DO

	DO I=1,N
		DO J=1,N
			Q=M(J)
			MM(Q)=A(I,J)
		END DO
		DO J=1,N
			A(I,J)=MM(J)
		END DO
	END DO
	
	RETURN
	
	END SUBROUTINE	MATRIX_INVERSE

*****************************************************************************
	! STEP 1
!	Subroutine for estimating the life cycle cost and optimal value of renewal time

	subroutine step1(alpha,m,c1,c2,r)
	implicit none
	integer, parameter:: n=3*10**4
	real*8:: alpha, m, r
	real*8:: c1, c2 ! c1 is social cost, c2 is direct repair cost
	real*8:: x1, x2, tp, lcc, minu, ci, dd, tt, JJ
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
			tp=tp+(x1+x2)*dd/2.0
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
	
	lcc =minu
	

	do i=1, n
	if (U(i).EQ.Minu) k=i
	

!	Print*, k
	end do

	Print*, "Mininimum LCC will be", Minu
	Print*, "Z vaule", k/100
	
	open(13, file='output\lcc.text')
	Write(13,*) "=====LCC========="
	write(13,*) "Minimum LCC will be=", Minu
	write(13,*) "Optimal renewal year", k/100


	end subroutine

*************************************************************************************************
	!	STEP 2
!	Subroutine for estimating the life cycle cost and optimal value of renewal time
!	The assmuption is to use the new type of pipe for renewal the old pipe 

	subroutine step2(alpha,m,c1,c2,r,lcc1)
	implicit none
	integer, parameter:: n=3*10**4
	real*8:: alpha, m, r
	real*8:: c1, c2, lcc1 ! c1 is social cost, c2 is direct repair cost
	real*8:: x1, x2, tp, lcc, minu, ci, dd, tt, JJJ
	real*8, dimension(n):: u
	integer:: i, k
	
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
			tp=tp+(x1+x2)*dd/2.0
			JJJ=(c1+c2+lcc1)*(1-x2-r*tp)+(c2+lcc1)*x2
			u(i) =JJJ
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
	
	do i=1, n
	if (U(i).EQ.Minu) k=i
	

!	Print*, k
	end do

	Print*, "Mininimum LCC will be", Minu
	Print*, "Z vaule", k/100
	

	end subroutine




*************************************************************************************************
!Note: If we use the subroutine Newton for estimating the value of alpha and m, it is possible but it has limitation that 
! the inverse maxtrix might not be estimated due to the unavailability of matrix rank after some iterations. Thus, it is 
! better to combine the pattern method and Newton Rapson as supplementaries for each other. Pattern method is used to 
! estimate alpha and m, then we use this result as the input for Newton Method just to estimate the t-value. 