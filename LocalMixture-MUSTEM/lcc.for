************************************************************************
	PROGRAM VARIABLE
		INTEGER,PARAMETER:: JMAX=5  !No.of condition states
		Integer, parameter:: ekmax =5 ! No of repair action
		Integer, parameter:: omax =3 ! No of repair action
		INTEGER,PARAMETER:: Year=20

	CHARACTER(LEN=19),PARAMETER:: FILE1="input\cost.csv"
	CHARACTER(LEN=19),PARAMETER:: FILE2="output\prob.csv"
	CHARACTER(LEN=19),PARAMETER:: FILE3="output\lcc.csv"
	CHARACTER(LEN=19),PARAMETER:: FILE4="output\sumlcc.csv"
	CHARACTER(LEN=19),PARAMETER:: FILE5="output\minlcc.csv"

	! Value of d can be changed upon the life expectancy of the ratings
	REAL*8:: R(jmax,jmax), Q(jmax,jmax), P(jmax,jmax), opd(year-1)
	REAL*8:: E(ekmax,Jmax), C(ekmax,Jmax), CK(Jmax) , LCC(Jmax-1,Year)
	REAL*8:: SUMLCC(ekmax,year), MINLCC(year), A(ekmax-1,year)

	REAL*8:: L(Jmax-1,Year), LIFE
	REAL*8,PARAMETER:: rate = 0.9
	INTEGER I,J, t, ek, o,m, sigma(ekmax,jmax), k, x,y,z, d

	open(10,file=file1)
	Open(20,file=file2)
	Open(30,file=file3)
	Open(40,file=file4)
	Open(50,file=file5)

	do i=1, jmax
		read(20,*) (P(i,j), j=1, jmax)
	end do

	! Repair matrix (Preventive action), number of preventive rule should be less than jmax

	do 44 ek=3, ekmax

	Do I=1, Jmax
		Do J=1, Jmax
			If (i==j .and. i < ek) Then 
				R(i,j)=1 
			Else
				R(i,j)=0 
				If (i/=j .and. i>=ek) Then
				R(i,1)=1

			End if; end if
		End Do
	End Do

	Q=matmul(P,R)

	E=0
	read(10,*) (c(ek,j), j=1, jmax)
			do i=1,Jmax-1
				do j=i, jmax ! 
					E(ek,i)=E(ek,i)+P(i,j)*C(ek,j) ! E_i
				end do
			end do
	! Life cycle cost estimation
	
	L=0
	LCC=0
	
	!SUMLCC(ek,year)=0
	
	do i=1, jmax-1
	L(i,year)=0
	LCC(i,year)=0
	end do
		do 45 t=year, 2, -1
			do i=1, jmax-1
				do	j=1, jmax-1
		
					L(i,t) = L(i,t)+Q(i,j)*LCC(j,t)
					LCC(i,t-1)=rate*(E(ek,i)+L(i,t))
				end do
			
	!		SUMLCC(ek,t)=SUMLCC(ek,t)+LCC(i,t)
			
			end do
	
45	continue
	
	Print*, "Value of matrix E(ek,i)"
	Print11, (E(ek,i), i=1, jmax-1)

	Print*, "Value of LCC with respect to condition j=" ,ek, 
     #"i=1"
		
	do i=1, jmax-1
	Print 11, (LCC(i,t), t=1,year)
	WRITE(30,11) (LCC(i,t), t=1,year)
	end do

******************************************************************
	! Caculating the SUM of repair cost in each year and make comparison
	
	SUMLCC=0
	do i=1, jmax-1
		do t=1, year
			SUMLCC(ek,t) = SUMLCC(ek,t)+LCC(i,t)
		end do
	end do
	

	!PRINT*, "SUM of repairing cost respect to each year"
	!PRINT 11, (SUMLCC(ek,t),t=1, year)
	WRITE(40,11) (SUMLCC(ek,t),t=1, year)
	
44	Continue

	CLOSE(40)
	Open(40,file=file4)

c	Finding the minimum value of LCC coressponding each year of calculation
	!Stop

	PRINT*, "Sumation of repairing cost with respect to each repair 
     #	strategy and year"

	do i=1, ekmax-2
		Read(40,*) (A(i,t), t=1, year)
		PRINT 11, (A(i,t), t=1, year)
	end do
	!Stop

	Do 47 t=1, year-1
		MINLCC(t)=A(1,t)
		do 48 i=1, ekmax-2
			if (MINLCC(t).GT.A(i,t)) MINLCC(t)=A(i,t)
	
48		Continue
	!	WRITE(50,11) MINLCC(t)

	PRINT*, "Minimum LCC value of year=",t ,"=" 
	 
	PRINT 11, MINLCC(t)
	
	do i=1, ekmax-2
		if (MINLCC(t)==A(i,t)) d=i
		end do
	opd(t)=d
	
	PRINT*, "Best repair action will be d =", d
	!WRITE(50,11) opd(t)

47	Continue
	
	WRITE(50,*) "********Minimum LCC and best repair action*******"
	do t=1, year-1
	WRITE(50,11) MINLCC(t), opd(t)
	end do


11	FORMAT(1X,6(F10.4,","),F10.4)
	stop
	CLOSE(10)
	CLOSE(20)	
	CLOSE(30)
	CLOSE(40)
	CLOSE(50)

	End 
**************************************************************
c	Result of program

