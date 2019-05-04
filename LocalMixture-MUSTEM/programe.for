
	PROGRAM estimate hazard and transition
	IMPLICIT NONE
	
	INTEGER,PARAMETER:: JMAX=5  !No.of condition states
	INTEGER,PARAMETER:: noepsilon=8  !Total number of sample

	INTEGER:: i, j
	REAL*8:: theta(jmax-1), optheta(jmax-1,noepsilon)
	REAL*8:: epk(noepsilon)

	!------------------------------------------
	CHARACTER(LEN=30),PARAMETER:: FILE1="epsilon.csv"
	CHARACTER(LEN=30),PARAMETER:: FILE2="hazard.csv"
	CHARACTER(LEN=30),PARAMETER:: FILE3="grouphazard.csv"
	
	OPEN(10,file=FILE1)
	OPEN(11,file=file2)
	Do i=1, jmax-1
	read(11,*) theta(i)
	end do
	Print*, "value of theta, everage hazard value"
	Do i=1, jmax-1
	print*, theta(i)
	end do

	Do j=1, noepsilon
	read(10,*) epk(j)
	end do
	Print*, "value of heterogeneity, groups"
	Do j=1, noepsilon
	print*, epk(j)
	end do
	do j=1, noepsilon
		do i=1, jmax-1
		optheta(i,j)=theta(i)*epk(j)
		end do
	end do
	
	Print*, "value of hazard rate according to each group"
	do j=1, noepsilon
		Print*, "value of hazard rate of group k=1", j
		
		print*, (optheta(i,j), i=1, jmax-1)
	end do
	OPEN(12,file=FILE3)
	
	do j=1, noepsilon
		write(12,20) (optheta(i,j), i=1, jmax-1)
	end do
	
	
20	FORMAT(1X,6(F10.4,","),F10.4)
	
	CLOSE(10)
	CLOSE(11)
	END PROGRAM 
	