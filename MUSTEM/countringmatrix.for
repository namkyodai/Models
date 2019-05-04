

	PROGRAM countingmatrix
	IMPLICIT NONE
	
	INTEGER,PARAMETER:: JMAX=5  !No.of condition states
	INTEGER,PARAMETER:: Totalnumber=2064  !Total number of sample
	!------------------------------------------
	INTEGER,PARAMETER:: TOTALDATA=2064	!No. of total data(for inverse)
	INTEGER:: i, j, CI, CJ
	INTEGER SCOUNT(JMAX,JMAX)
	REAL PROB(JMAX,JMAX)
	REAL A(JMAX)

	!------------------------------------------
	CHARACTER(LEN=30),PARAMETER:: FILE1="input\inputdata.csv"
	
	CHARACTER CH! !Useless character
	
	REAL*4, DIMENSION(Jmax,Jmax):: P
	
	OPEN(10,file=FILE1)
	READ(10,*) CH
	SCOUNT=0
	
	DO I=1,Totalnumber
		READ(10,*) CI, CJ
		SCOUNT(CI,CJ)=SCOUNT(CI,CJ)+1
	END DO

	!OPEN(11,file=FILE2)

	DO I=1, JMAX
	!	WRITE(11,*) (SCOUNT(I,J), J=1, JMAX)
		PRINT*, (SCOUNT(I,J),J=1,JMAX)
		A=SUM(SCOUNT, dim=2) ! This command will calculate the sum of row
	END DO
				
	DO I=1,JMAX
		PRINT*, A(I)
	END DO
		DO I=1,JMAX
			DO J=1, JMAX			
			IF (SCOUNT(I,J)==0) THEN
				PROB(I,J) = 0
			ELSE
				PROB(I,J)=real(SCOUNT(I,J))/A(I)
			END IF
		
		END DO
	END DO
	
	DO I=1,JMAX
	PROB(JMAX, JMAX) = 1
	PRINT*, (PROB(I,J), J=1, JMAX)
	END DO


	CLOSE(10)
	CLOSE(11)
	END PROGRAM countingmatrix
	