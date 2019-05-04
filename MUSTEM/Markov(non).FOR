************************************************************************
*	Estimation of Markov transition probabilities
************************************************************************
	MODULE VARIABLE
		INTEGER,PARAMETER:: JMAX=5  !No.of condition states
		INTEGER,PARAMETER:: MMAX=3  !No.of caracteridtic variable
		INTEGER,PARAMETER:: d=5  
		INTEGER,PARAMETER:: Year=8

	!------------------------------------------
		INTEGER,PARAMETER:: TOTALDATA=1237	!No. of total data(for inverse)

	!------------------------------------------
		INTEGER,PARAMETER:: TIMEMAX=100	!repeat No

		INTEGER,PARAMETER:: dimN=(JMAX-1)*MMAX !No. of simultaneous equations

	!---the case you want to change the No. of parameters-------------
		INTEGER,PARAMETER:: D1=6	!total No. of deleted parameters
		INTEGER,PARAMETER:: redimN=dimN-D1	!new No. of parameters
		INTEGER DEL(D1)
		DATA DEL/2,3,7,9,10,12/!No. of parameters to be deleted((J-1)*MMAX+M)
		!handling initial value of parameters (BETA(J)=J1*J+B1)

	
		!REAL*8,PARAMETER:: B1=1.3D-05,J1=5.1D-05	!30inv,15(8,11deleted)○
		REAL*8,PARAMETER:: B1=1.3D-04,J1=5.1D-04	!30inv,15(8,11deleted)○

		!-----------------------------------------
		CHARACTER(LEN=19),PARAMETER:: FILE1="input\inputdata.csv"
			!No. of total data	 MAX26748
		CHARACTER(LEN=19),PARAMETER:: FILE2="input\cost.csv"
		CHARACTER(LEN=13),PARAMETER:: FILE1000="output\df.csv"
		CHARACTER(LEN=17),PARAMETER:: FILE1001="output\beta-t.csv"
		CHARACTER(LEN=17),PARAMETER:: FILE1002="output\hazard.csv"
		CHARACTER(LEN=17),PARAMETER:: FILE1003="output\avlife.csv"
		CHARACTER(LEN=17),PARAMETER:: FILE1004="output\repair.csv"
		CHARACTER(LEN=17),PARAMETER:: FILE1005="output\prob.csv"
		
		REAL*8:: THETA(JMAX) !hazard function Θj
		REAL*8:: THETASA(JMAX,JMAX)	!Θj-Θi
		INTEGER:: MAXTHETA
		CHARACTER CH !useless character

		REAL*8 X1(MMAX)
		REAL Z1
		DATA X1/1.0, 0.337904, 0.133927/	!average
		DATA Z1/2/
	END MODULE
*************************
*	Main program
*************************
	PROGRAM ESTIMATION
	USE VARIABLE; IMPLICIT NONE
	REAL*8,EXTERNAL:: PROB,OBJ
	
	REAL*8:: OPTBETA(MMAX,JMAX-1),BETA(MMAX,JMAX-1)
	REAL*8:: G(dimN,dimN),A(redimN,redimN)
	REAL*8:: OPTOBJ,TESTOBJ,TOKEI
	REAL*8,PARAMETER:: EPS1=1.0D-06
	REAL TV(MMAX)
	REAL Z
	INTEGER I,J,L,M,TIME,EVAL

*------------------------------------------
*	Maximam likelihood estimation method
*------------------------------------------
!-----initial value of parameter β------------------------------
	DO J=1,JMAX-1;	DO M=1,MMAX
		BETA(M,J)=J*J1+B1	!better. 
				
	END DO;	END DO
	DO I=1,D1
		M=MOD(DEL(I)+(MMAX-1),MMAX)+1;J=(DEL(I)-M)/MMAX+1
	! Gia tri moi cua M va J se tuong ung voi vi tri ma gia tri beta(M,J) se bi xoa
		BETA(M,J)=0.0
	END DO

		!----------------------------------------------------
	OPTBETA=BETA ! Do gia tri beta khong thoa man da bi xoa nen co the doi
      !bien moi la OPTBETA, nhung co gia tri tuong duong duong voi Beta

	DO 100 TIME=1,TIMEMAX 
		!simultadneous nonlinear equations
		!determine approximate solution of β
		CALL NEWTON_RAPHSON(OPTBETA,G)
	
		EVAL=1	!EVAL=1 then it is convergent
		DO I=1,JMAX-1;DO M=1,MMAX
			!else, you have to recalcuration
		IF(DABS( (OPTBETA(M,I)-BETA(M,I))/OPTBETA(M,I) ) > EPS1) THEN
				EVAL=0;
		GOTO 200
			END IF
		END DO; END DO
 200		CONTINUE	!escape from IF sentence

		!renewal β
		BETA=OPTBETA
		IF(EVAL==1)	EXIT	!if β is convergent then escape
		
		IF(MOD(TIME,10)==0) THEN
		PRINT *, TIME ," times of loop"
		PRINT *, " Value of optimal betas"
			DO I=1,JMAX-1
				PRINT *, (BETA(M,I),M=1,MMAX)
			END DO
		END IF
		
		IF(TIME==TIMEMAX) THEN
			PRINT *, "not convergent, stop and increate the times of
     & looping";STOP
		END IF
 100	CONTINUE
	
	PRINT*, "Value of Hessian Matrix"
	Do i=1, redimN
	PRINT*, (G(i,j), j=1, redimN)
	End Do

*-------------------
*	Test
*-------------------
	CALL DEL_MATRIX(G,A)
	A=-A; CALL MATRIX_INVERSE(A,redimN)	!inverse matrix of Fisher Matrix
!	test
	OPTOBJ=OBJ(OPTBETA)
	PRINT *, "likelihood function"
	PRINT *, OPTOBJ
	OPEN(1001,file=FILE1001)
	
	!parameter β
	PRINT *, "unknown parameters β********************"
	WRITE(1001,*) "unknown parameters β***************"
	DO I=1,JMAX-1
		PRINT *, (OPTBETA(M,I),M=1,MMAX)
		WRITE(1001,20) (OPTBETA(M,I),M=1,MMAX)
	END DO 
 
	PRINT *, "t-value********"
	WRITE(1001,*) "t-value********"
	L=0
	DO J=1,JMAX-1
		DO 210 M=1,MMAX
			DO I=1,D1
				IF ((J-1)*MMAX+M==DEL(I))  THEN
					TV(M)=0.0; GOTO 210
				END IF
			END DO
			L=L+1; TV(M)=OPTBETA(M,J)/DSQRT(A(L,L))
210		CONTINUE
		PRINT *, (TV(M),M=1,MMAX)
		WRITE(1001,20) (TV(M),M=1,MMAX)
	END DO 
	CLOSE(1001)
 20	FORMAT(1X,2(F10.4,","),F10.4)

	!hazard function・transition prob matrix・unknown parameter β
	Z=Z1
	CALL PRINT1(OPTBETA,Z)

	STOP
	END PROGRAM

*****************************
*	Newton-Raphson method
*****************************
	SUBROUTINE NEWTON_RAPHSON(BETA,HES)
	USE VARIABLE; IMPLICIT NONE

	REAL*8:: BETA(MMAX,JMAX-1)
	REAL*8:: G(dimN,dimN),HES(dimN,dimN),F(dimN),C(redimN), b(redimN)
	REAL*8:: A(redimN,redimN),X(redimN)
	REAL*8:: df(JMAX-1), df2(JMAX-1,JMAX-1)
	REAL*8:: XK(MMAX,JMAX-1)
	INTEGER:: IK,JK
	REAL::ZK
	INTEGER:: I,J,K,L,M,M2,N
	
	OPEN(1,file=FILE1)
	READ(1,*) CH

	F=0.0;	G=0.0
	DO 100 K=1,TOTALDATA
		!condition states，inspection interval，characteristic vector
		READ(1,*) IK,JK,ZK,(XK(M,1), M=1,MMAX)
		DO J=2,JMAX-1;DO M=1,MMAX
			XK(M,J)=XK(M,1)
		END DO; END DO
		DO I=1,D1
			M=MOD(DEL(I)+(MMAX-1),MMAX)+1;J=(DEL(I)-M)/MMAX+1
			XK(M,J)=0.0
		END DO
		
		!Θvalue
		CALL THETAVALUE(IK,JK,XK,BETA)
		
		!calculate dlnπ_(ij)/dΘ_l
		CALL FIRST_DIFFERENTIAL(IK,JK,ZK,df)

			!Tinh toan cac gia tri cua phuong trinh dao ham cap 1 first order partial differential equations
		DO L=IK,JK
			IF(L==JMAX) CYCLE
			DO M=1,MMAX
				F((L-1)*MMAX+M)=F((L-1)*MMAX+M)+df(L)*XK(M,L)
 			END DO
		END DO
	
		!Dao ham cap 2 second order partial differential equations
		CALL HESSIAN(IK,JK,ZK,df2,df) ! Tinh toan Ma tran Heessian

		!Hessian matrix
		DO L=IK,JK;	DO N=IK,JK; DO M=1,MMAX; DO M2=1,MMAX
			IF(L/=JMAX .AND. N/=JMAX) THEN
				G((L-1)*MMAX+M,(N-1)*MMAX+M2)
     &		=G((L-1)*MMAX+M,(N-1)*MMAX+M2)+df2(L,N)*XK(M,L)*XK(M2,N)
			END IF
		END DO;END DO;END DO;END DO

 100	CONTINUE

	OPEN(1000,file=FILE1000)	
	DO J=1,JMAX-1
		WRITE(1000,30) (F((J-1)*MMAX+M),M=1,MMAX)
	END DO
 30	FORMAT(1X,2(F14.6,","),F14.6)
	CLOSE(1000)

	!resolve GC=F
	!if you delete some parameters
	!you have to delete row and colum vector of them
	L=0
	DO 410 I=1,dimN
		DO J=1,D1
			IF(I==DEL(J)) GOTO 410
		END DO
		L=L+1; X(L)=F(I)
 410	CONTINUE

	CALL DEL_MATRIX(G,A)
	HES=G

	!CALL MATRIX_1(redimN,A,X,C)
	do i=1,redimN
		b(i)=X(i) 
	end do ! Dinh nghia dau vao C la gia tri X
	! Sau khi goi Gaussj thi gia tri C chinh la nghiem cua Ax=b

	Call gaussj(A,redimN,b)
	do i=1,redimN
	C(i)=b(i)
	end do	
	
	I=0
	DO J=1,JMAX-1
		DO 500 M=1,MMAX
			DO N=1,D1
				IF((J-1)*MMAX+M==DEL(N)) GOTO 500
			END DO
			I=I+1
			BETA(M,J)=BETA(M,J)-C(I)
 500		CONTINUE
	END DO


	CLOSE(1)
	END SUBROUTINE NEWTON_RAPHSON


*************************************************
*	subroutine for  df(L)=dlnπ_(IK,JK)/dΘ_L
*************************************************

	SUBROUTINE FIRST_DIFFERENTIAL(IK,JK,ZK,df)
	USE VARIABLE ;IMPLICIT NONE
	REAL*8,EXTERNAL:: PROD1
	INTEGER,INTENT(IN):: IK,JK
	REAL,INTENT(IN):: ZK
	REAL*8:: df(JMAX-1)
	INTEGER E,H,L,M,P
	REAL*8 SUM1,SUM2,RESERVE
	df=0.0
	DO 100 L=IK,JK
		IF(L==JMAX) CYCLE !THETA(JMAX) is not defined
		SUM2=0.0	!denominator

		DO H=IK,JK
			SUM1=0.0	!numelator
			IF(H==L) THEN
				DO P=IK,JK
					IF (P/=L) SUM1=SUM1+1.0/THETASA(P,L)
				END DO
				SUM1=SUM1-ZK
			ELSE
				SUM1=SUM1-1.0/THETASA(L,H)
			END IF
			IF(L /= JK) SUM1=SUM1+1.0/THETA(L)

			RESERVE=DEXP(-THETASA(H,MAXTHETA)*ZK)/PROD1(IK,JK,H)
			SUM1=SUM1*RESERVE
			SUM2=SUM2+RESERVE
			df(L)=df(L)+SUM1
		END DO

		df(L)=df(L)/SUM2
 100	CONTINUE
	RETURN
	END SUBROUTINE FIRST_DIFFERENTIAL

********************************************
*	subroutine for Hessian matrix
********************************************
	SUBROUTINE HESSIAN(IK,JK,ZK,HES,df)
	USE VARIABLE; IMPLICIT NONE
	REAL*8,EXTERNAL:: PROD1

	REAL*8 HES(JMAX-1,JMAX-1)
	INTEGER,INTENT(IN):: IK,JK
	REAL,INTENT(IN):: ZK
	REAL*8,INTENT(IN):: df(JMAX-1)
	INTEGER E,H,L,N,P
	REAL*8 SUM1,SUM2,SUM3,RESERVE

	HES=0.0

	!off-diagonal element
	DO 100 L=IK,JK-1
	DO 200 N=L+1,JK
		IF(N==JMAX) CYCLE	!THETA(JMAX)is not defined
		SUM2=0.0	!denominator
		
		DO H=IK,JK
			SUM1=0.0	!numerator
			IF(H==L) THEN
				DO P=IK,JK
					IF(P/=L) SUM1=SUM1+1.0/THETASA(P,L)
				END DO
				SUM1=SUM1+1.0/THETA(L)+1.0/THETASA(N,L)-ZK
				SUM1=-SUM1/THETASA(N,L)
			ELSE IF(H==N) THEN
				DO P=IK,JK
					IF(P/=N) SUM1=SUM1+1.0/THETASA(P,N)
				END DO
				SUM1=(SUM1-ZK)*(1.0/THETA(L)-1.0/THETASA(L,N))
				SUM1=SUM1-(1.0/THETASA(L,N)**2)
			ELSE
				SUM1=1.0/THETA(L)-1.0/THETASA(L,H)
				SUM1=-SUM1/THETASA(N,H)
			END IF
			RESERVE=DEXP(-THETASA(H,MAXTHETA)*ZK)/PROD1(IK,JK,H)
			HES(L,N)=HES(L,N)+SUM1*RESERVE
			SUM2=SUM2+RESERVE
		END DO

		HES(L,N)=-df(L)*df(N)+HES(L,N)/SUM2
		IF(N /= JK) HES(L,N)=HES(L,N)+df(L)/THETA(N)

		HES(N,L)=HES(L,N)	!symmetric matrix
 200	CONTINUE
 100	CONTINUE
	
	!diagonal element
	DO 300 L=IK,JK
		IF(L==JMAX) CYCLE	!THETA(JMAX)is not defined
		SUM2=0.0	!denominator
		DO H=IK,JK
			SUM1=0.0	!numerator
			DO P=IK,JK
				IF(P/=L) SUM1=SUM1+1.0/THETASA(P,L)**2
			END DO
			RESERVE=-ZK
			IF(H==L) THEN
				DO P=IK,JK
					IF(P/=L) RESERVE=RESERVE+1.0/THETASA(P,L)
				END DO
				SUM1=SUM1+RESERVE**2
				IF(L/=JK) THEN
					SUM1=SUM1+(RESERVE-1.0/THETA(L))/THETA(L)
				END IF
			ELSE
				SUM1=2.0/THETASA(L,H)**2
				IF(L/=JK) THEN
					SUM1=SUM1-(1.0/THETA(L)**2)
					SUM1=SUM1-1.0/(THETA(L)*THETASA(L,H))
				END IF
			END IF
				
			RESERVE=DEXP(-THETASA(H,MAXTHETA)*ZK)/PROD1(IK,JK,H)	
			HES(L,L)=HES(L,L)+SUM1*RESERVE
			SUM2=SUM2+RESERVE
		END DO
		HES(L,L)=-df(L)**2+HES(L,L)/SUM2
		IF(L /= JK) HES(L,L)=HES(L,L)+df(L)/THETA(L)
 300	END DO


	END SUBROUTINE HESSIAN

***********************
*	Θi(i=IK,JK)
***********************
	SUBROUTINE THETAVALUE(IK,JK,XK,BETA) ! Khong co lien quan gi den gia tri
     &! IK va JK
	USE VARIABLE; IMPLICIT NONE
	REAL*8,INTENT(IN):: BETA(MMAX,JMAX-1),XK(MMAX,JMAX-1)
	INTEGER,INTENT(IN):: IK,JK
	INTEGER I,J,M
	REAL*8 R
	!calculate Θ
	THETA=0.0
	MAXTHETA=IK;R=0.0
	DO J=IK,JK
		IF(J==JMAX) CYCLE
		DO M=1,MMAX
			THETA(J)=THETA(J)+BETA(M,J)*XK(M,J)
		END DO
	!check maximam Θ
		IF(R < THETA(J)) THEN		
			R=THETA(J);MAXTHETA=J
		END IF
	END DO
	!Θi-Θj
	THETASA=0.0
	DO J=IK,JK; DO I=IK,JK
		THETASA(I,J)=THETA(I)-THETA(J)
	END DO;	END DO
	END SUBROUTINE
**********************************
*	function
**********************************
*(1)  PROD1 \prod_{e=a,\ne c}^{b} \frac{1}{\bar{theta}_{e,c}}
	REAL*8 Function PROD1(aa,bb,cc)
	USE VARIABLE;IMPLICIT NONE
	INTEGER,INTENT(IN):: aa,bb,cc
	INTEGER E

	PROD1=1.0
	DO E=aa,bb
		IF(E /= cc) PROD1=THETASA(E,cc)*PROD1
	END DO

	RETURN;	END

***************************
*	likelihood function
***************************
	REAL*8 FUNCTION OBJ(BETA)
	USE VARIABLE;IMPLICIT NONE
	REAL*8,EXTERNAL:: PROB
	REAL*8,INTENT(IN):: BETA(MMAX,JMAX-1)
	INTEGER I,J,K,M
	INTEGER IK,JK
	REAL ZK
	REAL*8 XK(MMAX,JMAX-1)
	REAL*8 PI

	OPEN(1,file=FILE1)
	READ(1,*) CH
	
	OBJ=0.0
	DO K=1,TOTALDATA
		!read from file 
		!condition states，inspection interval，characteristic vector
		READ(1,*) IK,JK,ZK,(XK(M,1), M=1,MMAX)
		DO J=2,JMAX-1;DO M=1,MMAX
			XK(M,J)=XK(M,1)
		END DO; END DO
		DO I=1,D1
			M=MOD(DEL(I)+(MMAX-1),MMAX)+1;J=(DEL(I)-M)/MMAX+1
			XK(M,J)=0.0
		END DO
		
		!Den day ta co the biet duoc gia tri cua ma tran XK
	! Giu lieu tu file input bao gom IK, JK, XK va gia tri Beta
		!Θ
		CALL THETAVALUE(IK,JK,XK,BETA)! Goi rountine nay de tinh gia tri Theta
	! Dau ra la gia tri Theta

		PI=PROB(IK,JK,ZK) ! PI la gia tri xac xuat
		!likelihood
		IF(PI<=1.0D-305) THEN	!the case Θ<0
								!there are some case π<0
			OBJ=OBJ-702.2884
		ELSE
			OBJ=OBJ+DLOG(PI)
		END IF
	END DO
	PRINT*, "Value of OBJ", OBJ

	CLOSE(1)
	RETURN
	END FUNCTION

***************************************
*	transition probabilities π_ij(Z)
***************************************
	REAL*8 FUNCTION PROB(I,J,Z) ! Cach tinh gia tri ma tran bien doi PROB
	USE VARIABLE;IMPLICIT NONE
	REAL*8,EXTERNAL:: PROD1 ! Ham nay la ham nao?
	INTEGER,INTENT(IN):: I,J
	REAL,INTENT(IN):: Z

	REAL*8 RESERVE
	INTEGER K
	
	PROB=0.0
	IF(Z <= 1.0D-06) RETURN
	RESERVE=1.0
	DO K=I,J-1
		RESERVE=RESERVE*THETA(K)
	END DO
	DO K=I,J
		PROB=PROB+DEXP(-THETA(K)*Z)/PROD1(I,J,K)
	END DO
	PROB=RESERVE*PROB

	RETURN; END

****************************
*	subroutine for print
****************************
	SUBROUTINE PRINT1(BETA,Z)		!for hazard function & trans prob
	USE VARIABLE;IMPLICIT NONE
	REAL*8,EXTERNAL:: PROB

	REAL*8 BETA(MMAX,JMAX-1)
	INTEGER I,J,M, IK,JK, ZK
	REAL Z
	REAL*8 XK(MMAX,JMAX-1)
	REAL*8 df(JMAX-1)

	DO J=1,JMAX-1;	DO M=1,MMAX
		XK(M,J)=X1(M) ! This will effec the value of life expectancy
	END DO;END DO

	CALL THETAVALUE(1,JMAX,XK,BETA)	
			
	OPEN(1002,file=FILE1002)

	!Θ
	PRINT *, "hazard function****************************"
	WRITE(1002,*) "hazard function***********"
	PRINT*, "Values of Theta"
	PRINT 10, (THETA(J),J=1,JMAX)
	WRITE(1002,10) (THETA(J),J=1,JMAX)
	PRINT*, "Value of THETASA"
	DO J=1, JMAX
	PRINT*, (THETASA(I,J), I=1, JMAX)
	END DO

	!π_ij
	PRINT *, "transition probabilities********************************"
	WRITE(1002,*) "transition probabilities**************************"
	DO I=1,JMAX
		PRINT 10, (PROB(I,J,Z),J=1,JMAX)
		WRITE(1002,10) (PROB(I,J,Z),J=1,JMAX)
		WRITE(1005,10) (PROB(I,J,Z),J=1,JMAX)
	END DO

 10	FORMAT(1X,6(F10.4,","),F10.4)
	CLOSE(1002)
	
	PRINT *,"expected life expectancy"
	PRINT 10, (1.0/THETA(J),J=1,JMAX)
	!if you call it then θ changes
		CALL ELIFE(BETA)
	
	CALL repairmatrix
	

	END SUBROUTINE
*****************************************
*	averaged expected life expectancy
*****************************************
	SUBROUTINE ELIFE(BETA)
	USE VARIABLE;IMPLICIT NONE
	REAL*8,INTENT(IN):: BETA(MMAX,JMAX-1)
	INTEGER IK,JK
	INTEGER J,K,M
	REAL*8 XK(MMAX,JMAX-1)
	REAL*8 LIFE(JMAX-1)
	REAL ZK

	OPEN(1,file=FILE1)
	OPEN(1003,file=FILE1003)
	READ(1,*) CH
	
	LIFE=0.0
	DO K=1,TOTALDATA
		!condition states，inspection interval，characteristic vector
		READ(1,*) IK,JK,ZK,(XK(M,1), M=1,MMAX)
		DO J=2,JMAX-1;DO M=1,MMAX
			XK(M,J)=XK(M,1)
		END DO; END DO
		CALL THETAVALUE(1,JMAX,XK,BETA)
		DO J=1,JMAX-1
			LIFE(J)=LIFE(J)+1.0/THETA(J)
		END DO
	END DO
	LIFE=LIFE/TOTALDATA
	PRINT *, "averaged expected life expectancy"
	PRINT *, (LIFE(J),J=1,JMAX-1)
	WRITE(1003,10) (LIFE(J),J=1,JMAX-1)
  10	FORMAT(1X,6(F10.4,","),F10.4)

	CLOSE(1)
	RETURN;	END SUBROUTINE
***************************************************************
*	making matrix（delete row or colum that all factor are 0）
***************************************************************
	SUBROUTINE DEL_MATRIX(A,B)
	USE VARIABLE;IMPLICIT NONE
	REAL*8,INTENT(IN):: A(dimN,dimN)
	REAL*8 B(redimN,redimN)
	INTEGER I,J,L,N,K

	L=0
	DO 300 I=1,dimN
		DO K=1,D1
			IF(I==DEL(K)) GOTO 300
		END DO
		L=L+1;N=0
		DO 400 J=1,dimN
			DO K=1,D1
				IF(J==DEL(K)) GOTO 400
			END DO
			N=N+1;	B(L,N)=A(I,J)
 400		CONTINUE
 300	CONTINUE
	END SUBROUTINE

**************************************************************
*	solution of simulutaneous equations Ax=b	 N demension
**************************************************************
      SUBROUTINE MATRIX_1(N,A,b,X)
      REAL*8 A(N,N),W(N),X(N)
	REAL*8,INTENT(IN):: b(N)
      INTEGER IP(N)
      DATA EPS/1.0D-06/

      DO K=1,N; IP(K)=K; END DO
      DO 10 K=1,N
        L=K; AL=DABS(A(IP(L),K))
        DO I=K+1,N
            IF (DABS(A(IP(I),K)).GT.AL) THEN
                L=I; AL=DABS(A(IP(L),K))
            END IF
        END DO
        IF (L.NE.K) THEN
            LV=IP(K); IP(K)=IP(L); IP(L)=LV
        END IF
        IF (DABS(A(IP(K),K)).LE.EPS) GO TO 900
        A(IP(K),K)=1.0/A(IP(K),K)
        DO I=K+1,N
            A(IP(I),K)=A(IP(I),K)*A(IP(K),K)
        END DO
        DO 40  J=K+1,N
            DO I=K+1,N
              W(I)=A(IP(I),J)-A(IP(I),K)*A(IP(K),J)
            END DO
            DO I=K+1,N
                A(IP(I),J)=W(I)
            END DO
 40     CONTINUE
 10   CONTINUE
      GO TO 999
*
 900    CONTINUE
      PRINT *, 'this matrix is not normal';STOP
*
 999  DO J=1,N; X(J)=b(IP(J)); END DO
      DO 620 J=1,N
            DO I=J+1,N
                X(I)=X(I)-A(IP(I),J)*X(J)
            END DO
 620    CONTINUE
      X(N)=X(N)*A(IP(N),N)
      DO 640 J=N,2,-1
        DO I=1,J-1
          X(I)=X(I)-A(IP(I),J)*X(J)
        END DO
        X(J-1)=X(J-1)*A(IP(J-1),J-1)
 640    CONTINUE
      RETURN; END SUBROUTINE MATRIX_1

*******************************************************
*	inverse matrix
*******************************************************
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
			PRINT *, "this matrix is not normal" ;STOP
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
 100	CONTINUE       !Gia tri code la 100

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
	
	RETURN;	END SUBROUTINE
************************************************
Chuong trinh tinh ma tran nghich dao cua A
************************************************
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
************************************************
!Changing the value of transition matrix P to Q, Q=P.R, R is repair matrix
!Matrix R is vector at a certain value of life expectancy of rating.
************************************************
	SUBROUTINE repairmatrix ! This subroutine is for repair strategy for LCC
	! Value of d can be changed upon the life expectancy of the ratings
	USE VARIABLE;IMPLICIT NONE
	REAL*8:: R(jmax,jmax), Q(jmax,jmax), P(jmax,jmax)
	REAL*8:: E(Jmax), C(Jmax) , LCC(Jmax-1,Year)
	REAL*8:: L(Jmax-1,Year)
	REAL*8,PARAMETER:: rate = 0.9
	INTEGER I,J, t
	
	Do I=1, Jmax
		Do J=1, Jmax
			If (i==j .and. i < d) Then 
				R(i,j)=1 
			Else
				R(i,j)=0 
				If (i/=j .and. i>=d) Then
				R(i,1)=1
			End if; end if
		End Do
	End Do
	PRINT*, "Repair matrix R"
	Do I=1, Jmax
		PRINT 11, (R(i,j), j=1,Jmax)
	End Do
	
	OPEN(1005,file=FILE1005)
	Do i = 1, Jmax
		Read(1005,*) (P(i,j), j=1, Jmax)
	End Do
	PRINT*, "transition matrix P"
	Do i=1, Jmax
		PRINT 11, (P(i,j), j=1, Jmax)
	End Do
	Q=matmul(P,R)	
	PRINT*, "transition matrix Q after repair"
	Do i=1, Jmax
		PRINT 11, (Q(i,j), j=1, Jmax)
	End Do
	OPEN(2, file=file2)
		Do J=1, Jmax
			READ(2,*) C(J)
		End Do
	E=0
	Do i=1,Jmax
		Do j=i, Jmax
		E(i)=E(i)+P(i,j)*C(j) ! e_i = SUM(p_ij*c_j)
		End Do
	End Do
	PRINT*, "Value of C(i) and E(i)"
	Do i=1, Jmax
		PRINT 11, C(i), E(i)
	End Do
***************************************************************************
! Start calculating the value of LCC(Year) = rate*(E_i + SUM(q_ij*LCC(Year+1)
*************************************************************************
	
	Do i=1, Jmax-1
		LCC(i,year) = 0
	End Do
	PRINT*, " LCC at the management term"
	Do i=1, Jmax-1
		PRINT*, LCC(i,year)
	End Do
	Do t = year, 2,-1
		Do i = 1, d-1
			Do j = 1, d-1
			L(i,t)=L(i,t)+Q(i,j)*LCC(j,t)
			LCC(i,t-1)=rate*(E(i)+L(i,t))
			End Do
		End Do
	End Do
	PRINT*, "Value of LCC"
	Do i=1, d-1
		PRINT 11, (LCC(i,t), t=1,year-1)
	End Do
					
11	FORMAT(1X,6(F10.4,","),F10.4)
	CLOSE(2)
	CLOSE(1005)
	End Subroutine