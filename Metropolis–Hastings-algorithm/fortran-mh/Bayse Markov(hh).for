*******************************************************
*	Hidden Markov Using Gibbs sampling for alpha and MH algorithm for Beta
*	Version for 5 condition states
*	True condition states are improved
*	True condition states generating from Markov transition probability

*******************************************************
	MODULE FILE
		CHARACTER(LEN=30),PARAMETER:: INPUT_F1="input\kenshoudata3.csv"
		CHARACTER(LEN=30),PARAMETER:: INPUT_F2="input\info1.csv"
		CHARACTER(LEN=30),PARAMETER:: OUTPUT_F1="output\point.csv"
		CHARACTER(LEN=30),PARAMETER:: OUTPUT_P="output\p.csv"
		CHARACTER(LEN=30),PARAMETER:: OUTPUT_Q="output\q.csv"
		CHARACTER(LEN=30),PARAMETER:: OUTPUT_F="output\f.csv"
		CHARACTER(LEN=30),PARAMETER:: OUTPUT_PI="output\pi.csv"
		CHARACTER(LEN=30),PARAMETER:: OUTPUT_F2="output\result.csv"
		CHARACTER(LEN=30),PARAMETER:: OUTPUTBETA="output\beta.csv"
		CHARACTER(LEN=30),PARAMETER:: AV_THETA="output\av_theta.csv"
		CHARACTER(LEN=30),PARAMETER:: SKSAMP="output\sksamp.csv"
		CHARACTER(LEN=17),PARAMETER:: FILE1005="output\likely.csv"
		CHARACTER(LEN=17),PARAMETER:: TEST1="output\test1.csv"
		CHARACTER(LEN=17),PARAMETER:: TEST2="output\test2.csv"
	END MODULE
!------------------------------------------------------------------------
	MODULE DIM
		INTEGER,PARAMETER:: JMAX=4	!Condition states
		INTEGER,PARAMETER:: MMAX=2	!No. of covariate
		INTEGER,PARAMETER:: JMJM=JMAX*MMAX
		INTEGER,PARAMETER:: NDATA=  17384  !17384   !17384	!Number of data
!		No. of data, route 15 is removed
		INTEGER,PARAMETER:: TMAX=2	!Maximum No. of inspection
!	Conditions for caculation
		INTEGER,PARAMETER:: NPOINT=1000 !Desired sampling
		INTEGER,PARAMETER:: IBURN=800	!Burning period for Markov chain convergent
		!For for Markov chain converge
		INTEGER,PARAMETER:: TIMEMAX=1
	!	1.Metropolis, 2.Independent MH, 3.AR-MH(QPATTERN=1,3)
	!  1.Metropolis law, 2. Independent MH Law, 3.AR-MH method (QPATTERN = 1,3 only)
		INTEGER,PARAMETER:: SMETHOD=2
	!	Creating a distribution method proposed
	
!1. Approximation Taylor, 2. The process of random walk, 3. Any normal distribution, 4.ARS
		INTEGER,PARAMETER:: QPATTERN=2
	!	TRUE: severe rejection, FALSE: the rejection sampling
		LOGICAL,PARAMETER:: PREJECT=.TRUE.
	! All prior information, 2. Partial prior information, 3. No Information
		INTEGER,PARAMETER:: INFOTYPE=3
		INTEGER,PARAMETER:: LIMIT=0	!収束しなかった場合の繰り返し回数/
		! Iterations did not converge when
	END MODULE
!---------------------------------------------
	MODULE VARIABLE
	USE DIM
		REAL*8 BETA(JMAX,MMAX)
		DATA BETA(1,1)/0.245/  !0.28/
		DATA BETA(2,1)/0.025/  !0.03/
		DATA BETA(3,1)/0.055/   !0.1/
		DATA BETA(4,1)/0.047/   !0.1/
		DATA BETA(1,2)/0.5942/   !0.5/
		DATA BETA(2,2)/0.1769/  !0.25/
		DATA BETA(3,2)/0.0454/   !0.0
		DATA BETA(4,2)/0.026/   !0.0
		INTEGER,PARAMETER:: D1=2	!削除パラメータ数/
		! Number of parameters removed
		INTEGER,PARAMETER:: D2=2	!削除パラメータ数(尤度比検定用)
	! Number of parameters removed (for the likelihood ratio test)
!---------機械誤差を表すパラメータParameters representing the mechanical error-----------------------------------------------------
!---------プログラムで無茶してるので要確認/I'm sure needed programs absurd
		REAL*8,PARAMETER:: PSIGMA(JMAX+1,JMAX+1)=0.0001
!---------乱数正規用パラメータParameters for the regular random---------------------------------------------------------
	REAL*8,PARAMETER:: MYUBAR=0.0		!乱数正規用/For the regular random
	REAL*8,PARAMETER:: SIGMABAR=0.08	!乱数正規用/For the regular random
!--------------------------------------------------------------------------------------
		INTEGER DEL(D1)
		INTEGER DEL2(D2)
		REAL*8:: LIKELY(D2),SUMLIKELY(D2)	!健全なstateからのハザード関数/
		!Hazard function in a healthy state
		DATA DEL/4,8/	!削除するパラメータ番号((J-1)*MMAX+M)
			!Parameter number to delete ((J-1) * MMAX + M)
		DATA DEL2/4,8/ 	!尤度比検定用削除パラメータ
		! Parameters for likelihood ratio test delete
		REAL*8 MYU0(MMAX,JMAX+1),ROH(MMAX,MMAX,JMAX+1)
!		INTEGER IK(TMAX,NDATA),JK(TMAX,NDATA)
		REAL*8 LAMB(JMAX,TMAX,NDATA)
		REAL*8 P(JMAX+1,JMAX+1), Q(JMAX+1,JMAX+1)
		REAL*8 ZETA(MMAX,JMAX), SIGMA(MMAX,MMAX,JMAX)
		REAL*8 MYU(JMAX+1,JMAX+1), NYU(JMAX+1,JMAX+1)
		REAL*8 AVX(MMAX)	!説明変数平均値(アウトプット用)
	!! Average explanatory variables (for output)
		DATA AVX(1)/1.0/
		DATA AVX(2)/0.207814/
		INTEGER SK(TMAX,NDATA), UK(TMAX,NDATA) !ベクトルs,u/! Vector s, u
		INTEGER PATTERN
		LOGICAL KENZENDO	!どの健全度分布からか(TRUE=真：FALSE=見掛け)
!! Or sound from any distribution (TRUE = true: FALSE = apparent)
	END MODULE
!------------------------------------------------------------------------
	MODULE BOUND
	USE DIM
		LOGICAL OB(JMJM),UB(JMJM)	!領域Dの上限下限の存在
	! Presence of upper and lower limits of D region
		DATA OB/JMJM*.FALSE./
		DATA UB/JMJM*.FALSE./
		REAL*8 OX(JMJM),UX(JMJM)	!領域Dの上限下限
	! Area upper and lower limits of D
		REAL*8 HTHETA(JMJM)	!事後分布のモード候補
	!Mode candidate of posteriori distribution
		REAL*8 SVS(JMJM)
		DATA SVS(1:4)/3.5, 3.5, 3.5, 3.5/
		DATA SVS(5:8)/3.5, 3.5, 3.5, 3.5/
	!	ARSを支配するパラメータ!	Parameter that rules ARS
		INTEGER,PARAMETER:: KMAX=2	!初期基準座標数
		!Number of initial standard coordinates
		INTEGER NT	!Tk：x座標の最大数
	!The maximum number of Tk:x coordinates
		REAL*8 SXK(JMAX*MMAX,KMAX)	!スタートポイント/!Start point
	END MODULE

	MODULE DATABASE
	USE DIM
		REAL*8 XK(MMAX,TMAX,NDATA), ZK(TMAX,NDATA)
		INTEGER IK(NDATA),JK(NDATA), MK(TMAX,NDATA)
	END MODULE
!---------------------------------------------
	MODULE MTHETA
	USE DIM
		REAL*8 THETA(JMAX+1),TSA(JMAX+1,JMAX+1),PSUM(JMAX+1),ZKK
		INTEGER IKK, JKK
		INTEGER RESI, RESJ
		REAL*8 RESZ,RESX(MMAX,JMAX)
	END MODULE
***************************************************************
*	MAIN PROGRAM
***************************************************************
	PROGRAM BAYES_MARKOV
	USE FILE; USE VARIABLE; USE FILE; USE BOUND; USE DATABASE;IMPLICIT NONE
	REAL*8 SAMP(0:NPOINT,JMJM)	!サンプリング集合 Desired samples
	REAL*8 PSAMP(NPOINT,JMAX+1,JMAX+1) !Pサンプリング集合/Sample set P
	REAL*8 QSAMP(NPOINT,JMAX+1,JMAX+1) !Qサンプリング集合Sample set Q
	REAL*8 FMI(NPOINT,JMAX+1,JMAX+1)   !Fサンプリング集合Sample set F
*	乱数の種/ Random seed, number)
	INTEGER IR1,IX,IY,IZ,TANE,D
	COMMON /IRAND1/IR1
	COMMON /IRAND2/IX,IY,IZ
	CHARACTER(LEN=30):: hour
	INTEGER I,J,M,ID,K
	REAL RES

	OPEN(11,file=TEST2)
!-----変更可能箇所  Points can change-----------------------------------
*	MHの初期値/ Initial value of MH
	SAMP=0.0
	DO J=1,JMAX; DO M=1,MMAX
		ID = (J-1)*MMAX+M
		SAMP(0,ID)=BETA(J,M)
		HTHETA(ID)=BETA(J,M)
	END DO; END DO
!	SAMP=4.0*SAMP/5.0
		KENZENDO=.TRUE.
!		KENZENDO=.FALSE.
	RES=0.0
	DO I=1,JMAX
	RES=RES+I
	END DO
	DO I=1,JMAX; DO K=1,JMAX
	P(I,K)=(JMAX+1-I)/RES; Q(I,K)=(JMAX+1-I)/RES
	END DO; END DO

	MYU=1.0; NYU=1.0

!------------------------------------------------
	CALL TIME(hour); PRINT *, hour
	MYU0=0.0; ROH=0.0
*	 Random seed
	CALL CLOCK (IR1,2)
	
	IX=16251; IY=21157; IZ=3583
	TANE=IR1
*	Data input from
	CALL INPUT()

	SK=MK  ! Assumption for initial of hidden vector is actually equal to
	! observed condition states
	UK=MK ! Observed condition states.
	SELECT CASE(SMETHOD)
		CASE(1); PRINT *, "Metropolis"
		CASE(2); PRINT *, "Independent MH"
		CASE(3); PRINT *, "ARMH"
	END SELECT


	
	DO 80 D=1,D2

	CALL MCMC_MH(SAMP,PSAMP,QSAMP,FMI,D)

80	END DO	
		LIKELY=SUMLIKELY/(NPOINT-IBURN)
!	OPEN(1005,file=FILE1005)
!	DO D=1,D2
!	WRITE(1005,15) D,LIKELY(D),LIKELY(1)-LIKELY(D)
!	END DO
!	CLOSE(1005)
15	FORMAT(1X,I5,",",(F20.8,","),F12.8)

!	Conversion to initial information
	IF(INFOTYPE/=3) CALL REVERSE(2)

	!CALL OUTPUT(TANE,SAMP,PSAMP,QSAMP,FMI)
	!CALL TIME(hour); PRINT *, hour
	END PROGRAM
**************************************************************
*	/ Law MH MCMC Algorithm
**************************************************************
	SUBROUTINE MCMC_MH(SAMP,PSAMP,QSAMP,FMI,DD)
	USE VARIABLE; USE BOUND; IMPLICIT NONE
	REAL,EXTERNAL:: GEWEKEZ
	REAL*8,INTENT(OUT):: SAMP(0:NPOINT,JMJM) !samplingsetβ
	REAL*8,INTENT(OUT):: PSAMP(NPOINT,JMAX+1,JMAX+1) !Pサンプリング集合
	REAL*8,INTENT(OUT):: QSAMP(NPOINT,JMAX+1,JMAX+1) !Qサンプリング集合
	REAL*8,INTENT(OUT):: FMI(NPOINT,JMAX+1,JMAX+1)   !Fサンプリング集合
	INTEGER,INTENT(IN):: DD	! Likelihood ratio test for
	REAL*8 MYU2(MMAX,JMAX+1),SS(MMAX,MMAX,JMAX+1),SM	! Mean/SS
	REAL*8 SAMP1(NPOINT)		!SUB sampling set
	REAL Z(JMJM); COMMON /GEWEKE/Z
	INTEGER I,J,M,T,ID,IPOINT,IREV,M1,M2,K,T2,II
	LOGICAL KENTEI

!*********************
*	Burn-in
!*********************

	DO 500 T2=1,1	!平均と分散を収束させる
	! To converge to the mean and variance
	DO J=1,JMAX; DO M=1,MMAX
		ID = (J-1)*MMAX+M
		SAMP(0,ID)=BETA(J,M)
		HTHETA(ID)=BETA(J,M)
	END DO; END DO
	IPOINT=1
	CALL SAMPLING_PROCESS(1,SAMP,PSAMP,QSAMP,FMI,IPOINT,DD)
	DO J=1,JMAX
	SVS(JMAX+J)=4.0/T2
	END DO
500	END DO

	DO ID=1,JMJM
	 SAMP(0,ID) = SAMP(IBURN,ID)
	END DO
	
!*********************
*	Sampling
!*********************
	IPOINT=1; IREV=0

 100	CONTINUE
	CALL SAMPLING_PROCESS(2,SAMP,PSAMP,QSAMP,FMI,IPOINT,DD)

!*********************
*	Hypothesis testing
!*********************
	IF(DD==100)THEN	!----Skip the test force----------
	KENTEI=.TRUE.
	I=1
		DO J=1,JMAX; DO M=1,MMAX
			ID = (J-1)*MMAX+M
			DO T=1,NPOINT
				SAMP1(T)=SAMP(T,ID)
			END DO
!	Parameters for delete (review required)
!			IF(DEL(I)==ID) THEN
!				I=I+1; CYCLE
!			END IF
			Z(ID)= GEWEKEZ(SAMP1,NPOINT)
		
			IF(Z(ID) < 1.96) THEN	!5% level of significance
				PRINT 90, ID,Z(ID), "Convergence"
			ELSE
				PRINT 90, ID,Z(ID), "NOT convergence"
				KENTEI=.FALSE.
			END IF
		END DO; END DO
 90	FORMAT(1X, I3,1X,F10.2, A12)

	IF(IREV >= LIMIT) RETURN
	IF(KENTEI==.FALSE.) THEN
		IREV = IREV + 1
		IPOINT=NPOINT*3/10
		DO J=1,JMAX; DO M=1,MMAX
			ID = (J-1)*MMAX+M
			DO T=1,IPOINT
				SAMP(T,ID) = SAMP(NPOINT-IPOINT+T,ID)
			END DO
		END DO; END DO
		IPOINT=IPOINT+1
		GOTO 100
	END IF
	END IF
	RETURN; END SUBROUTINE MCMC_MH
***************************************************
*	MH algorithm
***************************************************
	SUBROUTINE SAMPLING_PROCESS(ITYPE,SAMP,PSAMP,QSAMP,FMI,IPOINT,DD)
	USE VARIABLE; USE FILE; IMPLICIT NONE
	INTEGER,INTENT(IN):: ITYPE
	REAL*8,INTENT(OUT):: SAMP(0:NPOINT,JMJM)	 !βサンプリング集合
	REAL*8,INTENT(OUT):: PSAMP(NPOINT,JMAX+1,JMAX+1) !Pサンプリング集合
	REAL*8,INTENT(OUT):: QSAMP(NPOINT,JMAX+1,JMAX+1) !Qサンプリング集合
	REAL*8,INTENT(OUT):: FMI(NPOINT,JMAX+1,JMAX+1)   !Fサンプリング集合
	INTEGER,INTENT(IN):: DD	! For the likelihood ratio test
	REAL*8, EXTERNAL:: OBJ
	! Set sampling
	INTEGER,INTENT(IN):: IPOINT
	REAL*8 X	! Returns
	LOGICAL REJ
	INTEGER IREJ(JMJM),NTIM(JMJM); COMMON /ACEP/IREJ,NTIM
	CHARACTER(LEN=30):: hour
	INTEGER T,I,ID,J,M,B,TT,L,TNUM,II,C,TTT,K,D,S
	REAL*8 RESERVE,FMID,W
	REAL*8 AVTHETA(JMAX+1)
		INTEGER NEWU(TMAX,NDATA)!	 Subroutine forα
		REAL*8 NEWNYU(JMAX,JMAX),NEWMYU(JMAX,JMAX)
		INTEGER NEWS(TMAX,NDATA)	!	サブルーチンβ用

	IREJ=0; NTIM=0
	CALL TIME(hour); PRINT *, hour
		SELECT CASE(ITYPE)
			CASE(1); PRINT *, "Burn-in period"
	TNUM=IBURN
			CASE(2); PRINT *, "Sampling period"
	TNUM=NPOINT
		END SELECT

!	FMI=0.0
	DO 821 TT=IPOINT, TNUM

			IF(MOD(TT,100)==0.0) THEN
			PRINT*, SK
			END IF
			IF(MOD(TT,2)==0.0) THEN
			PRINT *, TT,"Round"
				DO J=1,JMAX
					PRINT*, (BETA(J,M),M=1,MMAX)
				END DO
				DO J=1,JMAX
					AVTHETA(J)=(BETA(J,1)+AVX(2)*BETA(J,2))
				END DO
				PRINT*,"AVG-THETA"
				PRINT*,AVTHETA
			END IF
	! Nhu vay, gia tri mac dinh NEWU ban dau se la = 0.0
!	NEWU=1.0 ! Ta co the de gia tri mac dinh NEWU=0.0, hay bat ki, vi se khong
	! co anh huong gi toi viec subroutine GIBBS_AlPHA Tim ket qua
		
!			CALL GIBBS_ALPHA(NEWU)
!			UK=NEWU
	!print*, uk
	!stop
	! Kiem tra gia tri thu duoc UK=Alpha bang cach duoi day
!	do k=1, ndata
!		Print*, (UK(i,k), i=1, tmax)
!		write(10, *) (UK(i,k), i=1, tmax)
!	end do
!	Stop

!	Den day, ta co gia tri moi cua Alpha dua tren Dirichlet distribution
!	Gia tri Alpha =UK nay co tac dung de tinh Beta chi tai vong lap nay.
			DO M=1,JMAX+1; DO J=1,JMAX+1
				PSAMP(TT,M,J)=P(M,J)
				QSAMP(TT,M,J)=Q(M,J)
			END DO; END DO
	!FMI samples (S can be derived in the process occurs)
			DO I=1,JMAX+1; DO M=1,JMAX+1
			FMID=0.0
			DO L=1,JMAX+1
			FMID=FMID+PSAMP(TT,M,L)*QSAMP(TT,L,I)
			END DO
			FMI(TT,M,I)=FMID
			END DO; END DO

	! Bat dau tinh toan gia tri Beta bang MCMC su dung Metropolis-Hasting
	CALL NORMAL_MH(DD)
!	Kiem tra gia tri cua Beta???
	Print*, DD

!	Stop

!Remove step parameter----------------------
		IF(D1/=0) THEN
		DO S=1,JMAX; DO M=1,MMAX
		DO D=1,D1
		ID = (S-1)*MMAX+M
		IF(MMAX*(S-1)+M==DEL(D))THEN
		BETA(S,M)=0.0
		END IF
		IF(MMAX*(S-1)+M==DEL2(DD)) THEN
		BETA(S,M)=0.0
		END IF
				SAMP(TT,ID)=BETA(S,M)
		END DO;END DO
		END DO
		ELSE
		DO S=1,JMAX; DO M=1,MMAX
		ID = (S-1)*MMAX+M
		IF(MMAX*(S-1)+M==DEL2(DD)) THEN
		BETA(S,M)=0.0
		END IF
				SAMP(TT,ID)=BETA(S,M)
		END DO;END DO
		END IF
!-------------------------------------------------
			CALL GIBBS_S(NEWS)
			SK=NEWS
			W=OBJ(BETA)
			IF(TT>=IPOINT) SUMLIKELY(DD)=SUMLIKELY(DD)+W

821	END DO
	
	OPEN(123,file=SKSAMP)
	DO K=1,NDATA
			WRITE(123,22) (SK(T,K),T=1,TMAX)
	END DO
	CLOSE(123)
22	FORMAT(1X,5(I5,","))
	

	CONTAINS
!---------------------------------
****************************************************************************************
*	α  Subroutine for generating value α using Gibbs-sampling with its prior
! dirichlet distribution density function.
!---------------------------------
	SUBROUTINE GIBBS_ALPHA(NEWU)
	USE VARIABLE; USE DATABASE; IMPLICIT NONE
		INTEGER,INTENT(OUT):: NEWU(TMAX,NDATA)
		REAL*8 NEWNYU(JMAX+1,JMAX+1),NEWMYU(JMAX+1,JMAX+1)
		REAL*8 NLI(JMAX+1,JMAX+1)
		INTEGER IFAULT,I,K,L,M
		REAL*8 P1(JMAX+1), Q1(JMAX+1), NYU1(JMAX+1) !, MYU1(JMAX+1)
!	 u for update
		REAL*8 PUTKL	!PROB(U_T^K=L)*RAN
		REAL*8 MAXP		!PROB(U_T^K=L)*RANのMAX値
		!PROB (U_T ^ K = L) * RAN for MAX value
		REAL*8,EXTERNAL::RANSU1,RSEIKI
		INTEGER MTK,STK,IDUM,TKL,B
		REAL*8 RES,W
		REAL*8 PKOUHO,PMYU
		INTEGER PINT(JMAX+1)
		REAL*8 DAMT,DAM(JMAX+1)
	
		CALL MAKE_NLI(UK,SK,NLI)
	! In this urgument, UK is observed data, and SK is true condition states, which is hidden
		NEWNYU=NYU+NLI	
!	  generate q
		DO I=1,JMAX+1
	Q1=0.0
!	print*, "value of Newnyu in DM"
			DO L=1,JMAX+1
			NYU1(L)=NEWNYU(L,I)
	!		print*, NYU1(L), L
	! Xac dinh gia tri cua vector v in dirichlet distribution

			END DO
!	stop
			CALL DM(Q1,NYU1,I)

			DO L=1,JMAX+1
			Q(L,I)=Q1(L)
			END DO
		END DO
	! Den day ta xac dinh duoc gia tri moi cua Alpha =Q(L,I)
!	Generate p

	DO M=1,JMAX+1
	DAM=0; DAMT=0
	DO L=1,JMAX+1
	DO K=1,NDATA
	PMYU=M
	PINT(M)=DINT(RSEIKI(PMYU,PSIGMA(M,L))+0.5)
! DINT is intrinsic function in F66 to convert double precision into integer value
! and truncate toward 0	
	
	IF(PINT(M)>JMAX+1) THEN
	CYCLE
	PINT(M)=JMAX+1
	ELSE IF(PINT(M)<1) THEN
	CYCLE
	PINT(M)=1
	ELSE IF(PINT(M)>M+1) THEN
	CYCLE
	PINT(M)=M+1
	ELSE IF(PINT(M)<M-1) THEN
	CYCLE
	PINT(M)=M-1
	END IF

	IF(PINT(M)==L) THEN
	DAM(L)=DAM(L)+1
	END IF
	END DO
	END DO

	DO L=1,JMAX+1
	DAMT=DAMT+DAM(L)
	END DO
	DO L=1,JMAX+1
	P(M,L)=DAM(L)/DAMT
	END DO
	END DO
	
	DO M=1,JMAX+1
!	PRINT*,(P(M,L),L=1,JMAX+1)
	END DO

!	a random number generator u
	B=0
	DO T=1,TMAX; DO K=1,NDATA
		MTK=MK(T,K)
		STK=SK(T,K)
		RES=0.0
		MAXP=0.0
		PUTKL=0.0
		W=RANSU1()
			IF(T==1) THEN
			DO L=1,STK
				RES=RES+P(MTK,L)*Q(L,STK)
			END DO
			DO L=1,STK
			PUTKL=PUTKL+P(MTK,L)*Q(L,STK)/RES
			IF(W < PUTKL) THEN
			TKL=L ;EXIT
			END IF
			END DO

	ELSE IF(T/=1) THEN

			DO L=1,STK
				RES=RES+P(MTK,L)*Q(L,STK)
			END DO

			DO L=1,STK
			PUTKL=PUTKL+P(MTK,L)*Q(L,STK)/RES
			IF(W < PUTKL) THEN
			TKL=L ;EXIT
			END IF
			END DO
	END IF
*******************Eliminate the force error**************************************
	TKL=MK(T,K)
*****************************************************************************

		NEWU(T,K)=TKL
		UK(T,K)=TKL
	IF(UK(T,K)==JMAX+1) B=B+1
	END DO; END DO
	PRINT*,B,"UK"
	END SUBROUTINE GIBBS_ALPHA

!---------------------------------
*S_T^K用  Internal procedures for S_T ^ K
!---------------------------------
	SUBROUTINE GIBBS_S(NEWS)
	USE VARIABLE; 	USE DATABASE ; IMPLICIT NONE
	INTEGER,INTENT(OUT):: NEWS(TMAX,NDATA)
	REAL*8,EXTERNAL:: PROB
	REAL*8,EXTERNAL:: RANSU1
	REAL*8 FJM(JMAX+1,JMAX+1),THETA(JMAX+1)
	REAL*8 WJT(JMAX+1,TMAX)
	REAL*8 RES,MAXP,PSTKI,W,A,B,C,D,E
	INTEGER I,L,M,TK,J,IDUM,TKL,K
	DO 300 K=1,NDATA
	RES=0.0 ; THETA=0.0
	DO I=1,JMAX; DO M=1,JMAX

	RES=0.0
		DO L=1,I
			RES=RES+P(M,L)*Q(L,I)
		END DO
		FJM(I,M)=RES
	END DO; END DO

	WJT=0.0
	DO TK=1,TMAX
	THETA=0.0 
	DO J=1,JMAX
	IF(J==JMAX+1) THEN
	THETA(J)=0.0
	ELSE
			DO M = 1, MMAX
				THETA(J)=THETA(J) + BETA(J,M)*XK(M,TK,K)
			END DO
	END IF
	END DO
	DO J=1,JMAX
		theta(jmax+1)=0.0
		IF(TK==1) THEN
!!!Check out!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			WJT(J,TK)=PROB(1,J,ZK(1,K),THETA)*PROB(J,SK(2,K),ZK(1,K),THETA)

		ELSE IF(2 <= TK .AND. TK <= TMAX-1) THEN
			WJT(J,TK)=PROB(SK(TK-1,K),J,ZK(TK,K),THETA)*
     & PROB(J,SK(TK+1,K),ZK(TK,K),THETA)
		ELSE
			WJT(J,TK)=PROB(SK(TMAX-1,K),J,ZK(TK,K),THETA)
		END IF
	END DO

	END DO
	DO TK=1,TMAX
	RES=0.0
! Stock room to reconsider (?)
		IF(TK==1)THEN	! TK = 1 is special
		DO J=1,SK(TK+1,K)
			RES=RES+WJT(J,TK)*FJM(J,MK(TK,K))
		END DO
		ELSE
		DO J=SK(TK-1,K),SK(TK,K)
			RES=RES+WJT(J,TK)*FJM(J,MK(TK,K))
		END DO
		END IF

	PSTKI=0.0
	W=RANSU1()

	IF(TK==1)THEN
	DO I=1,SK(TK+1,K)
	PSTKI=PSTKI+WJT(I,TK)*FJM(I,MK(TK,K))/RES
	IF(W < PSTKI) THEN
	TKL=I ;EXIT
	END IF
	END DO

	ELSE
	DO I=SK(TK-1,K),JMAX
	PSTKI=PSTKI+WJT(I,TK)*FJM(I,MK(TK,K))/RES
	IF(W < PSTKI) THEN
	TKL=I ;EXIT
	END IF
	END DO
	END IF

*******************Eliminate the force error************************************
	TKL=MK(TK,K)
*****************************************************************************

		NEWS(TK,K)=TKL
		SK(TK,K)=TKL

	END DO
300	END DO
	END SUBROUTINE GIBBS_S
!---------------------------------------------------------

	END SUBROUTINE SAMPLING_PROCESS
*******************************************************
*	MCMC su dung Metropolis-Hasting Algorithm
*******************************************************
	SUBROUTINE NORMAL_MH(DD)
	USE VARIABLE; USE FILE; USE BOUND; USE DATABASE;IMPLICIT NONE
	REAL*8, EXTERNAL:: OBJ,RANSU1,RSEIKI
	INTEGER,INTENT(IN):: DD	! For the likelihood ratio test
	REAL*8:: BETAA(JMAX,MMAX)
	REAL*8:: DBETA(JMAX,MMAX)
	REAL*8:: BETAD(JMAX,MMAX)
	REAL*8:: BETAE(JMAX,MMAX)
	REAL*8:: POS(JMAX,MMAX),W,H,LR(JMAX,MMAX)
	REAL*8:: OBJE(D2)
	REAL*8:: RNOR	!regular random
	REAL*8:: R
	INTEGER:: S,TIME,M,S1,M1,D,D22,D222,T
	REAL*8,PARAMETER:: EPS=1.0D-15
	OBJE=0.0
	
	! Initial value of Beta has been defined in MODULE VARIABLE

	DO 88 T=1,1
!	---Remove step parameter-------------------
	IF(D1/=0) THEN
	DO S=1,JMAX; DO M=1,MMAX
	DO D=1,D1

	IF(MMAX*(S-1)+M==DEL(D))THEN
	BETA(S,M)=0.0
	END IF
	IF(MMAX*(S-1)+M==DEL2(DD)) THEN
	BETA(S,M)=0.0
	END IF
	END DO;END DO
	 END DO
	ELSE
	DO S=1,JMAX; DO M=1,MMAX
	IF(MMAX*(S-1)+M==DEL2(DD)) THEN
	BETA(S,M)=0.0
	END IF
	END DO;END DO
	END IF
	
	! Initial value of Beta has been defined in MODULE VARIABLE

!-------------------------------------------------
	BETAD=BETA

	DO 22 TIME=1,TIMEMAX

	DO M=1,MMAX; DO S=1,JMAX

	IF(MMAX*(S-1)+M==DEL2(DD)) GO TO 55

	DO D=1,D1
	IF(MMAX*(S-1)+M==DEL(D)) GO TO 55
	END DO
	R=RANSU1()
	RNOR=RSEIKI(MYUBAR,SIGMABAR)

			BETAA=BETAD
			POS=BETAD
			BETAA(S,M)=BETAA(S,M)+RNOR

	H=OBJ(BETAD)
			IF(OBJ(BETAA)-H > 0.0) THEN
				POS(S,M)=BETAA(S,M)
			ELSE IF(OBJ(BETAA)-H > DLOG(R)) THEN
				POS(S,M)=BETAA(S,M)
			ELSE
			END IF
			BETAD=POS

55	END DO
	END DO

			BETA=BETAD

20	FORMAT(1X,3(F12.8,","),F12.8)
	
	W=OBJ(BETAD)
	PRINT*,"Log likelihood",W

22	END DO
88	END DO
      DBETA=BETA
!	Kiem tra gia tri moi cua Beta.
!	do k=1, jmax
		Print*, "Gia tri moi cua Beta"
		Print*, DBETA
!	end do
!	stop
99	CONTINUE

	END SUBROUTINE NORMAL_MH
****************************************************
*	Log-likelihood function / For complete posterior distribution
****************************************************
	REAL*8 FUNCTION OBJ(DBETA)
	USE VARIABLE; USE FILE; USE BOUND; USE DATABASE
	USE MTHETA;IMPLICIT NONE
	REAL*8,INTENT(IN):: DBETA(JMAX,MMAX)
	REAL*8, EXTERNAL:: PROB
	REAL*8 THETASA(JMAX+1,JMAX+1)
	REAL*8 XKK(MMAX,JMAX)
	INTEGER I,J,L,R,K,M,MINTHETA
	INTEGER EK,FK,GK ! Former condition states
	INTEGER V,X,Y,TK,NK
	REAL*8:: ZK1,ZK2,Z,W,U
	REAL*8 PI,QQ,AA
	OBJ=0.0
	V=0
	L=0
	DO 101 NK=1,NDATA ! Loading data
	DO TK=1,TMAX
	!! Between the two as a weak state, inspection intervals, characteristic vector
!	 Sampling or from anywhere
	IF(KENZENDO==.TRUE.)THEN
		!  T = 1 is special
		IF(TK==1) THEN
		IK(NK)=1; JK(NK)=SK(TK,NK)
		IKK=1; JKK=SK(TK,NK)
		ELSE
		IK(NK)=SK(TK-1,NK); JK(NK)=SK(TK,NK)	!SK
	!! SK soundness of each talk extracts from
		IKK=SK(TK-1,NK); JKK=SK(TK,NK)	!SK
	!! SK soundness of each talk extracts from
		END IF
	ELSE
		IF(TK==1) THEN
		IK(NK)=1; JK(NK)=MK(TK,NK)
		IKK=1; JKK=MK(TK,NK)
		ELSE
		IK(NK)=MK(TK-1,NK); JK(NK)=MK(TK,NK)	!MK
	! MK extracted from the soundness of each talk
		IKK=MK(TK-1,NK); JKK=MK(TK,NK)	!MK
		END IF
	END IF
	!! Between the two as a weak state, the inspection period, characteristic vector
		DO J=1,JMAX;DO M=1,MMAX
			XKK(M,J)=XK(M,TK,NK)
		END DO; END DO

		DO I=1,D1	!  Removed from the parameter space
			M = MOD(DEL(I)+(MMAX-1),MMAX) + 1; J = (DEL(I)-M)/MMAX+1
			XKK(M,J)=0.0
		END DO

		! θ values
		CALL THETAVALUE(IKK,JKK,XKK,DBETA,THETASA,MINTHETA)
		PI=PROB(IKK,JKK,ZK(TK,NK),THETA)
		!! Likelihood function
		IF(PI<=1.0D-305) THEN
			OBJ=OBJ-702.0
!	print*,"q",k
		ELSE IF(PI > 1.0000) THEN
	OBJ=OBJ
	ELSE
		QQ=DLOG(PI)
			OBJ=OBJ+QQ
		END IF
	END DO
101	END DO

	CLOSE(1)
	CLOSE(1009)
20	FORMAT(1X,5(I5,","),F12.8)
	RETURN
	END FUNCTION
*******************************************
*	Θi (i = IK, JK) find a value
* This subroutine can be refered to conventional Markov program
! This subroutine define value of Θi, Θi-Θj, and Reserve
*******************************************
	SUBROUTINE THETAVALUE(TIK,TJK,XKK,BETA,THETASA,MINTHETA)
	USE DIM; USE MTHETA; IMPLICIT NONE
	REAL*8,INTENT(IN):: XKK(MMAX,JMAX),BETA(JMAX,MMAX)
	REAL*8,INTENT(OUT):: THETASA(JMAX+1,JMAX+1)
	INTEGER,INTENT(IN):: TIK,TJK
	INTEGER,INTENT(OUT):: MINTHETA
	INTEGER I,J,M
	REAL*8 RESERVE

	!Caculating Θ
	THETA=0.0
	DO J=TIK,TJK
		IF(J==JMAX+1) THEN
	THETA(J)=0.0
	ELSE
		DO M=1,MMAX
			THETA(J)=THETA(J) + BETA(J,M)*XKK(M,J)
		END DO

	END IF
	END DO
	!	Θi-Θj
	THETASA=0.0
	DO J=TIK,TJK; DO I=TIK,TJK
		THETASA(I,J)=THETA(I)-THETA(J)
	END DO;	END DO
		
	!	Minimum to prevent overflow Θ
		MINTHETA=TIK; RESERVE=DEXP(600.0D0)
	DO J=TIK,TJK
		IF(J==JMAX+1) CYCLE
		IF(RESERVE > THETA(J)) THEN
			RESERVE=THETA(J); MINTHETA=J
		END IF
	END DO
	
	END SUBROUTINE


***********************************************************************
*	Input from the input file
***********************************************************************
	SUBROUTINE INPUT()
	USE DATABASE; USE VARIABLE; USE FILE; IMPLICIT NONE
	REAL*8 MAXX(MMAX),XK2(MMAX,NDATA),XKD(MMAX,NDATA)
	COMMON /ARR/MAXX
	CHARACTER CH
	INTEGER I,J,K,M,M1,M2,T
!	SCALE=1.0; 
	MAXX=0.0; MAXX(1)=1.0
	XK=1.0
	OPEN(1000,file=INPUT_F1)
		READ(1000,*) CH
		DO K=1,NDATA
	READ(1000,*) (MK(T,K),T=1,TMAX),(ZK(T,K),T=1,TMAX),
     & (XK2(M1,K),M1=2,MMAX)
			DO M1=2,MMAX; IF(XKD(M1,K) > MAXX(M1)) MAXX(M1)=XKD(M1,K)
	! The above line can be unnecessary.
			 END DO
!		スケール調整  Adjustment Scale
	DO T=1,TMAX; DO M=2,MMAX
	XK(M,T,K)=XK2(M,K)	
	END DO; END DO
	END DO
	CLOSE(1000)
	!----------------------------------------------------------------------------------
	! It is noted here that value os XK2 as covariate values can be more than 2
	! which corresponds to the number of intervals. Therefore we use following lines instead of above line
!	READ(1000,*) (MK(T,K),T=1,TMAX),(ZK(T,K),T=1,TMAX),
 !    & (XK2(M1,K),M1=1,MMAX)   ! here we change value M1=1, Mmax
!			DO M1=1,MMAX; IF(XKD(M1,K) > MAXX(M1)) MAXX(M1)=XKD(M1,K)
	! The above line can be unnecessary.
!			 END DO
!		スケール調整  Adjustment Scale
!	DO T=1,TMAX; DO M=1,MMAX
!	XK(M,T,K)=XK2(M,K)	
!	END DO; END DO
!	END DO
! Depending on the availability of data, we can adjust M1=2, Mmax...Then we need to add up one line
! with values of covariates into input file.
!---------------------------------------------------------------------------------------

	! After this subroutine, we can call value of MK as condition states values
	! ZK as time intervals values
	! XK2 as value of covarates.

!	Print*, MK(1,14889), MK(2,14889)
!	Print*, ZK(1,14889), ZK(2,14889)
!	Print*,  XK2(2,14889)

!	stop

	IF(INFOTYPE==3) RETURN	!! If no information

	OPEN(1001,file=INPUT_F2)
		DO J=1,JMAX
			READ(1001,*) CH
			READ(1001,*) CH,(MYU0(M1,J),M1=1,MMAX)
			DO M1=1,MMAX
				READ(1001,*) CH,(ROH(M1,M2,J),M2=1,MMAX)
			END DO
		END DO
	CLOSE(1001)

	PRINT *, "Prior information [Average]"
	DO J=1,JMAX
		PRINT 30, J,(MYU0(M1,J),M1=1,MMAX)
	END DO
30	FORMAT(1X,I2,6F10.4)

!	
!Scale adjustment to the inverse transformation
	CALL REVERSE(1)

	RETURN;	END SUBROUTINE INPUT
***************************************************
*	Subroutine OUTPUT
***************************************************
	SUBROUTINE OUTPUT(TANE,SAMP,PSAMP,QSAMP,FMI)
	USE VARIABLE; USE BOUND; USE FILE; USE MTHETA; IMPLICIT NONE
	REAL,EXTERNAL:: PROB, GEWEKEZ
	INTEGER,INTENT(IN):: TANE
	REAL*8,INTENT(IN):: SAMP(0:NPOINT,JMJM)	!サンプリング集合
	REAL*8,INTENT(IN):: PSAMP(NPOINT,JMAX+1,JMAX+1) !Pサンプリング集合
	REAL*8,INTENT(IN):: QSAMP(NPOINT,JMAX+1,JMAX+1) !Qサンプリング集合
	REAL*8,INTENT(IN):: FMI(NPOINT,JMAX+1,JMAX+1) !FMIサンプリング集合
	REAL*8 SAMP1(NPOINT)	!集約化
	REAL*8 MYU2(MMAX,JMAX),SS(MMAX,MMAX,JMAX),SM	!平均，分散
	REAL*8 UBAR(JMAX*MMAX),OBAR(JMAX*MMAX)	!信頼区間
	REAL Z(JMAX*MMAX); COMMON /GEWEKE/Z	!GEWEKE
	LOGICAL REJ
	INTEGER IREJ(JMJM),NTIM(JMJM); COMMON /ACEP/IREJ,NTIM
	INTEGER I,J,K,M1,M2,M
	REAL*8 PI(JMAX+1,JMAX+1),ZZ,DPI
	
!--------------------------
*	サンプリングの出力   Output sampling
!--------------------------
	OPEN(2000,file=OUTPUT_F1)
		WRITE(2000,*) "Random seed",",",TANE
		WRITE(2000,*) "Samples",",",NPOINT
		WRITE(2000,*) "All trial frequency"
		WRITE(2000,21)  (NTIM(I),I=1,JMJM)
		WRITE(2000,*) "Refusal frequency"
		WRITE(2000,21) (IREJ(I),I=1,JMJM)
		WRITE(2000,*) "state1",",,,","state2",",,,","state3",",,,","state4"
		WRITE(2000,*) ("B1",",","B2",",","B3",",",J=1,4)
		DO K=1,NPOINT
			WRITE(2000,20) (SAMP(K,J),J=1,JMAX*MMAX)
		END DO
	CLOSE(2000)
	!FMI
 	OPEN(7000,file=OUTPUT_F)
		WRITE(7000,*) "FMI"
		DO K=1,NPOINT
			WRITE(7000,20) (FMI(K,J,1),J=1,JMAX+1),(FMI(K,J,2),J=1,JMAX+1)
     & ,(FMI(K,J,3),J=1,JMAX+1),(FMI(K,J,4),J=1,JMAX+1),
     &(FMI(K,J,5),J=1,JMAX+1)
		END DO
	CLOSE(7000)


	!PML
	OPEN(5000,file=OUTPUT_P)
		WRITE(5000,*) "P","M","L"
		DO K=1,NPOINT
			WRITE(5000,20) (PSAMP(K,J,1),J=1,JMAX+1),(PSAMP(K,J,2),J=1,JMAX+1)
     & ,(PSAMP(K,J,3),J=1,JMAX+1),(PSAMP(K,J,4),J=1,JMAX+1),
     &(PSAMP(K,J,5),J=1,JMAX+1)
		END DO
	CLOSE(5000)

	!QLI
	OPEN(6000,file=OUTPUT_Q)
		WRITE(6000,*) "Q","L","I"
		DO K=1,NPOINT
			WRITE(6000,20) (QSAMP(K,J,1),J=1,JMAX+1),(QSAMP(K,J,2),J=1,JMAX+1)
     & ,(QSAMP(K,J,3),J=1,JMAX+1),(QSAMP(K,J,4),J=1,JMAX+1),
     &(QSAMP(K,J,5),J=1,JMAX+1)
		END DO
	CLOSE(6000)


20	FORMAT(1X,30(E15.8,","))
 21	FORMAT(1X,30(I7,","))
!--------------------------
*	 Output of statistic
!--------------------------
!	Sample average
	MYU2(1,:)= 0.0 ! 1.0
	MYU2(2,:)=0.0 !0.199
	DO J=1,JMAX; DO M1=1,MMAX
		DO K=1,NPOINT
			MYU2(M1,J)=MYU2(M1,J)+SAMP(K,(J-1)*MMAX+M1)
		END DO
	END DO;	END DO	
	MYU2=MYU2/NPOINT

!	Mean value calculation of θ
	THETA=0.0
		DO J=1,JMAX
		DO M=1,MMAX
			THETA(J)=THETA(J) + MYU2(M,J)*AVX(M)
		END DO
!		THETA(J) = DEXP(THETA(J))
		END DO

	OPEN(100,file=AV_THETA)
	DO J=1,JMAX
	WRITE(100,20) THETA(J)
	END DO
	CLOSE(100)

	!PI
	ZZ=1.0
		DO I=1,JMAX+1; DO J=1,JMAX+1
	DPI=PROB(I,J,ZZ,THETA)
	PI(I,J)=DPI
		END DO; END DO
 	OPEN(8000,file=OUTPUT_PI)
		WRITE(8000,*) "PI"
	DO I=1,JMAX+1
			WRITE(8000,20) (PI(I,J),J=1,JMAX+1)
	END DO
	CLOSE(8000)


!  Sample variance
	DO J=1, JMAX; DO M2=1,MMAX; DO M1=1,MMAX
		I = (J-1)*MMAX
		DO K=1,NPOINT
			SS(M1,M2,J)=SS(M1,M2,J)
     &			+(SAMP(K,I+M1)-MYU2(M1,J))*(SAMP(K,I+M2)-MYU2(M2,J))
		END DO
	END DO; END DO; END DO
	SS=SS/NPOINT

!	 Confidence interval
	DO J=1,JMAX*MMAX
		DO K=1,NPOINT; SAMP1(K)=SAMP(K,J);	END DO
		CALL BUBBLE(SAMP1,NPOINT) ! Rearrange the obtained data in order by 
	                               ! using BUBBLE sorting algorithm
		UBAR(J)=SAMP1(NPOINT*5/100)
		OBAR(J)=SAMP1(NPOINT*95/100)
		Z(J)= GEWEKEZ(SAMP1,NPOINT)
		END DO
 90	FORMAT(1X, F10.2, A12)

!	 File output
	OPEN(3000,file=OUTPUT_F2)
		SELECT CASE(INFOTYPE)
			CASE(1); WRITE(3000,*) "Prior all information"
			CASE(2); WRITE(3000,*) "Partial prior information"
			CASE(3); WRITE(3000,*) "No Information"
		END SELECT
		SELECT CASE(SMETHOD)
			CASE(1); WRITE(3000,*) "Sampling methods",",","Metropolis"
			CASE(2); WRITE(3000,*) "Sampling methods",",","Independent MH"
			CASE(3); WRITE(3000,*) "Sampling methods",",","ARMH"
		END SELECT
	

		SELECT CASE(QPATTERN)
	CASE(1); WRITE(3000,*) "Prodistribution",",","Taylor"
	CASE(2); WRITE(3000,*) "Prodistribution",",","Randomwalk"
			CASE(3); WRITE(3000,*) "Proposeddistribution",",","Regularproper"
			CASE DEFAULT; WRITE(3000,*) "Unmaking"
		END SELECT

		IF(INFOTYPE/=3) THEN
			WRITE(3000,*) "Initial information"
			WRITE(3000,*) 
     &		",", "state1",",,,","state2",",,,","state3",",,,","state4",",,,"
			WRITE(3000,10) "Initial average",((MYU0(M1,J),M1=1,MMAX),J=1,JMAX)
			DO M1=1,MMAX
				WRITE(3000,10) "Inidispersion",((ROH(M1,M2,J),M2=1,MMAX),J=1,JMAX)
			END DO
		END IF
				
		WRITE(3000,*) "Posterior information"
		WRITE(3000,*) 
     &		",", "state1",",,,","state2",",,,","state3",",,,","state4",",,,"
		WRITE(3000,10) "Sample average",((MYU2(M1,J),M1=1,MMAX),J=1,JMAX)
		DO M1=1,MMAX
			WRITE(3000,10) "Sample variance",((SS(M1,M2,J),M2=1,MMAX),J=1,JMAX)
		END DO

		WRITE(3000,10) "Confidence interval",(UBAR(J),J=1,JMAX*MMAX)
	! Upper level for testing (equivalent to 5%)
		WRITE(3000,10) "Confidence interval",(OBAR(J),J=1,JMAX*MMAX)
	! Lower level for testing (Equivalent to 95% )
		WRITE(3000,10) "Convergence criterion",(Z(J),J=1,JMAX*MMAX) 
	! Z(J) is result of Geweke test
 10	FORMAT(1X,A8,",",20(E15.8,","))
	CLOSE(3000)
	

	OPEN(1111,file=OUTPUTBETA)
		DO J=1,JMAX
			WRITE(1111,15) (MYU2(M1,J),M1=1,MMAX)
		END DO
	CLOSE(1111)

 15	FORMAT(2(F15.8,","))


	RETURN;	END SUBROUTINE OUTPUT
**********************************************
*	Scale-adjustment
**********************************************
	SUBROUTINE REVERSE(IP)
	USE VARIABLE; IMPLICIT NONE
	INTEGER,INTENT(IN):: IP
	REAL*8 A(MMAX,MMAX)
	REAL*8 MAXX(MMAX)	!Maxium numbers of covariates 
	!The maximum value of explanatory variable
	COMMON /ARR/MAXX
	INTEGER J,M1,M2
	
	IF(IP==1) THEN
		DO J=1,JMAX
			DO M1=2,MMAX; MYU0(M1,J)=MYU0(M1,J)*MAXX(M1); END DO
			DO M1=1,MMAX; DO M2=1,MMAX
				A(M2,M1) = ROH(M2,M1,J)*MAXX(M2)*MAXX(M1)
			END DO; END DO
			CALL MATRIX_INVERSE(A,MMAX)
			DO M1=1,MMAX; DO M2=1,MMAX
				ROH(M2,M1,J) = A(M2,M1)
			END DO; END DO
		END DO
	ELSE	! ! Reverse-adjustment
		DO J=1,JMAX
			DO M1=2,MMAX; MYU0(M1,J)=MYU0(M1,J)/MAXX(M1); END DO
			DO M1=1,MMAX; DO M2=1,MMAX
				A(M2,M1) = ROH(M2,M1,J)
			END DO; END DO
			CALL MATRIX_INVERSE(A,MMAX)
			DO M1=1,MMAX; DO M2=1,MMAX
				ROH(M2,M1,J) = A(M2,M1)/MAXX(M2)/MAXX(M1)
			END DO; END DO
		END DO
	END IF
	
	RETURN; END SUBROUTINE
*******************************************************
*  Convergence criterion to stationary distribution (testings of hypotheses of GEWEKE)
! Please refer to Geweke testing formula which is written in the paper.
! Geweke, 1996. Evaluating the accuracy of sampling-based approaches to the caculation of posterior moments.
! See the convergent condition in Geweke test at http://www.cryoung.org/www/software/MCMCthin/Rcode.txt
*******************************************************
	REAL FUNCTION GEWEKEZ(SAMP,NPOINT)
	IMPLICIT NONE
	REAL*8,EXTERNAL:: TVORA
	INTEGER,INTENT(IN):: NPOINT ! ! Number of samples
	REAL*8,INTENT(IN):: SAMP(NPOINT)	!!Specimen set
	INTEGER NP1, NP2
	REAL*8 SAMP1(NPOINT/10), SAMP2(NPOINT/50)
	REAL*8 MYU1,MYU2,S1,S2
	REAL*8 Z	! ! Test statistic
	INTEGER I
	CHARACTER CH

		NP1=NPOINT/10; NP2=NPOINT/50
		DO I=1, NP1; SAMP1(I)=SAMP(I); END DO
		DO I=1, NP2; SAMP2(I)=SAMP(I+NP2); END DO

		MYU1=0.0
		DO I=1,NP1; MYU1=MYU1+SAMP1(I);	END DO
		MYU1 = MYU1 / NP1
		S1=TVORA(SAMP1,MYU1,NP1) 
		
		MYU2=0.0
		DO I=1,NP2; MYU2=MYU2+SAMP2(I);	END DO
		MYU2 = MYU2 / NP2
		S2=TVORA(SAMP2,MYU2,NP2)
		
		GEWEKEZ=DABS((MYU1-MYU2) / SQRT(S1/NP1 + S2/NP2))

	RETURN
	END FUNCTION GEWEKEZ
*****************************************
* Caculating sample covariance matrix.
! Calculation of sample variance (time series) Reference literature
! A Simple, Positive Semi-Definite, Heteroskedasticity
! and Autocorrelation Consistent Covariance Matrix
! Read the material in reference forlder (whitney).
*****************************************
	REAL*8 FUNCTION TVORA(SAMPX,MYU,NN)
	IMPLICIT NONE
	INTEGER,INTENT(IN):: NN
	REAL*8,INTENT(IN):: SAMPX(NN),MYU
	INTEGER I,J,M
!	Please refer to simplified function for caculating S_T in theorem 1 of the paper cite{whitney}.
		M = 10	! Bandwidth
		TVORA= OMEGA(0)
		DO J=1, M
			TVORA = TVORA + (1.0 - DBLE(J)/(M+1.0))*2.0*OMEGA(J)
	! TVORA is S_T in the paper Whitney, it is asymptotic covariance matrix of parameter beta.
		END DO
	CONTAINS
!--------------------
* Internal procedure
!--------------------
		REAL*8 FUNCTION OMEGA(S)
		IMPLICIT NONE
		INTEGER,INTENT(IN):: S
		INTEGER K
			OMEGA=0.0
			DO K = S + 1, NN
				OMEGA=OMEGA +(SAMPX(K)-MYU)*(SAMPX(K-S)-MYU)
!	See Theorem 1 to understand formula (SAMPX(K)-MYU)*(SAMPX(K-S)-MYU)
!	It is actually deriving from formula (c'h_t)(c'h_{t-j})
			END DO
			OMEGA = OMEGA / NN
	! Then NN is total number of testing point.
		END FUNCTION
!---------------------------------------------------------
	END FUNCTION TVORA
*******************************************************

***************************************************
! This function RANSU1() is also used to generate random number in range (0,1)
***************************************************
	REAL*8 FUNCTION RANSU1()
	INTEGER IR1
	COMMON /IRAND1/IR1
		IR1=IAND(48828125*IR1, 2147483647)
		RANSU1=DBLE(IR1)/DBLE(2147483647)
	END FUNCTION
***************************************************
* Refer to code of Fortran：AS183-Generation of uniform pseudo-random number
!Obtained results will be random numbers in the range (0,1)
***************************************************
	REAL*8 FUNCTION RANSU2()
	INTEGER IX,IY,IZ
	COMMON /IRAND2/IX,IY,IZ

		IX=171*MOD(IX,177)-2*(IX/177)
		IY=172*MOD(IY,176)-35*(IY/176)
		IZ=170*MOD(IZ,178)-63*(IZ/178)

		IF(IX < 0) IX=IX+30269
		IF(IY < 0) IY=IY+30307
		IF(IZ < 0) IZ=IZ+30323
	   
		RANSU2=AMOD(FLOAT(IX)/30269.0+FLOAT(IY)/30307.0
     &		+FLOAT(IZ)/30323.0, 1.0)
	RETURN;	END FUNCTION
***************************************************
*	Bubble Sorting Algorithm/ Sap xep Noi Bot
* Refer to ready made functions and subroutine from
* http://en.wikibooks.org/wiki/Algorithm_Implementation/Sorting/Bubble_sort
* Following subroutine is defined in similarity with BASIC code
* For details of Bubble Sort definition, please refer to http://en.wikipedia.org/wiki/Bubble_sort

***************************************************
	SUBROUTINE BUBBLE(RR,NR)
	INTEGER,INTENT(IN):: NR
	REAL*8 RR(NR)
	INTEGER I,J
	REAL*8 RES

	DO I=NR,2,-1; DO J=1,I-1
		IF(RR(J) > RR(J+1)) THEN
			RES=RR(J);	RR(J)=RR(J+1);	RR(J+1)=RES
		END IF
	END DO;	END DO
	
	END SUBROUTINE
*******************************************************
! Caculating Invert Matrix, we can use the gamma subroutine in Fortran recipe
! to substitute for this Matrix_inverse (A,N) subroutine
! Using this subroutine, we can estimate inverst of A with N is matrix dimension.
*******************************************************
	SUBROUTINE MATRIX_INVERSE(A,N)
	IMPLICIT NONE
	INTEGER,INTENT(IN):: N
	REAL*8:: A(N,N),AW,AW2,MM(N)
	INTEGER:: K,I,J,P,Q
	REAL*8 EPS
	INTEGER LW,MW,L(N),M(N)

	EPS=1.0D-75

	DO J=1,N
		L(J)=J; M(J)=J
	END DO

	DO 100 K=1,N
		AW=DABS(A(K,K))
		P=K; Q=K
		DO J=K,N; DO I=K,N
			IF(AW < DABS(A(I,J))) THEN
				AW=DABS(A(I,J)); P=I; Q=J
			END IF
		END DO; END DO
		IF(AW < EPS) THEN
			PRINT *, "Singular matrix" ;STOP
		END IF
		IF(K /= P) THEN
			DO J=1,N
				AW=A(K,J); A(K,J)=A(P,J); A(P,J)=AW
			END DO
			MW=M(K);M(K)=M(P);M(P)=MW
		END IF
		IF(K /= Q) THEN
			DO I=1,N
				AW=A(I,K); A(I,K)=A(I,Q); A(I,Q)=AW
			END DO
			LW=L(K);L(K)=L(Q);L(Q)=LW
		END IF
		AW=1.0/A(K,K);	A(K,K)=1.0
		DO J=1,N; A(K,J)=A(K,J)*AW;	END DO
		DO I=1,N
			IF(I /= K) THEN
				AW=A(I,K); A(I,K)=0.0
				DO J=1,N; A(I,J) = A(I,J)-AW*A(K,J); END DO
			END IF
		END DO
 100	CONTINUE

	DO J=1,N
		DO I=1,N
			P=L(I);	MM(P)=A(I,J)
		END DO
		DO I=1,N; A(I,J)=MM(I); END DO
	END DO

	DO I=1,N
		DO J=1,N
			Q=M(J);	MM(Q)=A(I,J)
		END DO
		DO J=1,N; A(I,J)=MM(J);	END DO
	END DO
	RETURN;	END SUBROUTINE
***********************************************************************
*	Subroutine for defining  N^{i,(n-1)}_l calculation
*********************************************
	SUBROUTINE MAKE_NLI(UTK,STK,NLI)
	USE VARIABLE; IMPLICIT NONE
	REAL*8,INTENT(OUT):: NLI(JMAX+1,JMAX+1)
	INTEGER,INTENT(IN):: UTK(TMAX,NDATA) !Observed condition states M^k_t
	 ! Degree of real choice of health
	INTEGER,INTENT(IN):: STK(TMAX,NDATA) !True condition states
	INTEGER L,I
	REAL*8 RES
	INTEGER J,K

	NLI=0
	DO L=1,JMAX+1; DO I=1,JMAX+1
	RES=0	
	DO J=1,NDATA; DO K=1,TMAX	!
	IF(UTK(K,J)==L .AND. STK(K,J)==I) RES=RES+1
	END DO; END DO

	NLI(L,I)=RES	
	END DO; END DO
	END SUBROUTINE

*********************************************
*	M^{L(n-1)}_M
*********************************************

***************************************************
* Dirichlet distribution (multi-covariates)
! K is number of covariates taking integer value k must >1
! a is parameter  ν
! x is return values after running the subroutine. 
! Values of x(i) is actual probability
***************************************************
      subroutine Dm(x,a,k)
      INTEGER,INTENT(IN):: k		!number of covariates k
	REAL*8,INTENT(IN):: a(k)	!Initial value of ν
	REAL*8,INTENT(OUT):: x(k)		!Return value of x(i)
	INTEGER I
	REAL*8 y(k),b,c
	real*8 gene,sum
      sum=0.0
      do i=1,k
	 b=a(i)
	 c=1.0
       call gamma(gene,B,c)
       y(i)=gene
	 sum=sum+y(i)
	end do
!	Print*, "Value of x(i)"
      do i=1,k
	 x(i)=y(i)/sum ! This is for satisfying the sum of x(i)<=1
!	print*, x(i) 
!	write(11,*) x(i) We can check the result of summing up =1
	end do

	return
	end
****************************************************
*   gamma(x,k,te) is generation function with shape parameter k and location parameter te.
!   x is result after running this subroutine, parameter in the input is k and te.
! Refer ready made subroutine is algorithm 654 in the paper 
!FORTRAN subroutines for computing the incomplete gamma function ratios and their inverse
! Output of this subroutine is values of x, which is in relation with shape parameter k
! and location parameter te. Ranking value of x is positive values >0
****************************************************
      subroutine gamma(x,k,te)
	REAL*8,INTENT(IN):: K,TE
	REAL*8,INTENT(OUT):: X
	REAL*8,EXTERNAL::RAN1
      integer idum
	real*8 v0,del,eps,nm
      data e /2.71828/

      del=k-int(k)
      if(del==0.0) then
	 eps=0.0

	else

      do
       u=ran1(idum)	!
       v=ran1(idum)	!

      v0=exp(1.0)/(exp(1.0)+del)
      if(u<=v0) then
	 eps=(u/v0)**(1.0/del)
	 nm=v*eps**(del-1.0)
	else
	 eps=1-log((u-v0)/(1-v0))
	 nm=v*exp(-eps)
	end if

	if(nm<=eps**(del-1.0) * exp(-eps)) exit
      end do
      end if
	sum=0.0
	do i=1,int(k)
	 sum=sum+log(ran1(idum))
	end do
	x=te*(eps-sum)
	return
	end subroutine
***************************************************
! Random Number generators
! An example of a uniform (Quasi)-random number generator
! Genering random number in range (0,1) using integer seed , idum is integer number.
***************************************************
	FUNCTION ran1(idum)
	INTEGER idum,ia,im,iq,ir,ntab,ndiv
	REAL*8 ran1,am,eps,rnmx
	PARAMETER (ia=16807,im=2147483647,am=1./im,
     &iq=127773,ir=2836,ntab=32,ndiv=1+(im-1)/ntab,
     &eps=1.2e-7,rnmx=1.-eps)
	INTEGER j,k,iv(ntab),iy
	SAVE iv,iy
	DATA iv /ntab*0/, iy /0/
	if (idum<=0.or.iy==0) then
	  idum=max(-idum,1)
	  do j=ntab+8,1,-1
	    k=idum/iq
	    idum=ia*(idum-k*iq)-ir*k
	    if (idum<0) idum=idum+im
	    if (j<=ntab) iv(j)=idum
	  end do
	  iy=iv(1)
	endif
	k=idum/iq
	idum=ia*(idum-k*iq)-ir*k
	if (idum<0) idum=idum+im
	j=1+iy/ndiv
	iy=iv(j)
	iv(j)=idum
	ran1=min(am*iy,rnmx)
	END FUNCTION ran1
******************************************************
************************************
*	Element π_ij(Z) of transition probability matrix
!Function to caculate the fraction 
*	\prod_{e=a,\ne c}^{b} \frac{\theta_k}{\bar{theta}_{e}-\bar{theta}_{c}}
************************************

	REAL*8 FUNCTION PROB(I,J,Z,THETA)
	USE VARIABLE;IMPLICIT NONE
	INTEGER,INTENT(IN):: I,J
	REAL*8,INTENT(IN):: Z,THETA(JMAX+1)
	REAL*8,EXTERNAL:: PROD2	
	!For the time being in this (later union necessary)
	REAL*8 RESERVE,HH
	INTEGER K
	PROB=0.0
	IF(Z <= 1.0D-06) RETURN
	RESERVE=1.0
	DO K=I,J-1
		RESERVE=RESERVE*THETA(K)
	END DO
	DO K=I,J
		PROB=PROB+DEXP(-THETA(K)*Z)/PROD2(I,J,K,THETA,JMAX)
	END DO
	PROB=RESERVE*PROB
	RETURN; END
!------------------------------------------------------
* Function to caculate the denominator in the fraction 
*	\prod_{e=a,\ne c}^{b} \frac{1}{\bar{theta}_{e}-\bar{theta}_{c}}
!------------------------------------------------------
	REAL*8 FUNCTION PROD2(A,B,C,THETA,JMAX)
	IMPLICIT NONE
	INTEGER,INTENT(IN):: A,B,C,JMAX
	REAL*8,INTENT(IN):: THETA(JMAX+1)
	INTEGER E

	PROD2 = 1.0
	DO E = A, B
		IF(E /= C) PROD2 = (THETA(E)-THETA(C)) * PROD2
	END DO
	RETURN;	END
***************************************************
*Random number generation from normal distribution
*  This function is for generating random number using two other number random generator
! ransu2, which intends to reduce noise in randome number.
! This method is called Box-Mueller algorithm
! Another name is called Gauss=random
! Reference is http://www.fortran.com/gauss_random.html
! Below is just a slightly modification of function double precision FUNCTION cgauss()
! in program provided in http://www.fortran.com/gauss_random.html
! Results from running this function is random number with Mean = MYU and standard deviation SIGMA
***************************************************
	REAL*8 FUNCTION RSEIKI(MYU,SIGMA)
	IMPLICIT NONE
	REAL*8,INTENT(IN):: MYU, SIGMA
	INTEGER ISET
	REAL*8 FAC,GSET,RSQ,V1,V2
	SAVE ISET, GSET
	DATA ISET/0/
	REAL*8,EXTERNAL:: RANSU2

	IF (ISET==0) THEN
		DO
			V1 = 2.0*RANSU2()-1.0;	V2 = 2.0*RANSU2()-1.0
			RSQ = V1**2+V2**2
			IF(.NOT.(RSQ>=1.0 .OR. RSQ==0.)) EXIT
		END DO
	
		FAC=DSQRT(-2.0*DLOG(RSQ)/RSQ)
		GSET = V1*FAC
		RSEIKI = V2*FAC
		ISET = 1
	ELSE
		RSEIKI = GSET
		ISET = 0
	END IF
	RSEIKI=RSEIKI*SIGMA+MYU
	!Print*, rseiki
	!stop
	RETURN;	END FUNCTION
**********************************************