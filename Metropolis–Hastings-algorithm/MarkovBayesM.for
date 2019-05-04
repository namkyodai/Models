*******************************************************
*     Bayesian Inference of Markov transition Probability(Using MH algorithm)
*	(With normal distribution)
*	Caculating how fast the obtained data can be convegent
*******************************************************
	MODULE FILE
		CHARACTER(LEN=30),PARAMETER:: INPUT_F1="input\data10.csv"
		CHARACTER(LEN=30),PARAMETER:: INPUT_F2="input\info1.csv"
		CHARACTER(LEN=30),PARAMETER:: OUTPUT_F1="output\point.csv"
		CHARACTER(LEN=30),PARAMETER:: OUTPUT_F2="output\result.csv"
	END MODULE
!------------------------------------------------------------------------
	MODULE DIM
		INTEGER,PARAMETER:: JMAX=4	!Condition states (J-1)
		INTEGER,PARAMETER:: MMAX=1	!Numbers of covariates
		INTEGER,PARAMETER:: JMJM=JMAX*MMAX
		INTEGER,PARAMETER:: NDATA=845	!Number of data
!	Conditions for caculation
		INTEGER,PARAMETER:: NPOINT=30000 !Desired sample numbers
		INTEGER,PARAMETER:: IBURN=2000 !Convergent time of MCMC
	!	1.Metropolis alrorithm, 2.Independent MH alrorigthm, 3.AR-MH alrorightm(QPATTERN=1,3 only)
		INTEGER,PARAMETER:: SMETHOD=3
	!	Method of assuming proposal distribution
	!	1.Taylor approximation, 2.Random walk process, 3. Arbitrary normal distribution, 4.ARS
		INTEGER,PARAMETER:: QPATTERN=3
	!	TRUE:Accept，FALSE:Reject
		LOGICAL,PARAMETER:: PREJECT=.TRUE.
	!	1.All prior information, 2.Partial prior information, 3.No prior information
		INTEGER,PARAMETER:: INFOTYPE=3
		INTEGER,PARAMETER:: LIMIT=1	!Interation of program do not have conve
		!gence when
	END MODULE
!---------------------------------------------
	MODULE VARIABLE
	USE DIM
		REAL*8 BETA(JMAX,MMAX)
		DATA BETA(1:4,1)/-1.1011, -1.0,-1.0882,-1.0/
!		DATA BETA(2,1:3)/-1.5546, 0.2392, 3.0293/
!		DATA BETA(3,1:3)/-1.9727, 0.6957, 0.0/
!		DATA BETA(4,1:3)/-2.4399, 0.845, 0.513/
!		DATA BETA(5,1:3)/-2.3233, 0.0, 0.0/
!		DATA BETA(6,1:3)/-1.951, 1.5439, 0.0/
		INTEGER,PARAMETER:: D1=0	!Number of removed parameters
		INTEGER DEL(D1)
		DATA DEL/2,9,14,15,18/	
		!sequence number of deleted parameter ((J-1)*MMAX+M)
		REAL*8 LAMB(JMAX,NDATA)
		REAL*8 MYU0(MMAX,JMAX),ROH(MMAX,MMAX,JMAX)
	END MODULE
!------------------------------------------------------------------------
	MODULE BOUND
	USE DIM
		LOGICAL OB(JMJM),UB(JMJM)	!Upper and lower limits of region D
		DATA OB/JMJM*.FALSE./
		DATA UB/JMJM*.FALSE./
		REAL*8 OX(JMJM),UX(JMJM)	!Lower limits of region D
		REAL*8 HTHETA(JMJM)	!Proposed mode of posterior distribution
		REAL*8 SVS(JMJM)
		DATA SVS(1:4)/1.0, 1.0,1.0,1.0/
!		DATA SVS(4:6)/1.0, 0.1, 0.1/
!		DATA SVS(7:9)/1.0, 0.1, 0.1/
!		DATA SVS(10:12)/1.0, 0.1, 0.1/
!		DATA SVS(13:15)/1.0, 0.1, 0.1/
!		DATA SVS(16:18)/1.0, 0.1, 0.1/
	!	Parameters that rules out ARS
		INTEGER,PARAMETER:: KMAX=2	!Numbers of initial standard coordinates
		INTEGER NT	!Tk：Maximum number of x coordinates
		REAL*8 SXK(JMAX*MMAX,KMAX)	!Starting point
	END MODULE

	MODULE DATABASE
	USE DIM
		REAL*8 XK(MMAX,NDATA), ZK(NDATA)
		INTEGER IK(NDATA),JK(NDATA)
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
	USE VARIABLE; USE FILE; USE BOUND; IMPLICIT NONE
	REAL*8 SAMP(0:NPOINT,JMJM)	!Sample collection
*	Random number	
	INTEGER IR1,IX,IY,IZ,TANE
	COMMON /IRAND1/IR1
	COMMON /IRAND2/IX,IY,IZ
	CHARACTER(LEN=30):: hour
	INTEGER I,J,M,ID

!-----Changing point-----------------------------------
!	Initial value of ARS
	DO I=1,JMJM
		SXK(I,1)=-10.0;	SXK(I,2)=8.0
	END DO
	NT=30
*	Initial value of MH
	SAMP=0.0
	DO J=1,JMAX; DO M=1,MMAX
		ID = (J-1)*MMAX+M
		SAMP(0,ID)=BETA(J,M)
		HTHETA(ID)=BETA(J,M)
	END DO; END DO
	SAMP=4.0*SAMP/5.0
!------------------------------------------------
	CALL TIME(hour); PRINT *, hour
	MYU0=0.0; ROH=0.0
*     Random number	
	CALL CLOCK (IR1,2)
	IX=16251; IY=21157; IZ=3583
	TANE=IR1
*	Data input from WHERE
	CALL INPUT()
	
	SELECT CASE(SMETHOD)
		CASE(1); PRINT *, "Metropolis"
		CASE(2); PRINT *, "Independent MH"
		CASE(3); PRINT *, "ARMH"
	END SELECT
	
	CALL MCMC_MH(SAMP)
	
!	Conversion to initial information
	IF(INFOTYPE/=3) CALL REVERSE(2)

	CALL OUTPUT(TANE,SAMP)
!	CALL BEEPQQ (3000, 500); CALL BEEPQQ (4000, 500)
	CALL TIME(hour); PRINT *, hour
	END PROGRAM
**************************************************************
*	MCMC method using MH alrorigth
**************************************************************
	SUBROUTINE MCMC_MH(SAMP)
	USE VARIABLE; IMPLICIT NONE
	REAL,EXTERNAL:: GEWEKEZ
	REAL*8,INTENT(OUT):: SAMP(0:NPOINT,JMJM)	!Sampling set
	REAL*8 SAMP1(NPOINT)		!Sub-set of sampling
	REAL Z(JMJM); COMMON /GEWEKE/Z
	INTEGER I,J,M,T,ID,IPOINT,IREV
	LOGICAL KENTEI

!*********************
*	Burn-in
!*********************
	IPOINT=1
	
	CALL SAMPLING_PROCESS(1,SAMP,IPOINT)
	DO ID=1,JMJM; SAMP(0,ID) = SAMP(IBURN,ID); END DO
	
!*********************
*	Sampling
!*********************
	IPOINT=1; IREV=0
	
 100	CONTINUE
	CALL SAMPLING_PROCESS(2,SAMP,IPOINT)

!*********************
*	Hypotheses testing
!*********************
	KENTEI=.TRUE.
	I=1
		DO J=1,JMAX; DO M=1,MMAX
			ID = (J-1)*MMAX+M
			DO T=1,NPOINT
				SAMP1(T)=SAMP(T,ID)
			END DO
	
!			IF(DEL(I)==ID) THEN
!				I=I+1; CYCLE
!			END IF

			Z(ID)= GEWEKEZ(SAMP1,NPOINT)
			IF(Z(ID) < 1.96) THEN	!５％ Significant level
				PRINT 90, ID,Z(ID), "converge"
			ELSE
				PRINT 90, ID,Z(ID), "Not converge"
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
	RETURN; END SUBROUTINE MCMC_MH
***************************************************
*	MH algorithm
***************************************************
	SUBROUTINE SAMPLING_PROCESS(ITYPE,SAMP,IPOINT)
	USE VARIABLE; IMPLICIT NONE
	INTEGER,INTENT(IN):: ITYPE
	REAL*8,INTENT(OUT):: SAMP(0:NPOINT,JMJM)	!Sample collection
	INTEGER,INTENT(IN):: IPOINT
	REAL*8 X	!Return value
	LOGICAL REJ
	INTEGER IREJ(JMJM),NTIM(JMJM); COMMON /ACEP/IREJ,NTIM
	CHARACTER(LEN=30):: hour
	INTEGER T,I,ID,J,M
	IREJ=0; NTIM=0
	CALL TIME(hour); PRINT *, hour
		SELECT CASE(ITYPE)
			CASE(1); PRINT *, "Burn-in period"
			CASE(2); PRINT *, "Sampling period"
		END SELECT
		DO T=IPOINT,NPOINT
			I=1
			IF(MOD(T,1000)==0.0) PRINT *, T,"Next time around"
			DO J=1,JMAX; DO M=1,MMAX
				ID = (J-1)*MMAX+M
	
!				IF(DEL(I)==ID) THEN
!					I=I+1; CYCLE
!				END IF
 90	CONTINUE
	
				IF(INFOTYPE/=3) CALL INITIALIZE(ID)
	
				CALL INDMH_SAMPLE(X,SAMP(T-1,ID),ID,REJ)
	
				NTIM(ID) = NTIM(ID) +1
				BETA(J,M)=X
				IF (REJ==.TRUE.) THEN
					IREJ(ID) = IREJ(ID)+1
					IF(PREJECT==.TRUE.) GOTO 90
				END IF
				SAMP(T,ID)=X
			END DO; END DO
		END DO
	END SUBROUTINE
*****************************************************
	SUBROUTINE INDMH_SAMPLE(NX,PX,ID,REJ)
	USE BOUND; IMPLICIT NONE
	REAL*8,EXTERNAL:: LNFUNC, ARMHTYPE, RANSU1, RANSU2, ARS
	INTEGER,INTENT(IN):: ID
	REAL*8,INTENT(IN):: PX	!Previous parameters
	LOGICAL REJ
	REAL*8 MM, VS 
	!Mean/Variance of proposal distribution (normal distribution)
	REAL*8 NX, CX
	REAL*8 DH,DH2,Q,CW,PW,R
	REAL*8 ALPHA
			
!	Making proposal distribution
	SELECT CASE(QPATTERN)
		CASE(1:3) ! Normal distribution
			SELECT CASE(QPATTERN)
				CASE(1)	!Approximation
	stop
					CALL MAKEDH(HTHETA(ID),DH,DH2,ID)
	
					VS = -1.0/DH2	!Distributed
					IF(VS < 0.0) VS = SVS(ID)
					MM = HTHETA(ID) + VS*DH	!AVERAGE
				CASE(2) !Randome walk process
					MM = 0.0;	VS = SVS(ID)
				CASE(3) !Other (initial average belived)
					MM = HTHETA(ID); VS = SVS(ID)
			END SELECT
			VS = DSQRT(VS)	!Standar deviation
		CASE DEFAULT
	END SELECT

 !	STEP2-1
 100	CONTINUE
 !	Candidate generation
	SELECT CASE(QPATTERN)
		CASE(1); CX = RSEIKI(MM,VS)
		CASE(2); CX = PX + RSEIKI(MM,VS)
		CASE(3); CX = RSEIKI(MM,VS)
		CASE(4); CX = ARS(ID)
	END SELECT
	
	IF(UB(ID)==.TRUE. .AND. UX(ID)<=0.0) GOTO 100
	IF(OB(ID)==.TRUE. .AND. OX(ID)>=0.0) GOTO 100
	IF(SMETHOD==3) THEN	!A-R Step
		IF(QPATTERN == 1 .OR. QPATTERN==3) THEN
			ALPHA = ARMHTYPE(CX,HTHETA(ID),DH,DH2,ID)
			R=RANSU1()
			IF(ALPHA < DLOG(R)) GOTO 100
		ELSE
			PRINT *, "まだ対応していません" ; STOP
		END IF
	END IF

!	STEP2-2
	CW = LNFUNC(CX,ID)
	PW = LNFUNC(PX,ID)
	IF(SMETHOD/=1 .OR. QPATTERN/=2) THEN
		CW = CW-LNQ(CX,MM,VS)
		PW = PW-LNQ(PX,MM,VS)
	END IF
	
	ALPHA = CW - PW
	REJ=.FALSE.
	IF(ALPHA > 0.0) THEN
		NX=CX
	ELSE; ALPHA=DEXP(ALPHA)
!		STEP2-3
		R=RANSU1()
		IF(ALPHA >= R) THEN; NX=CX
		ELSE;	NX=PX; REJ=.TRUE.
		END IF
	END IF
	RETURN
!--------------------
*	Internal procedures
!--------------------
	CONTAINS
!----------------------------------
*	Logarithm value of proposal distribution
!----------------------------------
	REAL*8 FUNCTION LNQ(X,MYU,SIGMA)
	REAL*8,INTENT(IN):: X,MYU,SIGMA
	LNQ = LNSEIKI(X,MYU,SIGMA)
	END FUNCTION
!----------------------------------
*	The value of one dimension normal distribution (logarithm value) is partially requested. 
!----------------------------------
	REAL*8 FUNCTION LNSEIKI(X,MYU,SIGMA)
	IMPLICIT NONE
	REAL*8,INTENT(IN):: X, MYU, SIGMA
	
	LNSEIKI=-(X-MYU)**2/(2.0*SIGMA**2)

	END FUNCTION
!----------------------------------
*     From a normal distribution random number generator (Box-Muller, improved version of Shibuya)
!----------------------------------
	REAL*8 FUNCTION RSEIKI(MYU,SIGMA)
	IMPLICIT NONE
	REAL*8,INTENT(IN):: MYU, SIGMA
	INTEGER ISET
	REAL*8 FAC,GSET,RSQ,V1,V2
	SAVE ISET, GSET
	DATA ISET/0/

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

	RETURN;	END FUNCTION
!---------------------------------------------------------------
*	End internal procedures
!---------------------------------------------------------------
	END SUBROUTINE
********************************************************
*	For the AR-MH (so far only proposed distribution of the normal distribution)
********************************************************
	REAL*8 FUNCTION ARMHTYPE(X,M,D1,D2,ID)
	IMPLICIT NONE
		REAL*8,EXTERNAL:: LNFUNC
		INTEGER,INTENT(IN):: ID
		REAL*8,INTENT(IN):: X,M,D1,D2
		REAL*8 NRES

		NRES = LNFUNC(M,ID) + D1*(X-M) + D2/2.0*(X-M)**2
		ARMHTYPE = LNFUNC(X,ID) - NRES
	RETURN; END FUNCTION
***************************************************************
*	Logarithm of the conditional posterior distribution (partial)
***************************************************************
	REAL*8 FUNCTION LNFUNC(X,ID)
	USE MTHETA;	USE VARIABLE; IMPLICIT NONE
	REAL*8,EXTERNAL:: LH, PROD1
	REAL*8,INTENT(IN):: X
	INTEGER,INTENT(IN):: ID
	REAL*8 XKK(MMAX,JMAX)
	REAL*8 MYU1
	REAL*8 RES1
	COMMON /HEIKIN/MYU1
	INTEGER J1,M1, MINTHETA, CHECK
	INTEGER I,J,K,M,L

	M1 = MOD(ID,MMAX); IF(M1==0) M1=MMAX
	J1 = (ID-M1)/MMAX + 1

	BETA(J1,M1)=X

!	Likelihood	
	LNFUNC = 0.0
	DO 100 K = 1, NDATA
		CALL ARRANGE(K,XKK,J1,CHECK)
		SELECT CASE(CHECK)
			CASE(1); CONTINUE
			CASE(2); CONTINUE
!				LNFUNC = LNFUNC + RES1; CYCLE
			CASE(3); CYCLE
		END SELECT

		CALL THETAVALUE(IKK,JKK,XKK,BETA,TSA,MINTHETA)
		DO L=IKK,JKK	!Minimum will be adjusted to zero Θ
			PSUM(L)=DEXP(TSA(MINTHETA,L)*ZKK)/PROD1(IKK,JKK,L,JMAX,TSA)
		END DO
		RES1 = LH(MINTHETA)
	
		IF (J1==JKK) THEN
			LNFUNC = LNFUNC + RES1
		ELSE
			RES1 = RES1 + BETA(J1,M1)*XKK(M1,J1)
			LNFUNC = LNFUNC + RES1
		END IF
 100	CONTINUE
	
!	Prior distribution	
	IF(INFOTYPE==3) RETURN

	IF(INFOTYPE==1 .OR. M1==1) THEN
		IF(INFOTYPE==2) MYU1 = MYU0(M1,J1)
		RES1 = BETA(J1,M1) - MYU1
		LNFUNC = LNFUNC - ROH(M1,M1,J1) * RES1**2/2.0D0
	END IF

	RETURN; END FUNCTION
********************************************************
*The first partial differential of posteriori distribution and logarithm value of partial differential twice (partial)	
********************************************************
	SUBROUTINE MAKEDH(X,DH,DH2,ID)
	USE MTHETA;	USE VARIABLE;  IMPLICIT NONE
	REAL*8,EXTERNAL:: PROD1,F_DIF,S_DIF
	REAL*8,INTENT(OUT):: DH, DH2
	REAL*8,INTENT(IN):: X
	INTEGER,INTENT(IN):: ID
	REAL*8 XKK(MMAX,JMAX)
	REAL*8 MYU1
	COMMON /HEIKIN/MYU1
	INTEGER MINTHETA, M1, J1, CHECK
	REAL*8 DF, DF2, RESF, RESG
	INTEGER I,J,K,M
	
	M1 = MOD(ID,MMAX); IF(M1==0) M1=MMAX
	J1 = (ID-M1)/MMAX + 1
	BETA(J1,M1)=X

	DH=0.0;	DH2=0.0; THETA=0.0
	DO 100 K=1,NDATA
	
		CALL ARRANGE(K,XKK,J1,CHECK)
		SELECT CASE(CHECK)
			CASE(1); CONTINUE
			CASE(2); CONTINUE
!				DH = DH + RESF
!				DH2= DH2+ RESG
!				CYCLE
			CASE(3); CYCLE
		END SELECT
	
		!Value of Θ
		CALL THETAVALUE(IKK,JKK,XKK,BETA,TSA,MINTHETA)
				
		PSUM=0.0
		DO J=IKK,JKK
			PSUM(J)=DEXP(-TSA(J,MINTHETA)*ZKK)/PROD1(IKK,JKK,J,JMAX,TSA)
		END DO
		!The first floor partial differential dlnπ_(ij)/dΘ_l of π is obtained. 
		DF = F_DIF(J1)*THETA(J1)*XKK(M1,J1)
		!The second partial differential
		RESG = DF*XKK(M1,J1)*THETA(M1)
		RESG = RESG * (XKK(M1,J1)-RESG)
		RESG = RESG + S_DIF(J1)*(XKK(M1,J1)*THETA(M1))**2

		RESF = DF*XKK(M1,J1)*THETA(M1)
		IF(J1/=JKK) RESF = RESF + XKK(M1,J1)

		DH = DH + RESF
		DH2 = DH2 + RESG
 100	CONTINUE

	IF(INFOTYPE==3) RETURN

	IF(INFOTYPE==1 .OR. M1==1) THEN
		IF(INFOTYPE==2) MYU1 = MYU0(M1,J1)
		RESF = BETA(J1,M1) - MYU1
		DH = DH - ROH(M1,M1,J1) * RESF
		DH2 = DH2 - ROH(M1,M1,J1)
	END IF
	END SUBROUTINE MAKEDH
**********************************************
*	Evaluation H(x) = Ln G(x) of function H(x) and H '(x)
**********************************************
	SUBROUTINE MAKEH(V,H,DH,ID)
	USE MTHETA;	USE VARIABLE; IMPLICIT NONE
	REAL*8,EXTERNAL:: LH, F_DIF, PROD1
	REAL*8,INTENT(OUT):: H,DH
	REAL*8,INTENT(IN):: V
	INTEGER,INTENT(IN):: ID
	INTEGER MINTHETA
	REAL*8 MYU1
	COMMON /HEIKIN/MYU1
	INTEGER J1,M1
	REAL*8 RES, RES1, RES2, RES3, RESH,RESDH
	REAL*8 XKK(MMAX,JMAX)
	INTEGER CHECK	
	INTEGER I,J,K,M,L
	
	M1 = MOD(ID,MMAX); IF(M1==0) M1=MMAX
	J1 = (ID-M1)/MMAX + 1

	BETA(J1,M1)=V
	H = 0.0; DH = 0.0
	DO 100 K = 1, NDATA
		THETA=0.0

		CALL ARRANGE(K,XKK,J1,CHECK)
		SELECT CASE(CHECK)
			CASE(1); CONTINUE
			CASE(2)
				H = H + RESH;DH= DH+RESDH
				CYCLE
			CASE(3); CYCLE
		END SELECT

		CALL THETAVALUE(IKK,JKK,XKK,BETA,TSA,MINTHETA)
				
		DO L=IKK,JKK
			PSUM(L) = DEXP(-TSA(L,MINTHETA)*ZKK)/PROD1(IKK,JKK,L,JMAX,TSA)
		END DO
		
		RESH = LH(MINTHETA); RESDH = F_DIF(J1)*XKK(M1,J1)*THETA(M1)
		IF(J1/=JKK) THEN
			RESH = RESH + BETA(J1,M1)*XKK(M1,J1)
			RESDH = RESDH + XKK(M1,J1)
		END IF
		H = H + RESH
		DH = DH + RESDH
 100	CONTINUE

	IF(INFOTYPE==3) RETURN

	IF(INFOTYPE==1 .OR. M1==1) THEN
		IF(INFOTYPE==2) MYU1 = MYU0(M1,J1)
		RES = BETA(J1,M1) - MYU1
		H = H - ROH(M1,M1,J1) * RES**2/2.0D0
		DH = DH - ROH(M1,M1,J1) * RES
	END IF

	RETURN; END SUBROUTINE MAKEH
!----------------------------------------------------------
*     Logarithm value lnπ_(IK,JK) of posteriori distribution	
!----------------------------------------------------------
	REAL*8 FUNCTION LH(MINTHETA)
	USE MTHETA; IMPLICIT NONE
	INTEGER,INTENT(IN):: MINTHETA
	REAL*8 RES
	INTEGER L

	RES = 0.0
	DO L= IKK,JKK; RES = RES + PSUM(L);	END DO
	
	IF(RES <= 0.0) THEN
		PRINT *, "WHY ?", "Logarithm value"; STOP
	ELSE
	! Taking the place TSA(MINTHETA,JMAX+1) of minimum Θ
		LH = DLOG(RES) - TSA(MINTHETA,JMAX+1)*ZKK
	END IF

	END FUNCTION
!----------------------------------------------------------
*	The first dlnπ_(IK,JK)/partial differential df(L) = dΘ_L
!----------------------------------------------------------
	REAL*8 FUNCTION F_DIF(MJ)
	USE MTHETA; IMPLICIT NONE
	INTEGER,INTENT(IN):: MJ
	INTEGER L
	REAL*8 SUM1,SUM2

	SUM1 = 0.0; SUM2=0.0
	DO L= IKK,JKK
		IF(L==MJ) CYCLE
		SUM1 = SUM1 + PSUM(L)/TSA(L,MJ); SUM2 = SUM2 + TSA(L,MJ)
	END DO
	SUM1 = SUM1 + (SUM2-ZKK)*PSUM(MJ)

	SUM2 = 0.0
	DO L= IKK,JKK; SUM2 = SUM2 + PSUM(L);	END DO
!	When the overflow is generated, here is corrected. 
	F_DIF = SUM1/SUM2
	END FUNCTION F_DIF
!----------------------------------------------------------
*	The second partial differential
!----------------------------------------------------------
	REAL*8 FUNCTION S_DIF(MJ)
	USE MTHETA; IMPLICIT NONE
	REAL*8,EXTERNAL:: PROD1
	INTEGER,INTENT(IN):: MJ
	INTEGER H,P
	REAL*8 SUM1,SUM2,SUM3,RES

	S_DIF=0.0
	SUM1=0.0; SUM2=0.0; SUM3=0.0
	DO H=IKK,JKK
		IF(H==MJ) CYCLE
		RES = 1.0/TSA(H,MJ)**2
		SUM1 = SUM1 + RES
		SUM2 = SUM2 + 1.0/TSA(H,MJ)
		SUM3 = SUM3 + PSUM(H)
		S_DIF = S_DIF + 2.0*RES*PSUM(H)
	END DO
	
	RES = SUM1 + (SUM2-ZKK)**2
	S_DIF = S_DIF + RES*PSUM(MJ)
	SUM3 = SUM3 + PSUM(MJ)
	
	S_DIF = S_DIF/SUM3

	END FUNCTION S_DIF
********************************************************
*	Calculate labor-saving settings
********************************************************
	SUBROUTINE ARRANGE(NK,XKK,MJ,CHECK)
	USE MTHETA; USE VARIABLE; USE DATABASE; IMPLICIT NONE
	INTEGER,INTENT(IN):: NK,MJ
	REAL*8,INTENT(OUT):: XKK(MMAX,JMAX)
	LOGICAL,INTENT(OUT):: CHECK
	INTEGER I,J,M,N

		IF(IK(NK) > MJ .OR. JK(NK) < MJ) THEN
			CHECK=3; RETURN
		END IF

		!As conditions deteriorated between the two, the inspection period, characteristic vector
		DO J=1,JMAX;DO M=1,MMAX
			XKK(M,J)=XK(M,NK)
		END DO; END DO

!		DO I=1,D1	! Removed from the parameter space
!			M = MOD(DEL(I)+(MMAX-1),MMAX) + 1; J = (DEL(I)-M)/MMAX+1
!			XKK(M,J)=0.0
!		END DO

		IF (IK(NK)==RESI .AND. JK(NK)==RESJ .AND. ZK(NK)==RESZ) THEN
			CHECK=2
			DO J=1,JMAX; DO M=1,MMAX
				IF(XKK(M,J) /= RESX(M,J)) THEN
					CHECK=1; GOTO 100
				END IF
			END DO; END DO

 100			CONTINUE
			IF (CHECK==2) RETURN
		END IF

		IKK=IK(NK); JKK=JK(NK); ZKK=ZK(NK)
		RESI=IKK; RESJ=JKK; RESZ=ZKK; RESX=XKK
	END SUBROUTINE
*******************************************
*	The value of Θi(i=IK,JK) is requested. 
*******************************************
	SUBROUTINE THETAVALUE(IK,JK,XKK,BETA,THETASA,MINTHETA)
	USE DIM; USE MTHETA; IMPLICIT NONE
	REAL*8,INTENT(IN):: XKK(MMAX,JMAX),BETA(JMAX,MMAX)
	REAL*8,INTENT(OUT):: THETASA(JMAX+1,JMAX+1)
	INTEGER,INTENT(IN):: IK,JK
	INTEGER,INTENT(OUT):: MINTHETA
	INTEGER I,J,M
	REAL*8 RESERVE

	!Calculation of Θ
	THETA=0.0
	DO J=IK,JK
		IF(J==JMAX+1) CYCLE
		DO M=1,MMAX
			THETA(J)=THETA(J) + BETA(J,M)*XKK(M,J)
		END DO
		THETA(J) = DEXP(THETA(J))
	END DO
	!Difference of Θi-Θj
	THETASA=0.0
	DO J=IK,JK; DO I=IK,JK
		THETASA(I,J)=THETA(I)-THETA(J)
	END DO;	END DO
		
	!	To prevent the overflow, minimum Θ is examined. 
	MINTHETA=IK; RESERVE=DEXP(600.0D0)
	DO J=IK,JK
		IF(J==JMAX+1) CYCLE
		IF(RESERVE > THETA(J)) THEN
			RESERVE=THETA(J); MINTHETA=J
		END IF
	END DO
	
	END SUBROUTINE
********************************************************
*	\prod_{e=a,\ne c}^{b} \frac{1}{\bar{theta}_{e,c}}
********************************************************
	REAL*8 FUNCTION PROD1(A,B,C,JMAX,THETASA)
	IMPLICIT NONE
	INTEGER,INTENT(IN):: A,B,C,JMAX
	REAL*8,INTENT(IN):: THETASA(JMAX+1,JMAX+1)
	INTEGER E

	PROD1 = 1.0
	DO E = A, B
		IF(E /= C) PROD1 = THETASA(E,C) * PROD1
	END DO
	RETURN;	END
***********************************************
	SUBROUTINE INITIALIZE(ID)
***********************************************
	USE VARIABLE; IMPLICIT NONE
	INTEGER,INTENT(IN):: ID
	REAL*8 MYU1
	COMMON /HEIKIN/MYU1
	INTEGER M,M1,J1

	M1 = MOD(ID,MMAX); IF(M1==0) M1=MMAX
	J1 = (ID-M1)/MMAX + 1

	MYU1=0.0
	DO M=1,MMAX
		IF(M/=M1) MYU1 = MYU1 + (BETA(J1,M)-MYU0(M,J1))*ROH(M1,M,J1)
	END DO
	MYU1=MYU0(M1,J1) - MYU1/ROH(M1,M1,J1)
	RETURN; END SUBROUTINE
***********************************************************************
*	Input from input file
***********************************************************************
	SUBROUTINE INPUT()
	USE DATABASE; USE VARIABLE; USE FILE; IMPLICIT NONE
	REAL*8 MAXX(MMAX)
	COMMON /ARR/MAXX
	CHARACTER CH
	INTEGER I,J,K,M,M1,M2

	MAXX=0.0; MAXX(1)=1.0
	XK=1.0
	OPEN(1000,file=INPUT_F1)
		READ(1000,*) CH
		DO K=1,NDATA
			READ(1000,*) IK(K),JK(K),ZK(K),(XK(M1,K),M1=2,MMAX)
			DO M1=2,MMAX; IF(XK(M1,K) > MAXX(M1)) MAXX(M1)=XK(M1,K); END DO

		END DO
	CLOSE(1000)
!		Adjustment Scale
	DO K=1,NDATA
		DO M1=2,MMAX; XK(M1,K)=XK(M1,K)/MAXX(M1); END DO
	END DO
	XK=XK

	IF(INFOTYPE==3) RETURN	!If no information

	OPEN(1001,file=INPUT_F2)
		DO J=1,JMAX
			READ(1001,*) CH
			READ(1001,*) CH,(MYU0(M1,J),M1=1,MMAX)
			DO M1=1,MMAX
				READ(1001,*) CH,(ROH(M1,M2,J),M2=1,MMAX)
			END DO
		END DO
	CLOSE(1001)

	PRINT *, "Prior information [average]"
	DO J=1,JMAX
		PRINT 30, J,(MYU0(M1,J),M1=1,MMAX)
	END DO
30	FORMAT(1X,I2,6F10.4)

!		It is conversion and a scale adjustment to the matrix inverse. 
	CALL REVERSE(1)

	RETURN;	END SUBROUTINE INPUT
***************************************************
*	Output
***************************************************
	SUBROUTINE OUTPUT(TANE,SAMP)
	USE VARIABLE; USE BOUND; USE FILE; IMPLICIT NONE
	REAL,EXTERNAL:: GEWEKEZ
	INTEGER,INTENT(IN):: TANE
	REAL*8,INTENT(IN):: SAMP(0:NPOINT,JMJM)	!Sampling set
	REAL*8 SAMP1(NPOINT)	!Consolidating
	REAL*8 MYU(MMAX,JMAX),SS(MMAX,MMAX,JMAX),SM	!Mean, Variance
	REAL*8 UBAR(JMAX*MMAX),OBAR(JMAX*MMAX)	!Confidence
	REAL Z(JMAX*MMAX); COMMON /GEWEKE/Z	!Hypotheses test GEWEKE
	LOGICAL REJ
	INTEGER IREJ(JMJM),NTIM(JMJM); COMMON /ACEP/IREJ,NTIM
	INTEGER I,J,K,M1,M2
	
!--------------------------
*	Output sampling
!--------------------------
	OPEN(2000,file=OUTPUT_F1)
		WRITE(2000,*) "Randome seed",",",TANE
		WRITE(2000,*) "Samples",",",NPOINT
		WRITE(2000,*) "All trial frequency"
		WRITE(2000,21)  (NTIM(I),I=1,JMJM)
		WRITE(2000,*) "Reject frequency"
		WRITE(2000,21) (IREJ(I),I=1,JMJM)
		WRITE(2000,*) "State 1",",,,","State ",",,,","State 3",",,,","State 4"
		WRITE(2000,*) ("B1",",","B2",",","B3",",",J=1,4)
		DO K=1,NPOINT
			WRITE(2000,20) (SAMP(K,J),J=1,JMAX*MMAX)
		END DO
	CLOSE(2000)
 20	FORMAT(1X,30(E15.8,","))
 21	FORMAT(1X,30(I7,","))
!--------------------------
*	Output Statistics
!--------------------------
!	Sample mean
	MYU=0.0
	DO J=1,JMAX; DO M1=1,MMAX
		DO K=1,NPOINT
			MYU(M1,J)=MYU(M1,J)+SAMP(K,(J-1)*MMAX+M1)
		END DO
	END DO;	END DO	
	MYU=MYU/NPOINT

!	Sample variance
	DO J=1, JMAX; DO M2=1,MMAX; DO M1=1,MMAX
		I = (J-1)*MMAX
		DO K=1,NPOINT
			SS(M1,M2,J)=SS(M1,M2,J)
     &			+(SAMP(K,I+M1)-MYU(M1,J))*(SAMP(K,I+M2)-MYU(M2,J))
		END DO
	END DO; END DO; END DO
	SS=SS/NPOINT

!	Confidence	
	DO J=1,JMAX*MMAX
		DO K=1,NPOINT; SAMP1(K)=SAMP(K,J);	END DO
		CALL BUBBLE(SAMP1,NPOINT)
		UBAR(J)=SAMP1(NPOINT*5/100)
		OBAR(J)=SAMP1(NPOINT*95/100)
	END DO
 90	FORMAT(1X, F10.2, A12)

!	Output file
	OPEN(3000,file=OUTPUT_F2)
		SELECT CASE(INFOTYPE)
			CASE(1); WRITE(3000,*) "All prior information"
			CASE(2); WRITE(3000,*) "Partial prior information"
			CASE(3); WRITE(3000,*) "No prior information"
		END SELECT
		SELECT CASE(SMETHOD)
			CASE(1); WRITE(3000,*) "Sampling methods",",","Metropolis"
			CASE(2); WRITE(3000,*) "Sampling methods",",","independent MH"
			CASE(3); WRITE(3000,*) "Sampling methods",",","ARMH"
		END SELECT
	

		SELECT CASE(QPATTERN)
	CASE(1); WRITE(3000,*) "Proposed distribution",",","Taylor"
	CASE(2); WRITE(3000,*) "Proposed distribution",",","Randomwalk"
		CASE(3); WRITE(3000,*) "Proposed distribution",",","Regular proper"
			CASE DEFAULT; WRITE(3000,*) "Refusal frequency"
		END SELECT

		IF(INFOTYPE/=3) THEN
			WRITE(3000,*) "Initial information"
			WRITE(3000,*) 
     &		",", "state1",",,,","state2",",,,","state3",",,,","state4",",,,"
			WRITE(3000,10) "Average initial",((MYU0(M1,J),M1=1,MMAX),J=1,JMAX)
			DO M1=1,MMAX
	WRITE(3000,10) "Initial dispersion",((ROH(M1,M2,J),M2=1,MMAX),J=1,JMAX)
			END DO
		END IF
				
		WRITE(3000,*) "Post information"
		WRITE(3000,*) 
     &		",", "state1",",,,","state2",",,,","state3",",,,","state4",",,,"
		WRITE(3000,10) "sample mean",((MYU(M1,J),M1=1,MMAX),J=1,JMAX)
		DO M1=1,MMAX
			WRITE(3000,10) "sample variance",((SS(M1,M2,J),M2=1,MMAX),J=1,JMAX)
		END DO

		WRITE(3000,10) "confidence",(UBAR(J),J=1,JMAX*MMAX)
		WRITE(3000,10) "confidence",(OBAR(J),J=1,JMAX*MMAX)
		WRITE(3000,10) "convergence creteria",(Z(J),J=1,JMAX*MMAX)
 10	FORMAT(1X,A8,",",20(E15.8,","))
	CLOSE(3000)

	RETURN;	END SUBROUTINE OUTPUT
**********************************************
*	Scale adjustment
**********************************************
	SUBROUTINE REVERSE(IP)
	USE VARIABLE; IMPLICIT NONE
	INTEGER,INTENT(IN):: IP
	REAL*8 A(MMAX,MMAX)
	REAL*8 MAXX(MMAX)	!The maximum value of explanatory variable
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
	ELSE	! Reverse Adjustment
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
*     Stationary distribution of convergence creteria (GEWEKE test hypotheses)
*******************************************************
	REAL FUNCTION GEWEKEZ(SAMP,NPOINT)
	IMPLICIT NONE
	REAL*8,EXTERNAL:: TVORA
	INTEGER,INTENT(IN):: NPOINT ! number of samples
	REAL*8,INTENT(IN):: SAMP(NPOINT)	!Sample sets (specimem collection)
	INTEGER NP1, NP2
	REAL*8 SAMP1(NPOINT/10), SAMP2(NPOINT/50)
	REAL*8 MYU1,MYU2,S1,S2
	REAL*8 Z	! statistic test
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
*	Calculation of sample variance (time series)
! Reference paper
! A Simple, Positive Semi-Definite, Heteroskedasticity
!	and Autocorrelation Consistent Covariance Matrix
*****************************************
	REAL*8 FUNCTION TVORA(SAMPX,MYU,NN)
	IMPLICIT NONE
	INTEGER,INTENT(IN):: NN
	REAL*8,INTENT(IN):: SAMPX(NN),MYU
	INTEGER I,J,M

		M = 10	! Bandwith
		TVORA= OMEGA(0)
		DO J=1, M
			TVORA = TVORA + (1.0 - DBLE(J)/(M+1.0))*2.0*OMEGA(J)
		END DO
	CONTAINS
!--------------------
*	Internal procedures
!--------------------
		REAL*8 FUNCTION OMEGA(S)
		IMPLICIT NONE
		INTEGER,INTENT(IN):: S
		INTEGER K
			OMEGA=0.0
			DO K = S + 1, NN
				OMEGA=OMEGA +(SAMPX(K)-MYU)*(SAMPX(K-S)-MYU)
			END DO
			OMEGA = OMEGA / NN
		END FUNCTION
!---------------------------------------------------------
	END FUNCTION TVORA
*******************************************************
*	Adaptive rejection sampling
*******************************************************
	REAL*8 FUNCTION ARS(ID)
	USE BOUND; IMPLICIT NONE
		INTEGER,INTENT(IN):: ID
		INTEGER IFAULT,I,K
		REAL*8 XX(KMAX)	!Initial coordinates
		REAL*8 RES

		DO K=1,KMAX; XX(K)=SXK(ID,K); END DO
		DO I=1,1000
			IFAULT=0
			CALL REJECT_SAMPLE(ARS,NT,KMAX, 
     &			UB(ID),OB(ID),UX(ID),OX(ID),XX,ID,IFAULT)
			SELECT CASE(IFAULT)
				CASE(0); EXIT
				CASE(1); XX(1)=XX(1)-(DABS(XX(KMAX))-DABS(XX(1)))/10.0
				CASE(2); XX(KMAX)=XX(KMAX)+(DABS(XX(KMAX))-DABS(XX(1)))/10.0
				CASE(3)
					RES=(DABS(XX(KMAX))-DABS(XX(1)))/10.0
					XX(1)=XX(1)+RES; XX(KMAX)=XX(KMAX)-RES
				CASE(6)
					RES=(DABS(XX(KMAX))-DABS(XX(1)))/10.0
					XX(1)=XX(1)+RES; XX(KMAX)=XX(KMAX)-RES
			END SELECT
			IF(I==1000) THEN
				PRINT *, "WHY?"; STOP
			END IF
		END DO
		SXK(ID,1)=XX(1); SXK(ID,KMAX)=XX(KMAX)

	END FUNCTION
!----------------------------------------------------------------------
	SUBROUTINE REJECT_SAMPLE(XPOINT,NT,NK,UB,OB,UX,OX,XX,NO,IFAULT)
	IMPLICIT NONE
	INTEGER,INTENT(IN):: NO	
	!Which conditional probability distribution is generated?
	INTEGER,INTENT(IN):: NK
	INTEGER IFAULT	!Error adjustment
	REAL*8,EXTERNAL:: RANSU1
	REAL*8,PARAMETER:: EPS10=1.0D-16	!ε
	REAL*8 XPOINT, XX(NK)
	INTEGER,INTENT(IN):: NT	!maximum number of coordinates x
	LOGICAL,INTENT(IN):: OB, UB	!upper and lower limits of the region D
	REAL*8,INTENT(IN):: UX,OX	!Lower limit
	REAL*8 HK(NT),DHK(NT),ZK(0:NT),HZK(0:NT),XKK(NT)
	INTEGER,PARAMETER:: IEMAX=700	!The maximum argument of EXP
	INTEGER,PARAMETER:: NIT=1000	!Number of cycles
	REAL*8 RX	!Generation coordinates
	REAL*8 W	!Uniform random
	REAL*8 LK,UK,HX,RH,RDH,RZ1,RZ2,RZH1,RZH2
	INTEGER RJ	!Insert location
	INTEGER DR	!TEST results
	REAL*8 HZMAX	!scale adjustment of H 
	REAL*8 RESERVE
	INTEGER K
	INTEGER I,IT

!---------------------------
!	Initialization step
!---------------------------
	K=NK
	DO I=1,K
		XKK(I)=XX(I)	!Initial standard coordinates Tk
		XPOINT=XKK(I)
		CALL MAKEH(XPOINT,HK(I),DHK(I),NO)
	END DO

	IF(UB==.FALSE. .AND. DHK(1) < 0.0) THEN
		PRINT *, "ID=",NO,"Inclination",DHK(1),"Left the range more widely"
		IFAULT=1; RETURN
	END IF

	IF(OB==.FALSE. .AND. DHK(K) > 0.0) THEN
		PRINT *, "ID=",NO,"inclination",DHK(K),"Left the range more widely"
		IFAULT=2; RETURN
	END IF
	
	HZK(0)=HK(1); HZK(K)=HK(K)
	IF (UB==.TRUE.) THEN
		ZK(0)=UX
		HZK(0)=HK(1)+DHK(1)*(UX-XKK(1))
	END IF
	IF (OB==.TRUE.) THEN
		ZK(K)=OX
		HZK(K)=HK(K)+DHK(K)*(OX-XKK(K))
	END IF

	DO I=1,K-1
		CALL MAKEZ(ZK(I),HZK(I),XKK(I),XKK(I+1),HK(I),HK(I+1),DHK(I),DHK(I+1))
	END DO

	DO 100 IT=1,NIT
		IF(UB==.TRUE.) THEN; HZMAX=DMAX1(HZK(0),HZK(1))
		ELSE; HZMAX=HZK(1)
		END IF
		DO I=2,K-1
			HZMAX=DMAX1(HZMAX,HZK(I))
		END DO
		IF(OB==.TRUE.) HZMAX=DMAX1(HZMAX,HZK(K))
!---------------------------
!	Sampling step
!---------------------------
!	Randomness from Sk
		CALL SKSAMPLE(RX,RJ,K,IFAULT)
		IF(IFAULT==6) RETURN
!Generating uniform random	
		W=RANSU1()

		UK=U(RX,XKK(RJ),HK(RJ),DHK(RJ))
		IF(RX<XKK(RJ)) RJ=RJ-1
		IF(RJ==0 .OR. RJ==K) THEN
			LK=-1.0D20	!Infinity
		ELSE
			LK=L(RX,XKK(RJ),XKK(RJ+1),HK(RJ),HK(RJ+1))
		END IF

!	SQUEEZING TEST
		IF(W <= DEXP(LK-UK)) THEN
			DR=0
!	REJECTION TEST
		ELSE
			CALL MAKEH(RX,RH,RDH,NO) 
			IF(W <= DEXP(RH-UK)) THEN;	DR=1
			ELSE; DR=2
			END IF
		END IF
!---------------------------
!	Updating step
!---------------------------
		IF(DR==0 .OR. DR==1) THEN
			XPOINT=RX; RETURN
		END IF
		IF(DR==2) THEN
			IF(RJ/=0) THEN
				CALL MAKEZ(RZ1,RZH1,XKK(RJ),RX,HK(RJ),RH,DHK(RJ),RDH)
			END IF
			IF(RJ/=K) THEN
     				CALL MAKEZ(RZ2,RZH2,RX,XKK(RJ+1),RH,HK(RJ+1),RDH,DHK(RJ+1))
			END IF
			IF(K==NT) THEN
				PRINT *, "ID=",NO, "Tk is full"; IFAULT=3; RETURN
			END IF
!	Tk,Zj is sorted in ascending order
			DO I=K,RJ+1,-1
				XKK(I+1)=XKK(I); ZK(I+1)=ZK(I); HZK(I+1)=HZK(I)
				HK(I+1)=HK(I); DHK(I+1)=DHK(I)
			END DO
			XKK(RJ+1)=RX; HK(RJ+1)=RH; DHK(RJ+1)=RDH
			IF(RJ==0) THEN
				ZK(0)=UX
				HZK(0)=HK(1)+DHK(1)*(UX-XKK(1))
				HZMAX=DMAX1(HZMAX,HZK(0))
			ELSE; ZK(RJ)=RZ1; HZK(RJ)=RZH1; HZMAX=DMAX1(HZMAX,RZH1)
			END IF
			IF(RJ==K) THEN
				ZK(K+1)=OX
				HZK(K+1)=HK(K)+DHK(K)*(OX-XKK(K))
				HZMAX=DMAX1(HZMAX,HZK(K+1))
			ELSE; ZK(RJ+1)=RZ2; HZK(RJ+1)=RZH2; HZMAX=DMAX1(HZMAX,RZH2)
			END IF
			K=K+1
		END IF
 100	CONTINUE

	PRINT *, "Necessary numbers of sampling"
	RETURN
!---------------------------------------------------------
*	Internal procedure
!---------------------------------------------------------
	CONTAINS
	SUBROUTINE MAKEZ(Z,ZH,X1,X2,H1,H2,DH1,DH2)
		REAL*8,INTENT(IN):: X1,X2,H1,H2,DH1,DH2
		REAL*8,INTENT(OUT):: Z,ZH
		REAL*8 D

		D=DH2-DH1
		IF(DABS(D) <= EPS10) THEN
			Z=0.5*(X1+X2)
			ZH=0.5*(H1+H2)
		ELSE IF(DABS(H1) < DABS(H2)) THEN
			Z=X2+(H1-H2+DH1*(X2-X1))/D
			ZH=DH1*(Z-X1)+H1
		ELSE
			Z=X1+(H1-H2+DH2*(X2-X1))/D
			ZH=DH2*(Z-X2)+H2
		END IF
		RETURN
	END SUBROUTINE
!------------
	REAL*8 FUNCTION U(X,XJ,H,DH)	!UPPER HULL
		REAL*8,INTENT(IN):: X,XJ,H,DH
		
		U=H+(X-XJ)*DH; RETURN
	END FUNCTION
!------------
	REAL*8 FUNCTION L(X,XJ,XJ1,H,H1)	!LOWER HULL
		REAL*8,INTENT(IN):: X,XJ,XJ1,H,H1
		
		L=H1+(X-XJ1)*(H1-H)/(XJ1-XJ); RETURN
	END FUNCTION
**************************************************
*	Randome number genration from Sk
!		Refer to the program of Gilk
**************************************************
	SUBROUTINE SKSAMPLE(X,JJ,KMAX,IFAULT)
	REAL*8,EXTERNAL:: RANSU2
!	Values to be estimated
	INTEGER,INTENT(OUT):: JJ	!J:Insert location
	REAL*8,INTENT(OUT):: X	!X:Sample coordinate
	INTEGER,INTENT(IN):: KMAX
	INTEGER,INTENT(OUT):: IFAULT
	REAL*8 SSUM(KMAX)	!Total areas under Uk
	REAL*8 LSUM
	REAL*8 EH,SIGN	!Value for storage
	LOGICAL HORIZ
	REAL*8 R	!Random
	INTEGER J

!	Uk: Caculating areas under each lines
	CALL UAREA(KMAX,SSUM)
	IF(SSUM(KMAX)<=0.0) THEN
		IFAULT=6; PRINT *, IFAULT,"Area below 0"; RETURN
	END IF
	LSUM=DLOG(SSUM(KMAX))
!	Transform the probability distribution
	DO J=1,KMAX-1
		R=SSUM(J)/SSUM(KMAX)
		SSUM(J)=R
	END DO
	SSUM(KMAX)=1.0

!	(0,1) will generate a uniform random number, determine the Area J
	R=RANSU2()
	IF(R==0.0) THEN
		PRINT *, "Random number is 0 to 2"; STOP
	END IF
	DO J=1,KMAX
		IF(SSUM(J) > R) EXIT	!エリアj(Zj-1 < X* < Zj)
	END DO
	JJ=J

	IF(JJ /= 1) R=R-SSUM(JJ-1)
		
	IF(JJ==1 .AND. UB==.FALSE.) THEN
		X=(DLOG(DHK(J)*R)+LSUM-HK(J)+XKK(J)*DHK(J)+HZMAX)/DHK(J)
	ELSE
		J=JJ
		EH=HZK(J-1)-HZMAX-LSUM
		HORIZ=(DABS(DHK(J)) < EPS10)
		IF(HORIZ==.TRUE.) THEN
			X=ZK(J-1)+R*DEXP(-EH)
		ELSE
			SIGN=DABS(DHK(J))/DHK(J)
			EH=DLOG(DABS(DHK(J))*R)-EH
			IF(EH < IEMAX) THEN
				X=ZK(J-1)+DLOG(1.0+SIGN*DEXP(EH))/DHK(J)
			ELSE
				X=ZK(J-1)+EH/DHK(J)
			END IF
		END IF
	END IF
	END SUBROUTINE SKSAMPLE
!-----------------------------------------------------------------------------
	SUBROUTINE UAREA(KMAX,SSUM)	! Uk: caculating areas under each line
	INTEGER,INTENT(IN):: KMAX
	REAL*8,INTENT(OUT):: SSUM(KMAX)	!Total areas under Uk
	REAL*8 SUM1
	INTEGER J
	LOGICAL HORIZ
	REAL*8 D
	
	J=1
	HORIZ=(DABS(DHK(J)) < EPS10)
	IF(UB==.FALSE. .AND. HORIZ==.FALSE.) THEN
		SSUM(1)=DEXP(HZK(J)-HZMAX)/DHK(J)
	ELSE IF(UB==.TRUE. .AND. HORIZ==.TRUE.) THEN
		SSUM(1)=(ZK(J)-UX)*DEXP(HZK(0)-HZMAX)
	ELSE IF(UB==.TRUE. .AND. HORIZ==.FALSE.) THEN
		D=HZK(0)-HZK(J)
		IF(D > IEMAX) THEN;	SSUM(1)=-DEXP(HZK(0)-HZMAX)/DHK(J)
		ELSE;	SSUM(1)=DEXP(HZK(J)-HZMAX)*(1.0-DEXP(D))/DHK(J)
		END IF
	ELSE
		SSUM(1)=0.0
	END IF

	DO J=2,KMAX-1
		HORIZ=(DABS(DHK(J)) < EPS10)
		D=HZK(J-1)-HZK(J)
		IF(HORIZ==.TRUE.) THEN
			SUM1=(ZK(J)-ZK(J-1))*DEXP((HZK(J)+HZK(J-1))*0.5-HZMAX)
		ELSE
			IF(D < IEMAX) THEN;	SUM1=DEXP(HZK(J)-HZMAX)*(1.0-DEXP(D))/DHK(J)
			ELSE; SUM1=-DEXP(HZK(J-1)-HZMAX)/DHK(J)
			END IF
		END IF
		SSUM(J)=SSUM(J-1)+SUM1
	END DO

	J=KMAX
	HORIZ=(DABS(DHK(J)) < EPS10)
	IF(OB==.FALSE. .AND. HORIZ==.FALSE.) THEN
		SUM1=-DEXP(HZK(J-1)-HZMAX)/DHK(J)
	ELSE IF(OB==.TRUE. .AND. HORIZ==.TRUE.) THEN
		SUM1=(OX-XKK(J))*DEXP((HZK(KMAX)+HK(J))*0.5-HZMAX)
	ELSE IF(OB==.TRUE. .AND. HORIZ==.FALSE.) THEN
		D=HZK(J-1)-HZK(KMAX)
		IF(D > IEMAX) THEN;	SUM1=-DEXP(HZK(J-1)-HZMAX)/DHK(J)
		ELSE; SUM1=DEXP(HZK(KMAX)-HZMAX)*(1.0-DEXP(D))/DHK(J)
		END IF
	END IF
	SSUM(KMAX)=SSUM(KMAX-1)+SUM1

	END SUBROUTINE
!---------------------------------------------------------
	END SUBROUTINE REJECT_SAMPLE
***************************************************
*Generation of pseudouniform random number (multiplying joint method)
***************************************************
	REAL*8 FUNCTION RANSU1()
	INTEGER IR1
	COMMON /IRAND1/IR1
		IR1=IAND(48828125*IR1, 2147483647)
		RANSU1=DBLE(IR1)/DBLE(2147483647)
	END FUNCTION
***************************************************
*	The generation of a pseudouniform random number (Multi multiplication joint method: AS183) The cycle is long. 
***************************************************
	REAL*8 FUNCTION RANSU2()
!It is necessary to give both of the seed of random numbers within the range from 1 to 30000. 
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
*Bubble sorting
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
*Matrix inverse (complete pivot method)
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
!(N-K+1)It searches for the one with the maximum absolute value from among the small next procession. 
		AW=DABS(A(K,K))
		P=K; Q=K
		DO J=K,N; DO I=K,N
			IF(AW < DABS(A(I,J))) THEN
				AW=DABS(A(I,J)); P=I; Q=J
			END IF
		END DO; END DO
		IF(AW < EPS) THEN
			PRINT *, "この行列特異" ;STOP
		END IF
!The one with the maximum absolute value is brought on the left. 
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
!	Operation to K line
		AW=1.0/A(K,K);	A(K,K)=1.0
		DO J=1,N; A(K,J)=A(K,J)*AW;	END DO
!	Operation excluding K line
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