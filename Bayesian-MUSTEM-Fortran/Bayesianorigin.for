*******************************************************
*     マルコフ推移確率のベイズ推定(Gibbsサンプリング)
*	(正規分布)
*******************************************************
	MODULE FILE
	CHARACTER(LEN=30),PARAMETER:: INPUT_F1="input\data.csv"
		CHARACTER(LEN=30),PARAMETER:: INPUT_F2="input\noinfo.csv"
		CHARACTER(LEN=30),PARAMETER:: OUTPUT_F1="output\point.csv"
		CHARACTER(LEN=30),PARAMETER:: OUTPUT_F2="output\result.csv"
	END MODULE
!------------------------------------------------------------------------
	MODULE VARIABLE
		INTEGER,PARAMETER:: NDATA=880	!データ数
		INTEGER,PARAMETER:: JMAX=6	!状態ランク数(J-1)
		INTEGER,PARAMETER:: MMAX=3	!説明変数の数
		REAL*8 BETA(JMAX,MMAX)
		INTEGER,PARAMETER:: D1=2	!削除パラメータ数
		INTEGER DEL(D1)
		DATA DEL/2,6/	!削除するパラメータ番号((J-1)*MMAX+M)
		INTEGER IK(NDATA),JK(NDATA)
		REAL*8 XK(MMAX,NDATA), ZK(NDATA), LAMB(JMAX,NDATA)
		REAL*8 MYU0(MMAX,JMAX),ROH(MMAX,MMAX,JMAX)
		INTEGER PATTERN
		INTEGER,PARAMETER:: LIMIT=1
	END MODULE
**********************************************
*	Gibbsサンプリングに必要な値
**********************************************
	SUBROUTINE FIRST(IBURN,KMAX,NT,SXK,JM,MM)
	INTEGER,INTENT(OUT):: IBURN	!マルコフチェーン収束期間
*	ARSを支配するパラメータ
	INTEGER,INTENT(IN):: KMAX,JM,MM	!初期基準座標数
	INTEGER,INTENT(OUT):: NT	!Tk：x座標の最大数
	REAL*8,INTENT(OUT):: SXK(JM*MM,KMAX)	!スタートポイント
	INTEGER I,K

		IBURN=3000;	NT=30
	DO I=1,JM*MM
		SXK(I,1)=-20.0;	SXK(I,2)=20.0
	END DO
!	SXK(3,1)=-20.0;	SXK(3,2)=35.0
!	SXK(6,1)=-30.0;		SXK(6,2)=20.0
!	SXK(11,1)=-30.0;	SXK(11,2)=15.0
!	SXK(12,1)=-30.0;	SXK(12,2)=5.0
!	SXK(16,1)=-300.0;	SXK(16,2)=10.0
!	SXK(17,1)=-20.0;	SXK(17,2)=50.0
	END SUBROUTINE
***************************************************************
*	MAIN PROGRAM
***************************************************************
	PROGRAM BAYES_MARKOV

	USE VARIABLE; USE FILE; IMPLICIT NONE
	INTEGER,PARAMETER:: NPOINT=10000 !希望サンプリング
	REAL*8 SAMP(NPOINT,JMAX*MMAX)	!サンプリング集合
!---------------------------------------------------------
*	乱数の種
	INTEGER IR1,IX,IY,IZ,TANE
	COMMON /IRAND1/IR1
	COMMON /IRAND2/IX,IY,IZ
!---------------------------------------------------------
	CHARACTER(LEN=30):: hour
	INTEGER I

	CALL TIME(hour); PRINT *, hour
	MYU0=0.0; ROH=0.0

*	乱数の種
	CALL CLOCK (IR1,2)
	IX=16251; IY=21157; IZ=3583
	TANE=IR1

*	ギブスサンプリング初期値
	BETA(1,1)=-1.0	; BETA(1,2)=0.0	;	BETA(1,3)=0.015
	BETA(2,1)=-2.0	; BETA(2,2)=0.0	;	BETA(2,3)=0.016
	BETA(3,1)=-2.1	; BETA(3,2)=0.0	;	BETA(3,3)=0.0
	BETA(4,1)=-1.1	; BETA(4,2)=0.0	;	BETA(4,3)=0.018
	BETA(5,1)=-1.2	; BETA(5,2)=0.0	;	BETA(5,3)=0.0
	BETA(6,1)=-1.5	; BETA(6,2)=0.0	;	BETA(6,3)=0.0

	DO I=1,1000
	PRINT *,"全事前情報(1) or 部分的事前情報(2) or 無情報(3)"
		READ *, PATTERN
		SELECT CASE(PATTERN)
			CASE(1:3); EXIT
		END SELECT
	END DO
!	PATTERN=1
*	データファイルからの入力
	CALL INPUT()
	CALL MCMC_G(SAMP,NPOINT)
	
!	初期情報への変換
	IF(PATTERN/=3) CALL REVERSE(2)

	CALL OUTPUT(TANE,SAMP,NPOINT)
!	CALL BEEPQQ (3000, 500); CALL BEEPQQ (4000, 500)
	CALL TIME(hour); PRINT *, hour
	STOP;	END PROGRAM
**********************************************
*	関数H(x)，H'(x)の評価  H(x)=Ln G(x)
**********************************************
	SUBROUTINE MAKEH(V,H,DH,ID)
	USE VARIABLE; IMPLICIT NONE
	REAL*8,EXTERNAL:: PROD1
	REAL*8,INTENT(OUT):: H,DH
	REAL*8,INTENT(IN):: V
	INTEGER,INTENT(IN):: ID
	REAL*8 THETA(JMAX+1),RSUM(JMAX+1)
	REAL*8 MYU1
	COMMON /HEIKIN/MYU1
	INTEGER J1,M1,IKK,JKK
	REAL*8 RES, RES1, RES2, MT, RESH,RESDH
	INTEGER RESI, RESJ
	REAL*8 XKK(MMAX,JMAX)
	LOGICAL CHECK
	REAL*8 RESZ,RESX(MMAX,JMAX)
	INTEGER I,J,K,M,L

	M1 = MOD(ID,MMAX); IF(M1==0) M1=MMAX
	J1 = (ID-M1)/MMAX + 1

	BETA(J1,M1)=V
	H = 0.0; DH = 0.0
	DO 100 K = 1, NDATA
		THETA=0.0
		IKK=IK(K); JKK=JK(K)
		IF(IKK > J1 .OR. JKK<J1) CYCLE
		DO J=1,JMAX;DO M=1,MMAX
			XKK(M,J)=XK(M,K)
		END DO; END DO

		DO I=1,D1	! パラメータ空間から除去
			M = MOD(DEL(I)+(MMAX-1),MMAX) + 1; J = (DEL(I)-M)/MMAX+1
			XKK(M,J)=0.0
		END DO

!----------------------------------------
*	特別設定（計算省力化）
!----------------------------------------
		IF (IK(K)==RESI .AND. JK(K)==RESJ .AND. ZK(K)==RESZ) THEN
			CHECK=.TRUE.
			DO J=1,JMAX; DO M=1,MMAX
				IF(XKK(M,J) /= RESX(M,J)) THEN
					CHECK=.FALSE.; GOTO 41
				END IF
			END DO; END DO

 41			CONTINUE
			IF (CHECK==.TRUE.) THEN
				H = H + RESH
				DH = DH + RESDH
				CYCLE
			END IF
		END IF
!----------------------------------------
		RESI=IK(K); RESJ=JK(K); RESZ=ZK(K); RESX=XKK

		DO J=IKK,JKK
			IF(J==JMAX+1) EXIT
			DO M = 1, MMAX
				THETA(J)=THETA(J) + BETA(J,M)*XKK(M,J1)
			END DO
			THETA(J)=DEXP(THETA(J))
		END DO
		! 最小値の探索
		MT = THETA(IKK)
		DO J=IKK+1,JKK
			IF(THETA(J) < MT) MT=THETA(J)
		END DO
		DO L=IKK,JKK	!最小のΘで0になるように調整
			RSUM(L) = DEXP((MT-THETA(L))*ZK(K)) / PROD1(IKK,JKK,L,THETA,JMAX)
		END DO
		
		RES1 = 0.0; RES2=0.0
		DO L= IKK,JKK
			IF(L==J1) CYCLE
			RES = 1.0/(THETA(L)-THETA(J1))
			RES1 = RES1 + RES*RSUM(L);	RES2 = RES2 + RES
		END DO
		RES1 = RES1 + (RES2-ZK(K))*RSUM(J1)

		RES = 0.0
		DO L= IKK,JKK; RES = RES + RSUM(L);	END DO

		IF(RES == 0.0) THEN
			RESH = -DEXP(300.0D0) - MT*ZK(K) 
		ELSE IF(RES<0.0) THEN
			PRINT *, K, "値がマイナス．WHY ? "; STOP
		ELSE
			RESH = DLOG(RES) - MT*ZK(K)
		END IF
		RESDH = RES1/RES*THETA(J1)*XKK(M1,J1)

		IF(J1/=JKK) THEN
			RESH = RESH + BETA(J1,M1)*XKK(M1,J1)
			RESDH = RESDH + XKK(M1,J1)
		END IF
		H = H + RESH
		DH = DH + RESDH
 100	CONTINUE

	IF(PATTERN==3) RETURN

	IF(PATTERN==1 .OR. M1==1) THEN
		IF(PATTERN==2) MYU1 = MYU0(M1,J1)
		RES = BETA(J1,M1) - MYU1
		H = H - ROH(M1,M1,J1) * RES**2/2.0D0
		DH = DH - ROH(M1,M1,J1) * RES
	END IF

	RETURN; END SUBROUTINE
!------------------------------------------------------
*	\prod_{e=a,\ne c}^{b} \frac{1}{\bar{theta}_{e,c}}
!------------------------------------------------------
	REAL*8 FUNCTION PROD1(A,B,C,THETA,JMAX)
	IMPLICIT NONE
	INTEGER,INTENT(IN):: A,B,C,JMAX
	REAL*8,INTENT(IN):: THETA(JMAX+1)
	INTEGER E

	PROD1 = 1.0
	DO E = A, B
		IF(E /= C) PROD1 = (THETA(E)-THETA(C)) * PROD1
	END DO
	RETURN;	END
!------------------------------------------------------
	SUBROUTINE INITIALIZE(ID)
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
*	インプットファイルからの入力
***********************************************************************
	SUBROUTINE INPUT()
	USE VARIABLE; USE FILE; IMPLICIT NONE
	REAL*8 MAXX(MMAX), SCALE
	COMMON /ARR/MAXX,SCALE
	CHARACTER CH
	INTEGER J,K,M1,M2

	SCALE=1.0; MAXX=0.0; MAXX(1)=1.0
	XK=1.0
	OPEN(1000,file=INPUT_F1)
		READ(1000,*) CH
		DO K=1,NDATA
			READ(1000,*) IK(K),JK(K),ZK(K),(XK(M1,K),M1=2,MMAX)
			DO M1=2,MMAX; IF(XK(M1,K) > MAXX(M1)) MAXX(M1)=XK(M1,K); END DO
		END DO
	CLOSE(1000)
!		スケール調整
	DO K=1,NDATA
		DO M1=2,MMAX; XK(M1,K)=XK(M1,K)/MAXX(M1); END DO
	END DO
	XK=XK*SCALE

	IF(PATTERN==3) RETURN	!無情報の場合

	OPEN(1001,file=INPUT_F2)
		DO J=1,JMAX
			READ(1001,*) CH
			READ(1001,*) CH,(MYU0(M1,J),M1=1,MMAX)
			DO M1=1,MMAX
				READ(1001,*) CH,(ROH(M1,M2,J),M2=1,MMAX)
			END DO
		END DO
	CLOSE(1001)

	PRINT *, "事前情報[平均]"
	DO J=1,JMAX
		PRINT 30, J,(MYU0(M1,J),M1=1,MMAX)
	END DO
30	FORMAT(1X,I2,6F10.4)

!		逆行列へ変換・スケール調整
	CALL REVERSE(1)

	RETURN
	END SUBROUTINE INPUT
***************************************************
*	アウトプット
***************************************************
	SUBROUTINE OUTPUT(TANE,SAMP,NPOINT)
	USE VARIABLE; USE FILE; IMPLICIT NONE
	INTEGER,INTENT(IN):: TANE,NPOINT
	REAL*8,INTENT(IN):: SAMP(NPOINT,JMAX*MMAX)	!サンプリング集合
	REAL*8 SAMP1(NPOINT)	!集約化
	REAL*8 MYU(MMAX,JMAX),SS(MMAX,MMAX,JMAX),SM	!平均，分散
	REAL*8 UBAR(JMAX*MMAX),OBAR(JMAX*MMAX)	!信頼区間
	REAL Z(JMAX*MMAX); COMMON /GEWEKE/Z	!GEWEKE仮説検定統計量
	INTEGER I,J,K,M1,M2
	
!--------------------------
*	サンプリングの出力
!--------------------------
	OPEN(2000,file=OUTPUT_F1)
		WRITE(2000,*) "乱数の種",",",TANE
		WRITE(2000,*) "サンプル数",",",NPOINT
		WRITE(2000,*) "状態1",",,,","状態2",",,,","状態3",",,,","状態4"
		WRITE(2000,*) ("B1",",","B2",",","B3",",",J=1,4)
		DO K=1,NPOINT
			WRITE(2000,20) (SAMP(K,J),J=1,JMAX*MMAX)
		END DO
	CLOSE(2000)
 20	FORMAT(1X,30(E15.8,","))

!--------------------------
*	統計量の出力
!--------------------------
!	標本平均
	MYU=0.0
	DO J=1,JMAX; DO M1=1,MMAX
		DO K=1,NPOINT
			MYU(M1,J)=MYU(M1,J)+SAMP(K,(J-1)*MMAX+M1)
		END DO
	END DO;	END DO	
	MYU=MYU/NPOINT

!	標本分散
	DO J=1, JMAX; DO M2=1,MMAX; DO M1=1,MMAX
		I = (J-1)*MMAX
		DO K=1,NPOINT
			SS(M1,M2,J)=SS(M1,M2,J)
     &			+(SAMP(K,I+M1)-MYU(M1,J))*(SAMP(K,I+M2)-MYU(M2,J))
		END DO
	END DO; END DO; END DO
	SS=SS/NPOINT

!	信頼区間
	DO J=1,JMAX*MMAX
		DO K=1,NPOINT
			SAMP1(K)=SAMP(K,J)
		END DO
		CALL BUBBLE(SAMP1,NPOINT)
		UBAR(J)=SAMP1(NPOINT*5/100)
		OBAR(J)=SAMP1(NPOINT*95/100)
	END DO
 90	FORMAT(1X, F10.2, A12)

!	ファイル出力	
	OPEN(3000,file=OUTPUT_F2)
		SELECT CASE(PATTERN)
			CASE(1); WRITE(3000,*) "全事前情報"
			CASE(2); WRITE(3000,*) "部分事前情報"
			CASE(3); WRITE(3000,*) "無情報"
		END SELECT
	
		IF(PATTERN/=3) THEN
			WRITE(3000,*) "初期情報"
			WRITE(3000,*) 
     &		",", "状態1",",,,","状態2",",,,","状態3",",,,","状態4",",,,"
			WRITE(3000,10) "初期平均",((MYU0(M1,J),M1=1,MMAX),J=1,JMAX)
			DO M1=1,MMAX
				WRITE(3000,10) "初期分散",((ROH(M1,M2,J),M2=1,MMAX),J=1,JMAX)
			END DO
		END IF
				
		WRITE(3000,*) "事後情報"
		WRITE(3000,*) 
     &		",", "状態1",",,,","状態2",",,,","状態3",",,,","状態4",",,,"
		WRITE(3000,10) "標本平均",((MYU(M1,J),M1=1,MMAX),J=1,JMAX)
		DO M1=1,MMAX
			WRITE(3000,10) "標本分散",((SS(M1,M2,J),M2=1,MMAX),J=1,JMAX)
		END DO

		WRITE(3000,10) "信頼区間",(UBAR(J),J=1,JMAX*MMAX)
		WRITE(3000,10) "信頼区間",(OBAR(J),J=1,JMAX*MMAX)
		WRITE(3000,10) "収束判定",(Z(J),J=1,JMAX*MMAX)
 10	FORMAT(1X,A8,",",20(E15.8,","))
	CLOSE(3000)

	RETURN;	END SUBROUTINE OUTPUT
**********************************************
*	スケール調整
**********************************************
	SUBROUTINE REVERSE(IP)
	USE VARIABLE; IMPLICIT NONE
	INTEGER,INTENT(IN):: IP
	REAL*8 A(MMAX,MMAX)
	REAL*8 MAXX(MMAX)	!説明変数の最大値
	REAL*8 SCALE
	COMMON /ARR/MAXX,SCALE
	INTEGER J,M1,M2
	
	IF(IP==1) THEN
		DO J=1,JMAX
			DO M1=2,MMAX; MYU0(M1,J)=MYU0(M1,J)*MAXX(M1); END DO
			DO M1=1,MMAX; DO M2=1,MMAX
				A(M2,M1) = ROH(M2,M1,J)*MAXX(M2)*MAXX(M1)/SCALE**2
			END DO; END DO
			CALL MATRIX_INVERSE(A,MMAX)
			DO M1=1,MMAX; DO M2=1,MMAX
				ROH(M2,M1,J) = A(M2,M1)
			END DO; END DO
		END DO
		MYU0 = MYU0/SCALE
	ELSE	! 逆調整
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
		MYU0 = MYU0*SCALE
		ROH = ROH * SCALE**2
	END IF
	RETURN; END SUBROUTINE
**************************************************************
*	ギブスサンプリングを用いたMCMC法
**************************************************************
	SUBROUTINE MCMC_G(SAMP,NPOINT)
	USE VARIABLE; IMPLICIT NONE
	INTEGER,INTENT(IN):: NPOINT ! 希望サンプリング数
	REAL*8,INTENT(OUT):: SAMP(NPOINT,JMAX*MMAX)	!サンプリング集合
	REAL*8 SAMP1(NPOINT)		!SUBサンプリング集合
	INTEGER IBURN	
*	ARSを支配するパラメータ
	INTEGER,PARAMETER:: KMAX=2	!初期基準座標数
	INTEGER NT	!Tk：x座標の最大数
	REAL*8 SXK(JMAX*MMAX,KMAX)	!スタートポイント
!-----------------------------------------
	LOGICAL OB(JMAX,MMAX),UB(JMAX,MMAX)	!領域Dの上限下限の存在
	REAL*8 OX(JMAX,MMAX),UX(JMAX,MMAX)	!領域Dの上限下限
	REAL,EXTERNAL:: GEWEKEZ
	REAL Z(JMAX*MMAX); COMMON /GEWEKE/Z	!GEWEKE仮説検定統計量
	REAL*8 MAXX(MMAX), SCALE
	COMMON /ARR/MAXX,SCALE
	LOGICAL KENTEI
	CHARACTER(LEN=30):: hour
!-----------------------------------------
	REAL*8 X	!戻り値
!-----------------------------------------
	INTEGER IPOINT
	INTEGER I,J,M,T,ID, IREV,d
	LOGICAL STAY

*	ARS条件の設定
	OB=.FALSE.; UB=.FALSE.
	CALL FIRST(IBURN,KMAX,NT,SXK,JMAX,MMAX)
!*********************
*	Burn-in
!*********************
		PRINT *, "Burn-in期間"
		DO T=1,IBURN
		I=1
			IF(MOD(T,1000)==0.0) PRINT *, T,"回目"
			DO J=1,JMAX; DO M=1,MMAX
				ID = (J-1)*MMAX+M
				
				IF(DEL(I)==ID) THEN
					I=I+1; CYCLE
				END IF

				IF(PATTERN/=3) CALL INITIALIZE(ID)
				STAY=.TRUE.
				CALL GIBBS(X,UB(J,M),OB(J,M),UX(J,M),OX(J,M))
				IF(STAY==.FALSE.) STOP
			END DO;	END DO
		END DO
!*********************
*	Sampling
!*********************
		SAMP=0.0; IREV=0; IPOINT=1
 100	CONTINUE
	CALL TIME(hour); PRINT *, hour
		PRINT *, "Sampling期間"
		DO T=IPOINT,NPOINT
		I=1
			IF(MOD(T,1000)==0.0) PRINT *, T,"回目"
			DO J=1,JMAX; DO M=1,MMAX
				ID = (J-1)*MMAX+M

				IF(DEL(I)==ID) THEN
					I=I+1; CYCLE
				END IF

				IF(PATTERN/=3) CALL INITIALIZE(ID)
				
				CALL GIBBS(X,UB(J,M),OB(J,M),UX(J,M),OX(J,M))
				IF(STAY==.FALSE.) STOP
				BETA(J,M)=X; SAMP(T,ID)=X/MAXX(M)*SCALE
			END DO;	END DO
		END DO

!*********************
*	仮説検定
!*********************
	KENTEI=.TRUE.
	I=1
		DO J=1,JMAX; DO M=1,MMAX
			ID = (J-1)*MMAX+M
			DO T=1,NPOINT
				SAMP1(T)=SAMP(T,ID)
			END DO
	
			IF(DEL(I)==ID) THEN
				I=I+1; CYCLE
			END IF

			Z(ID)= GEWEKEZ(SAMP1,NPOINT)
			IF(Z(ID) < 1.96) THEN	!５％有意水準
				PRINT 90, ID,Z(ID), "収束"
			ELSE
				PRINT 90, ID,Z(ID), "収束せず"
				KENTEI=.FALSE.
			END IF
		END DO; END DO
 90	FORMAT(1X, I3,1X,F10.2, A12)

	IF(IREV>=LIMIT) RETURN
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

!---------------------------------
*	内部手続き
!---------------------------------
	CONTAINS
	SUBROUTINE GIBBS(XJ,UBJ,OBJ,UXJ,OXJ)
		LOGICAL,INTENT(IN):: UBJ,OBJ
		REAL*8,INTENT(IN):: UXJ,OXJ
		REAL*8,INTENT(OUT):: XJ
		INTEGER IFAULT,I,K
		REAL*8 XX(KMAX)	!初期座標
		REAL*8 RES

		DO K=1,KMAX; XX(K)=SXK(ID,K); END DO
		DO I=1,1000
			IFAULT=0
			CALL REJECT_SAMPLE(XJ,NT,KMAX,UBJ,OBJ,UXJ,OXJ,XX,ID,IFAULT)
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
				PRINT *, "WHY?"; STAY=.FALSE.
			END IF
		END DO
		SXK(ID,1)=XX(1); SXK(ID,KMAX)=XX(KMAX)
	END SUBROUTINE GIBBS
!---------------------------------------------------------
	END SUBROUTINE MCMC_G
*******************************************************
*	適応的棄却サンプリング
*******************************************************
	SUBROUTINE REJECT_SAMPLE(XPOINT,NT,NK,UB,OB,UX,OX,XX,NO,IFAULT)
	IMPLICIT NONE
	INTEGER,INTENT(IN):: NO	!どの条件付き確率分布を発生させるか
	INTEGER,INTENT(IN):: NK
	INTEGER IFAULT	!エラー判定
	REAL*8,EXTERNAL:: RANSU1
	REAL*8,PARAMETER:: EPS10=1.0D-16	!ε
	REAL*8 XPOINT, XX(NK)
	INTEGER,INTENT(IN):: NT	!x座標の最大数
	LOGICAL,INTENT(IN):: OB, UB	!領域Dの上限下限の存在
	REAL*8,INTENT(IN):: UX,OX	!下限，上限
	REAL*8 HK(NT),DHK(NT),ZK(0:NT),HZK(0:NT),XKK(NT)
	INTEGER,PARAMETER:: IEMAX=700	!EXPの最大引数
	INTEGER,PARAMETER:: NIT=1000	!繰返し数
	REAL*8 RX	!発生座標
	REAL*8 W	!一様乱数
	REAL*8 LK,UK,HX,RH,RDH,RZ1,RZ2,RZH1,RZH2
	INTEGER RJ	!挿入場所
	INTEGER DR	!TEST結果
	REAL*8 HZMAX	!Hのスケール調整
	REAL*8 RESERVE
	INTEGER K
	INTEGER I,IT

!---------------------------
!	Initialization step
!---------------------------
	K=NK
	DO I=1,K
		XKK(I)=XX(I)	!初期基準座標 Tk
		XPOINT=XKK(I)
		CALL MAKEH(XPOINT,HK(I),DHK(I),NO)
	END DO

	IF(UB==.FALSE. .AND. DHK(1) < 0.0) THEN
		PRINT *, "ID=",NO,"傾き",DHK(1),"左範囲をもっと広く"
		IFAULT=1; RETURN
	END IF

	IF(OB==.FALSE. .AND. DHK(K) > 0.0) THEN
		PRINT *, "ID=",NO,"傾き",DHK(K),"右範囲をもっと広く"
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
!	Skから乱数発生
		CALL SKSAMPLE(RX,RJ,K,IFAULT)
		IF(IFAULT==6) RETURN
!	一様乱数を発生
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
				PRINT *, "ID=",NO, "Tkがいっぱいです"; IFAULT=3; RETURN
			END IF
!	Tk,Zjを昇順に並べ替え
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

	PRINT *, "必要なサンプリング数が集まらず"
	RETURN
!---------------------------------------------------------
*	内部手続き
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
*	Skからの乱数発生
!	遺憾ながらGilksのプログラムを大いに参照
**************************************************
	SUBROUTINE SKSAMPLE(X,JJ,KMAX,IFAULT)
	REAL*8,EXTERNAL:: RANSU2
!	求めたい値
	INTEGER,INTENT(OUT):: JJ	!J:挿入場所
	REAL*8,INTENT(OUT):: X	!X:サンプル座標
	INTEGER,INTENT(IN):: KMAX
	INTEGER,INTENT(OUT):: IFAULT
	REAL*8 SSUM(KMAX)	!Uk下の合計面積
	REAL*8 LSUM
	REAL*8 EH,SIGN	!値保管用
	LOGICAL HORIZ
	REAL*8 R	!乱数
	INTEGER J

!	直線Uk下の各面積を求める
	CALL UAREA(KMAX,SSUM)
	IF(SSUM(KMAX)<=0.0) THEN
		IFAULT=6; PRINT *, IFAULT,"面積が0以下"; RETURN
	END IF
	LSUM=DLOG(SSUM(KMAX))
!	分布確率に変換
	DO J=1,KMAX-1
		R=SSUM(J)/SSUM(KMAX)
		SSUM(J)=R
	END DO
	SSUM(KMAX)=1.0

!	（0,1）の一様乱数を発生させ，エリアJを決定
	R=RANSU2()
	IF(R==0.0) THEN
		PRINT *, "乱数2が0に"; STOP
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
	SUBROUTINE UAREA(KMAX,SSUM)	! 直線Uk下の各面積を求める
	INTEGER,INTENT(IN):: KMAX
	REAL*8,INTENT(OUT):: SSUM(KMAX)	!Uk下の合計面積
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
*******************************************************
*     定常分布への収束判定（GEWEKEの仮説検定）
*******************************************************
	REAL FUNCTION GEWEKEZ(SAMP,NPOINT)
	IMPLICIT NONE
	REAL*8,EXTERNAL:: TVORA
	INTEGER,INTENT(IN):: NPOINT ! 標本数
	REAL*8,INTENT(IN):: SAMP(NPOINT)	!標本集合
	INTEGER NP1, NP2
	REAL*8 SAMP1(NPOINT/10), SAMP2(NPOINT/50)
	REAL*8 MYU1,MYU2,S1,S2
	REAL*8 Z	! 検定統計量
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
*	標本分散の計算（時系列）
! 参考文献
! A Simple, Positive Semi-Definite, Heteroskedasticity
!	and Autocorrelation Consistent Covariance Matrix
*****************************************
	REAL*8 FUNCTION TVORA(SAMPX,MYU,NN)
	IMPLICIT NONE
	INTEGER,INTENT(IN):: NN
	REAL*8,INTENT(IN):: SAMPX(NN),MYU
	INTEGER I,J,M

		M = 10	! バンド幅
		TVORA= OMEGA(0)
		DO J=1, M
			TVORA = TVORA + (1.0 - DBLE(J)/(M+1.0))*2.0*OMEGA(J)
		END DO
	CONTAINS
!--------------------
*	内部手続き
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
***************************************************
*	疑似一様乱数の発生（乗算合同式法）
***************************************************
	REAL*8 FUNCTION RANSU1()
	INTEGER IR1
	COMMON /IRAND1/IR1
		IR1=IAND(48828125*IR1, 2147483647)
		RANSU1=DBLE(IR1)/DBLE(2147483647)
	END FUNCTION
***************************************************
*	疑似一様乱数の発生（マルチ乗算合同式法：AS183）
!	周期が長い
***************************************************
	REAL*8 FUNCTION RANSU2()
!	乱数の種はいずれも1から30000の範囲で与えなければいけない
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
*	バブルソート
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
*	逆行列 (完全ピボット法)
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
!	(N-K+1)次小行列の中から絶対値が最大のものを探索
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
!	絶対値が最大のものを左上に持ってくる
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
!	Ｋ行に対する操作
		AW=1.0/A(K,K);	A(K,K)=1.0
		DO J=1,N; A(K,J)=A(K,J)*AW;	END DO
!	Ｋ行以外への操作
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