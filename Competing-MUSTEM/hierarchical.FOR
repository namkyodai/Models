************************************************************************
*	階層型指数劣化ハザードモデル(面が吸収状態バージョン)
************************************************************************
	MODULE VARIABLE
		INTEGER,PARAMETER:: IMAX=3,JMAX=3,MMAX=3 !劣化ランク、タイプ、特性数
		INTEGER,PARAMETER:: EMAX=9 !実際のレーティング
		!------------------------------------------
		INTEGER,PARAMETER:: TOTALDATA=4721	!使用データ数
	!------------------------------------------
		INTEGER,PARAMETER:: TIMEMAX=1	!繰り返し回数
		INTEGER,PARAMETER:: NMAX=(1+(IMAX-1)*JMAX)*JMAX		!ρの数
		INTEGER,PARAMETER:: dimN=NMAX*MMAX	!連立方程式の数
		INTEGER,PARAMETER:: XMAX=9 !逆行列Ｂの最大行数
		INTEGER,PARAMETER:: YMAX=9 !逆行列Ｂの最大列数

	!---推定パラメータ数を変更したい場合-------------
		INTEGER,PARAMETER:: D1=49	!削除パラメータ数
		INTEGER,PARAMETER:: D2=1	!削除パラメータ数
		INTEGER,PARAMETER:: redimN=dimN-D1	!修正後パラメータ数
		INTEGER DEL(D1)
		INTEGER DEL2(D2)
		
	!削除するパラメータ番号((J-1)*MMAX+M)
	DATA DEL/63,62,60,59,57,56,54,53,51,50,48,47,45,44,42,41,39,38,34,
     &36,35,33,32,30,29,27,26,24,23,21,20,18,17,15,14,12,11,9,8,6,5,3,2
     &,61,58,55/
	DATA DEL2/18/
!	DATA DEL2/66,62,22/ 	!尤度比検定用削除パラメータ
	!-----------------------------------------
		!データのレコード数	 MAX10231
!		CHARACTER(LEN=19),PARAMETER:: FILE1="input\inputdata.csv"
		CHARACTER(LEN=30),PARAMETER:: FILE1="input\inputdata(h).csv"
!		CHARACTER(LEN=20),PARAMETER:: FILE1="input\inputdata3.csv"
		CHARACTER(LEN=13),PARAMETER:: FILE1000="output\df.csv"
		CHARACTER(LEN=17),PARAMETER:: FILE1001="output\beta-h.csv"
		CHARACTER(LEN=17),PARAMETER:: FILE1002="output\hazard.csv"
		CHARACTER(LEN=15),PARAMETER:: FILE1003="output\prob.csv"
		CHARACTER(LEN=16),PARAMETER:: FILE1004="output\delta.csv"
		CHARACTER(LEN=17),PARAMETER:: FILE1005="output\likely.csv"
		CHARACTER(LEN=18),PARAMETER:: FILE1006="output\t-value.csv"
		REAL*8:: THETA(0:IMAX,0:JMAX) 
		!ハザード関数Θij ランク最大は０として存在
		REAL*8:: ETA(0:IMAX,0:JMAX,0:IMAX,0:JMAX)	
		!Θij-Θlrランク最大は０として存在
		REAL*8:: RHO(0:IMAX,0:JMAX,1:JMAX) !推移強度ρijr
		INTEGER:: MAXTHETA
		INTEGER:: MAXRHO
		CHARACTER CH
		REAL Z1
		DATA Z1/1.0/
	END MODULE
***********************************************************************
*	メインプログラム
***********************************************************************
	PROGRAM ESTIMATION
	USE VARIABLE; IMPLICIT NONE
	REAL*8,EXTERNAL:: PROB,OBJ
	REAL*8:: TOPTBETA(MMAX,NMAX),TBETA(MMAX,NMAX),TBETAS(MMAX,NMAX)
	REAL*8:: OPTBETA(1:MMAX,0:IMAX-1,0:JMAX,1:JMAX)
	REAL*8:: BETA(1:MMAX,0:IMAX-1,0:JMAX,1:JMAX)
	REAL*8:: PI(IMAX*JMAX+1,IMAX*JMAX+1)
	REAL*8:: BETAS(1:MMAX,0:IMAX-1,0:JMAX,1:JMAX)
	REAL*8:: BETAA(1:MMAX,0:IMAX-1,0:JMAX,1:JMAX)
	REAL*8:: DELTA(1:MMAX,0:IMAX-1,0:JMAX,1:JMAX)
	REAL*8:: TBETAA(MMAX,NMAX),TDELTA(MMAX,NMAX),ALPHA(MMAX,NMAX)
	REAL*8:: G(dimN,dimN),A(redimN,redimN),OPTOBJ,TESTOBJ,TOKEI,QQ
	REAL*8,PARAMETER:: EPS1=1.0D-08
	REAL*8:: DF2(dimN,dimN)
	REAL*8 DF2_0(JMAX*MMAX,JMAX*MMAX),LR(MMAX,NMAX)
	REAL*8 DF2_A(JMAX*MMAX,(IMAX-1)*JMAX*JMAX*MMAX)
	REAL*8 DF2_B((IMAX-1)*JMAX*JMAX*MMAX,JMAX*MMAX)
	REAL*8 DF2_C((IMAX-1)*JMAX*JMAX*MMAX,(IMAX-1)*JMAX*JMAX*MMAX)
	REAL TV(MMAX)
	REAL Z
	INTEGER I,J,R,M,TIME,EVAL,L,S,P,Q
	INTEGER W,X,Y,V,TIMER
*----------------------------
*	最尤推定法
*----------------------------
!-----初期パラメータβ------------------------------

	OPEN(1001,file=FILE1001)
	READ(1001,*) CH
	DO S=1,NMAX
	READ(1001,*) (TBETA(M,S),M=1,MMAX)
	print*, (TBETA(M,S),M=1,MMAX)
	END DO
	CLOSE(1001)
	CALL TRANSITION_2(TBETA,BETA)
	DELTA=0.09
!----------------------------------------------------

!-----擬似パターン法-----------------
	CALL PATTERN(BETA,DELTA)
!------------------------------------
	OPTBETA=BETA
*-------------------
*	検定
*-------------------
	CALL TRANSITION_1(OPTBETA,TOPTBETA)
	!ハザード関数・推移確率行列・未知パラメータβ
	CALL PRINT1(OPTBETA,Z)
****************************************************
*	Ｔ値算出用
****************************************************
!	OPEN(1004,file=FILE1004)
!	DO S=1,NMAX
!	READ(1004,*) (TDELTA(M,S),M=1,MMAX)
!	END DO
!	CLOSE(1004)
!	CALL TRANSITION_2(TDELTA,DELTA)
!		CALL HESSIAN_0(DELTA,BETA,DF2_0)
!		CALL HESSIAN_A(DELTA,BETA,DF2_A)
!		CALL HESSIAN_B(DELTA,BETA,DF2_B)
!		CALL HESSIAN_C(DELTA,BETA,DF2_C)
!	DO P=1,JMAX*MMAX; DO Q=1,JMAX*MMAX
!	DF2(P,Q)=DF2_0(P,Q)
!	END DO; END DO
	
!	DO P=1,JMAX*MMAX; DO Q=1,(IMAX-1)*JMAX*JMAX*MMAX
!	DF2(P,Q+JMAX*MMAX)=DF2_A(P,Q)
!	END DO; END DO
!	DO P=1,(IMAX-1)*JMAX*JMAX*MMAX; DO Q=1,JMAX*MMAX
!	DF2(P+JMAX*MMAX,Q)=DF2_B(P,Q)
!	END DO; END DO
	
!	DO P=1,(IMAX-1)*JMAX*JMAX*MMAX; DO Q=1,(IMAX-1)*JMAX*JMAX*MMAX
!	DF2(P+JMAX*MMAX,Q+JMAX*MMAX)=DF2_C(P,Q)
!	END DO; END DO
!	G=DF2
!	CALL DEL_MATRIX(G,A)
!	A=-A; CALL MATRIX_INVERSE(A,redimN)	!フィッシャー行列の逆行列
!	OPEN(1006,file=FILE1006)
!	PRINT *, "t値********"
!	WRITE(1006,*) "t値********"
!	L=0
!	DO J=1,NMAX
!		DO 210 M=1,MMAX
!			DO I=1,D1
!				IF ((J-1)*MMAX+M==DEL(I))  THEN
!					TV(M)=0.0; GOTO 210 !削除されてる場合t値は0
!				END IF
!			END DO
!			L=L+1; TV(M)=TOPTBETA(M,J)/DSQRT(A(L,L))
!210		CONTINUE
!		PRINT *, (TV(M),M=1,MMAX)
!		WRITE(1006,*) (TV(M),M=1,MMAX)
!	END DO
!	CLOSE(1006)
 11	FORMAT(1X,3(F12.8,","))
20	FORMAT(1X,2(F12.8,","),F12.8)
***********************************************
	PI=0.0
****************************************************
*	推移行列
****************************************************
	PI(1,1)=PROB(0,0,0,0,Z1)
	DO L=1,IMAX; DO R=1,JMAX
	I=0; J=0
	PI(1,3*(L-1)+1+R)=PROB(I,J,L,R,Z1)
	END DO; END DO
	DO I=1,IMAX; DO J=1,JMAX; DO L=I,IMAX; DO R=1,JMAX
	PI(3*(I-1)+1+J,3*(L-1)+1+R)=PROB(I,J,L,R,Z1)
	END DO; END DO; END DO; END DO
	!π_ij
	OPEN(1003,file=FILE1003)
	PRINT *, "推移確率"
	WRITE(1003,*) "推移確率"
	DO I=1,IMAX*JMAX+1
		PRINT 10, (PI(I,J),J=1,IMAX*JMAX+1)
		WRITE(1003,10) (PI(I,J),J=1,IMAX*JMAX+1)
	END DO
 10	FORMAT(1X,10(F12.8,","),F12.8)
	CLOSE(1003)
	STOP
	END PROGRAM

*******************************************************
*	擬似パターン法による直接探索法
*******************************************************
	SUBROUTINE PATTERN(BETA,DELTA)
	USE VARIABLE; IMPLICIT NONE
	REAL*8, INTENT(IN):: BETA(1:MMAX,0:IMAX-1,0:JMAX,1:JMAX)
	REAL*8, INTENT(IN):: DELTA(1:MMAX,0:IMAX-1,0:JMAX,1:JMAX)
	REAL*8, EXTERNAL:: OBJ
	REAL*8:: BETAA(1:MMAX,0:IMAX-1,0:JMAX,1:JMAX)
	REAL*8:: BETAB(1:MMAX,0:IMAX-1,0:JMAX,1:JMAX)
	REAL*8:: BETAD(1:MMAX,0:IMAX-1,0:JMAX,1:JMAX)
	REAL*8:: BETAE(1:MMAX,0:IMAX-1,0:JMAX,1:JMAX)
	REAL*8:: POS(1:MMAX,1:NMAX),W,H,LR(MMAX,NMAX)
	REAL*8:: TBETA(1:MMAX,1:NMAX)
	REAL*8:: TBETAA(1:MMAX,1:NMAX)
	REAL*8:: TBETAB(1:MMAX,1:NMAX)
	REAL*8:: TBETAD(1:MMAX,1:NMAX)
	REAL*8:: TDELTA(1:MMAX,1:NMAX)
	REAL*8:: OBJE(D2)
	INTEGER:: S,TIME,M,P,Q,S1,M1,D,D22,D222,T
	REAL*8,PARAMETER:: EPS=1.0D-15
	OBJE=0.0
	CALL TRANSITION_1(DELTA,TDELTA)
	CALL TRANSITION_1(BETA,TBETA)
	
	DO 80 D22=1,D2
	DO 88 T=1,1
	DO D=1,D1; DO S=1,NMAX; DO M=1,MMAX
	IF(MMAX*(S-1)+M==DEL(D))THEN
	TBETA(M,S)=0.0
	END IF
	IF(MMAX*(S-1)+M==DEL2(D22)) THEN
	TBETA(M,S)=0.0
	END IF
	END DO;END DO; END DO
	TBETAD=TBETA
	CALL TRANSITION_2(TBETAD,BETAD)
	TDELTA=0.09
	DO 22 TIME=1,1
	PRINT*,TIME
	CALL TRANSITION_1(BETAD,TBETAD)
	DO M=1,MMAX; DO S=1,NMAX
	IF(MMAX*(S-1)+M==DEL2(D22)) GO TO 55
	DO D=1,D1
	IF(MMAX*(S-1)+M==DEL(D)) GO TO 55
	END DO
			TBETAA=TBETAD
			TBETAB=TBETAD
			POS=TBETAD
			TBETAA(M,S)=TBETAA(M,S)+TDELTA(M,S)
			TBETAB(M,S)=TBETAB(M,S)-TDELTA(M,S)
		CALL TRANSITION_2(TBETAA,BETAA)
		CALL TRANSITION_2(TBETAB,BETAB)
		CALL TRANSITION_2(TBETAD,BETAD)
	H=OBJ(BETAD)
			IF(OBJ(BETAA)-H > 0.0) THEN
				POS(M,S)=TBETAA(M,S)
			ELSE IF(OBJ(BETAB)-H > 0.0) THEN
				POS(M,S)=TBETAB(M,S)
			ELSE
			END IF
			TBETAD=POS
55	END DO
	END DO
		DO 33 S=1,NMAX
		DO M=1,MMAX
			IF(DABS(TBETAD(M,S)-TBETA(M,S)) < EPS) THEN
				TDELTA(M,S)=0.7*TDELTA(M,S)
			END IF
!			IF(TDELTA(M,S) < EPS) GO TO 99
		END DO
33		END DO
			TBETA=TBETAD

	CALL TRANSITION_2(TBETAD,BETAD)
	OPEN(1001,file=FILE1001)
	WRITE(1001,*) "未知パラメータβ***************"
	DO S=1,NMAX
		WRITE(1001,20) (TBETA(M,S),M=1,MMAX)
	END DO 
	CLOSE(1001)

	OPEN(1004,file=FILE1004)
	WRITE(1001,*) "ステップ幅δ***************"
	DO S=1,NMAX
		WRITE(1004,20) (TDELTA(M,S),M=1,MMAX)
	END DO 
	CLOSE(1004)
 		PRINT *, TBETA
20	FORMAT(1X,2(F12.8,","),F12.8)
	
	CALL TRANSITION_2(TBETA,BETAD)
	W=OBJ(BETAD)
	OBJE(D22)=W
	PRINT*,"対数尤度",W
	OPEN(1005,file=FILE1005)
	DO D222=1,D2
	WRITE(1005,*) OBJE(D222)
	END DO
	CLOSE(1005)
22	END DO
88	END DO
80	END DO
99	CALL TRANSITION_2(TBETA,BETA)
	CALL TRANSITION_2(TDELTA,DELTA)

	END SUBROUTINE PATTERN

*******************************************
*	ρijr(i=IK,LK),θの値を求める
*******************************************
	SUBROUTINE THETAVALUE(IK,JK,LK,RK,XK,BETA)
	USE VARIABLE; IMPLICIT NONE
	REAL*8,INTENT(IN):: BETA(1:MMAX,0:IMAX-1,0:JMAX,1:JMAX)
	REAL*8,INTENT(IN):: XK(1:MMAX,0:IMAX-1,0:JMAX,1:JMAX)
	INTEGER,INTENT(IN):: IK,JK,LK,RK
	INTEGER I,J,L,M,R
	REAL*8 RESERVE,RESERVE2

	!ρの計算
	RHO=0.0
	MAXRHO=IK;RESERVE=0.0
	DO 100 I=0,IMAX
	DO J=0,JMAX                                 
	DO R=1,JMAX									
		IF(I==IMAX) CYCLE !ランク3は推移強度無し
		IF(J==JMAX) CYCLE !面はすべて吸収状態
		DO M=1,MMAX
			RHO(I,J,R)=RHO(I,J,R)+XK(M,I,J,R)*BETA(M,I,J,R)
		END DO
		RHO(I,J,R)=EXP(RHO(I,J,R))
		!オーバーフローを防ぐため最大のρを調べておく
		IF(RESERVE < RHO(I,J,R)) THEN
			RESERVE=RHO(I,J,R);MAXRHO=I
		END IF
	END DO
	END DO
100	CONTINUE
	!θの計算
	THETA=0.0;RESERVE2=0.0
	DO 200 I=0,IMAX
	DO J=0,JMAX
		IF(I==IMAX) CYCLE !計算の都合上θ=0としている
		IF(J==JMAX) CYCLE !面はすべて吸収状態
		DO R=1,JMAX                              
			THETA(I,J)=THETA(I,J)+RHO(I,J,R)
		END DO
!	THETA(I,J)=DEXP(THETA(I,J))
			
		!オーバーフローを防ぐため最大のθを調べておく(?)
		IF(RESERVE2 < THETA(I,J)) THEN
			RESERVE2=THETA(I,J);MAXTHETA=I
		END IF
	END DO
200	END DO
	!各Θij-Θlrの差(η)
	ETA=0.0
	DO I=0,IMAX; DO J=0,JMAX
	DO L=0,IMAX; DO R=0,JMAX 
		ETA(I,J,L,R)=THETA(I,J)-THETA(L,R)
	END DO;	END DO; END DO; END DO
	END SUBROUTINE THETAVALUE

****************************************************
*	対数尤度関数
****************************************************
	REAL*8 FUNCTION OBJ(BETA)
	USE VARIABLE; IMPLICIT NONE
	REAL*8, EXTERNAL:: PROB
	REAL*8, INTENT(IN):: BETA(1:MMAX,0:IMAX-1,0:JMAX,1:JMAX)
	INTEGER I,J,L,R,K,M !ランク、タイプ、データ数
	INTEGER EK,FK !元レーティング
	INTEGER IK,JK,LK,RK,V,X,Y
	REAL:: ZK,Z,W,U
	REAL*8 XK(1:MMAX,0:IMAX-1,0:JMAX,1:JMAX)
	REAL*8 PI,Q,QQ
	OPEN(1,file=FILE1)
	READ(1,*) CH
	OBJ=0.0
	Q=0.0
	V=0
	L=0
	DO 101 K=1,TOTALDATA !データの読み込み
		!２時点間の劣化状態、点検間隔、特性ベクトル
		READ(1,*) EK,FK,ZK,(XK(M,0,0,1), M=1,MMAX)
	IF(X==EK .AND. Y==FK .AND. Z==ZK .AND. W==XK(2,0,0,1)
     & .AND.V==XK(3,0,0,1) .AND.U==XK(MMAX,0,0,1)) THEN
	OBJ=OBJ+QQ
	GO TO 101
	END IF
	X=EK; Y=FK; Z=ZK; W=XK(2,0,0,1) 
	V=XK(3,0,0,1); U=XK(MMAX,0,0,1) !説明変数
		DO M=1,MMAX
			XK(M,0,0,2)=XK(M,0,0,1)
			XK(M,0,0,3)=XK(M,0,0,1)
		END DO
		DO I=1,IMAX-1; DO J=1,JMAX; DO R=1,JMAX; DO M=1,MMAX
			XK(M,I,J,R)=XK(M,0,0,1)
		END DO; END DO; END DO; END DO

!--------------------------------------------
!レーティングを劣化ランクとタイプに変換する
!--------------------------------------------
	IF(EK==0)THEN !ランク０は特別
			IK=0; JK=0
		ELSE
			JK=MOD(EK+JMAX-1,JMAX)+1
			IK=(EK-JK)/(JMAX)+1
	END IF
	IF(FK==0)THEN
			LK=0; RK=0
		ELSE
			RK=MOD(FK+JMAX-1,JMAX)+1
			LK=(FK-RK)/(JMAX)+1
		END IF

		!θの値
		CALL THETAVALUE(IK,JK,LK,RK,XK,BETA)
		PI=PROB(IK,JK,LK,RK,ZK)
		!尤度関数
		IF(PI<=1.0D-305) THEN
			OBJ=OBJ-702.0
!	print*,k,"less"
	L=L+1
		ELSE IF(PI>1) THEN
!	print*,k,"more"
		ELSE
		QQ=DLOG(PI)
			OBJ=OBJ+QQ
		END IF
101	END DO
	CLOSE(1)
	RETURN
	END FUNCTION

************************************
*	推移確率行列の要素π_ij,lr(Z)
************************************
	REAL*8 FUNCTION PROB(IK,JK,LK,RK,ZK)
	USE VARIABLE; IMPLICIT NONE
	REAL*8,EXTERNAL:: AIJ
	REAL*8,EXTERNAL:: CIJ
	INTEGER,INTENT(IN):: IK,JK,LK,RK
	REAL,INTENT(IN):: ZK
	INTEGER R,P
	INTEGER REDXMAX
	REAL*8 PROBA,Q,AA
!	REAL*8 BIJ(REDXMAX,REDXMAX)
	REAL*8 BIJ0(XMAX,XMAX)
	REAL*8 BIJ1(XMAX-1*JMAX,XMAX-1*JMAX)
	REAL*8 BIJ2(XMAX-2*JMAX,XMAX-2*JMAX)
		REDXMAX=XMAX-3*IK !行列Ｂの次元を与えておく

	PROB=0.0
	PROBA=0.0
	IF(ZK <= 1.0D-06) RETURN
	!ランク０
	IF(IK==0) THEN
		IF(IK==LK) THEN
	!推移しない場合
		PROBA=0.0
			IF(RK==JK) THEN
				PROBA=DEXP(-THETA(IK,JK)*ZK)
			END IF
		ELSE
	!推移する場合
		CALL TRANS0(AIJ,BIJ0,IK,JK)
		CALL MATRIX_INVERSE(BIJ0,REDXMAX)
	DO P=IK+1,JMAX; DO R=1,JMAX
	PROBA=PROBA+CIJ(P,R,IK,JK,ZK)
     & *BIJ0(3*P+R-3*IK-3,3*LK+RK-3*IK-3)
	END DO; END DO
	END IF
	
	!ランク１
	ELSE IF(IK==1) THEN
		IF(IK==LK) THEN
	!推移しない場合
		PROBA=0.0
			IF(RK==JK) THEN
				PROBA=DEXP(-THETA(IK,JK)*ZK)
			ELSE
			PROBA=0.0
			END IF
		ELSE
	!推移する場合
		CALL TRANS1(AIJ,BIJ1,IK,JK)
		CALL MATRIX_INVERSE(BIJ1,REDXMAX)
	DO P=IK+1,IMAX; DO R=1,JMAX
	PROBA=PROBA+CIJ(P,R,IK,JK,ZK)
     & *BIJ1(3*P+R-3*IK-3,3*LK+RK-3*IK-3)
	END DO; END DO
	END IF

	!ランク２
	ELSE IF(IK==2) THEN
		IF(IK==LK) THEN
	!推移しない場合
		PROBA=0.0
			IF(RK==JK) THEN
				PROBA=DEXP(-THETA(IK,JK)*ZK)
			ELSE
			PROBA=0.0
			END IF
		ELSE
	!推移する場合
		CALL TRANS2(AIJ,BIJ2,IK,JK)
		CALL MATRIX_INVERSE(BIJ2,REDXMAX)
	
	DO P=IK+1,IMAX; DO R=1,JMAX
	PROBA=PROBA+CIJ(P,R,IK,JK,ZK)
     & *BIJ2(3*P+R-3*IK-3,3*LK+RK-3*IK-3)
	END DO; END DO
	END IF

	!ランク３
	ELSE IF(IK==3)THEN
		PROBA=0.0
			IF(RK==JK) THEN
				PROBA=1.0
			END IF
	END IF
	PROB=PROBA
	RETURN; END FUNCTION

**************************************************
*	a^mr_pq(I,J)を求める
**************************************************
	REAL*8 FUNCTION AIJ(P,Q,M,R,I,J)
	USE VARIABLE; IMPLICIT NONE
	INTEGER, INTENT(IN):: P,Q,R,M,I,J
	INTEGER R1,R2
	REAL*8 KAWARI
	AIJ=0.0
	KAWARI=0.0
	IF(Q-P==1) THEN 
	KAWARI=RHO(P,M,R)/ETA(P,M,Q,R)
	END IF

	IF(Q-P==2) THEN  !実際は１パターンしかないが一般化を意識
	DO R2=1,JMAX-1	!これをJMAXに変えると面も考慮できる
	KAWARI=KAWARI+RHO(P,M,R2)*RHO(P+1,R2,R)/ETA(P,M,P+2,R)
     & /ETA(P+1,R2,P+2,R)
	END DO
	END IF
	AIJ=KAWARI
	RETURN; END FUNCTION
	
**************************************************
*	C^r_p(I,J)を求める
**************************************************
	REAL*8 FUNCTION CIJ(P,R,I,J,Z)
	USE VARIABLE; IMPLICIT NONE
	INTEGER, INTENT(IN):: P,R,I,J
	REAL, INTENT(IN):: Z
	INTEGER R1,R2
	REAL*8 KAWARI
	CIJ=0.0
	KAWARI=0.0

	IF(P-I==1) THEN
	KAWARI=RHO(I,J,R)/ETA(I,J,P,R)
	END IF

	IF(P-I==2) THEN
	DO R1=1,JMAX-1	!これをJMAXに変えると面も考慮できる
	KAWARI=KAWARI+RHO(I,J,R1)*RHO(I+1,R1,R)/ETA(I,J,P,R)
     & /ETA(I+1,R1,P,R)
	END DO
	END IF

	IF(P-I==3) THEN
	DO R1=1,JMAX-1; DO R2=1,JMAX-1	!これをJMAXに変えると面も考慮できる
	KAWARI=KAWARI+RHO(I,J,R1)*RHO(I+1,R1,R2)*RHO(I+2,R2,R)
     & /ETA(I,J,P,R)/ETA(I+1,R1,P,R)/ETA(I+2,R2,P,R)
	END DO; END DO
	END IF	
	CIJ=KAWARI*(DEXP(-THETA(P,R)*Z)-DEXP(-THETA(I,J)*Z))
	RETURN; END FUNCTION
	
**************************************************
*	AIJからB(IJ)を求める (IK=0)
**************************************************
	SUBROUTINE TRANS0(AIJ,BIJ0,IK,JK)
	USE VARIABLE; IMPLICIT NONE
	REAL*8,EXTERNAL:: AIJ      !AIJ(P,Q,M,R,I,J)
	INTEGER,INTENT(IN):: IK,JK
	INTEGER REDXMAX
!	REAL*8 BIJ(REDXMAX,REDXMAX) !MAXで９×９の行列
	REAL*8 BIJ0(XMAX,XMAX)
	INTEGER P,Q,M,R,I,J			!IK=3はないので問題なし
	INTEGER X,Y
	BIJ0=0.0
	REDXMAX=XMAX-3*IK !行列Ｂの次元を与えておく
	
	DO P=IK+1,IMAX-1; DO Q=IK+2,IMAX
	DO M=1,JMAX-1; DO R=1,JMAX	!MをJMAXに変えると面も考慮できる
		BIJ0(3*P+M-3-3*IK,3*Q+R-3-3*IK)=AIJ(P,Q,M,R,IK,JK)
	END DO; END DO; END DO; END DO
	DO X=1,REDXMAX
	BIJ0(X,X)=1.0
	END DO
	END SUBROUTINE TRANS0

****************************************
*	AIJからB(IJ)を求める !ちょっと怪しい(IK=1)
****************************************
	SUBROUTINE TRANS1(AIJ,BIJ1,IK,JK)
	USE VARIABLE; IMPLICIT NONE 
	REAL*8,EXTERNAL:: AIJ      !AIJ(P,Q,M,R,I,J)
	INTEGER,INTENT(IN):: IK,JK
	INTEGER REDXMAX
!	REAL*8 BIJ(REDXMAX,REDXMAX) !MAXで９×９の行列
	REAL*8 BIJ1(XMAX-1*IMAX,XMAX-1*JMAX)
	INTEGER P,Q,M,R,I,J			!IK=3はないので問題なし
	INTEGER X,Y
	BIJ1=0.0

	REDXMAX=XMAX-3*IK !行列Ｂの次元を与えておく
	
	DO P=IK+1,IMAX-1; DO Q=IK+2,IMAX
	DO M=1,JMAX-1; DO R=1,JMAX	!MをJMAXに変えると面も考慮できる(?)
		BIJ1(3*P+M-3-3*IK,3*Q+R-3-3*IK)=AIJ(P,Q,M,R,IK,JK)
	END DO; END DO; END DO; END DO
	DO X=1,REDXMAX
	BIJ1(X,X)=1.0
	END DO

	END SUBROUTINE TRANS1

****************************************
*	AIJからB(IJ)を求める !ちょっと怪しい(IK=2)
****************************************
	SUBROUTINE TRANS2(AIJ,BIJ2,IK,JK)
	USE VARIABLE; IMPLICIT NONE 
	REAL*8,EXTERNAL:: AIJ      !AIJ(P,Q,M,R,I,J)
	INTEGER,INTENT(IN):: IK,JK
	INTEGER REDXMAX
!	REAL*8 BIJ(REDXMAX,REDXMAX) !MAXで９×９の行列
	REAL*8 BIJ2(XMAX-2*JMAX,XMAX-2*JMAX)
	INTEGER P,Q,M,R,I,J			!IK=3はないので問題なし
	INTEGER X,Y
	BIJ2=0.0

	REDXMAX=0
	REDXMAX=XMAX-3*IK !行列Ｂの次元を与えておく
	
	DO P=IK+1,IMAX-1; DO Q=IK+2,IMAX
	DO M=1,JMAX-1; DO R=1,JMAX	!MをJMAXに変えると面も考慮できる(?)
		BIJ2(3*P+M-3-3*IK,3*Q+R-3-3*IK)=AIJ(P,Q,M,R,IK,JK)
	END DO; END DO; END DO; END DO

	DO X=1,REDXMAX
	BIJ2(X,X)=1.0
	END DO

	END SUBROUTINE TRANS2

********************************************************
*	BETA(M,I,J,R)をTBETA(M,N)に変換する
********************************************************	
	SUBROUTINE TRANSITION_1(BETA,TBETA)
	USE VARIABLE; IMPLICIT NONE
	REAL*8,INTENT(IN):: BETA(1:MMAX,0:IMAX-1,0:JMAX,1:JMAX)
	REAL*8 TBETA(1:MMAX,1:NMAX)
	INTEGER M,I,J,R

	!ランクが０の場合は特別
	DO M=1,MMAX
	TBETA(M,1)=BETA(M,0,0,1)
	TBETA(M,2)=BETA(M,0,0,2)
	TBETA(M,3)=BETA(M,0,0,3)
	END DO

	DO M=1,MMAX; DO I=1,IMAX-1; DO J=1,JMAX; DO R=1,JMAX

	TBETA(M,9*(I-1)+3*(J-1)+R+3)=BETA(M,I,J,R)

	END DO; END DO; END DO; END DO

	END SUBROUTINE TRANSITION_1

********************************************************
*	TBETA(M,N)をBETA(M,I,J,R)に変換する
********************************************************	
	SUBROUTINE TRANSITION_2(TBETA,BETA)
	USE VARIABLE; IMPLICIT NONE
	REAL*8 BETA(1:MMAX,0:IMAX-1,0:JMAX,1:JMAX)
	REAL*8,INTENT(IN)::TBETA(MMAX,NMAX)
	INTEGER M,I,J,R

	!ランクが０の場合は特別
	DO M=1,MMAX
	BETA(M,0,0,1)=TBETA(M,1)
	BETA(M,0,0,2)=TBETA(M,2)
	BETA(M,0,0,3)=TBETA(M,3)
	END DO

	DO M=1,MMAX; DO I=1,IMAX-1; DO J=1,JMAX; DO R=1,JMAX

	BETA(M,I,J,R)=TBETA(M,9*(I-1)+3*(J-1)+R+3)
	END DO; END DO; END DO; END DO

	END SUBROUTINE TRANSITION_2

*********************************************************
*	プリント用サブルーチン
*********************************************************
	SUBROUTINE PRINT1(BETA,Z)		!ハザード関数，推移確率プリント用
	USE VARIABLE;IMPLICIT NONE
	REAL*8,EXTERNAL:: PROB
	REAL*8,INTENT(IN):: BETA(1:MMAX,0:IMAX-1,0:JMAX,1:JMAX)
	REAL*8 XK(1:MMAX,0:IMAX-1,0:JMAX,1:JMAX)
	INTEGER I,J,M,L,R
	REAL Z
	XK(1,0,0,1)=1.0
!	XK(2,0,0,1)=0.676864432	!構造平均
!	XK(3,0,0,1)=0.488026586	!車線平均
!	XK(3,0,0,1)=0.621239028	!交通量平均
!	XK(2,0,0,1)=0.317173297	!上下平均
!	XK(2,0,0,1)=0.0
!	XK(3,0,0,1)=1.0
!	XK(4,0,0,1)=0.370308829
!	XK(4,0,0,1)=1.0

		DO M=1,MMAX
			XK(M,0,0,2)=XK(M,0,0,1)
			XK(M,0,0,3)=XK(M,0,0,1)
		END DO
		DO I=1,IMAX-1; DO J=1,JMAX; DO R=1,JMAX; DO M=1,MMAX
			XK(M,I,J,R)=XK(M,0,0,1)
		END DO; END DO; END DO; END DO

	CALL THETAVALUE(I,J,L,R,XK,BETA)
			
	OPEN(1002,file=FILE1002)
	!ρ
		PRINT 10, (RHO(0,0,R),R=1,JMAX)
		WRITE(1002,10) (RHO(0,0,R),R=1,JMAX)

	DO I=1,IMAX; DO J=1,JMAX 
		PRINT 10, (RHO(I,J,R),R=1,JMAX)
		WRITE(1002,10) (RHO(I,J,R),R=1,JMAX)
	END DO; END DO

	!Θ
	PRINT *, "ハザード関数"
	WRITE(1002,*) "ハザード関数"
	PRINT 10, (THETA(0,0))
	WRITE(1002,10) (THETA(0,0))
	DO I=1,IMAX
	PRINT 10, (THETA(I,J),J=1,JMAX)
	WRITE(1002,10) (THETA(I,J),J=1,JMAX)
	END DO
	
 10	FORMAT(1X,6(F10.4,","),F10.4)
	CLOSE(1002)
	
	PRINT *,"期待寿命"
	DO I=0,IMAX
	PRINT 10, (1.0/THETA(I,J),J=1,JMAX)
	END DO
	!呼び出すとΘの値が変化してしまう
!	CALL ELIFE(BETA)
	END SUBROUTINE
	

!*************************************************************
!*	平均レーティング期待寿命							まだやってない
!*************************************************************
!	SUBROUTINE ELIFE(BETA)
!	USE VARIABLE;IMPLICIT NONE
!	REAL*8,EXTERNAL:: PROB
!	REAL*8,INTENT(IN):: BETA(MMAX,JMAX-1)
!	INTEGER IK,JK
!	INTEGER J,K,M
!	REAL*8 XK(MMAX,JMAX-1)
!	REAL*8 LIFE(JMAX-1)
!	REAL ZK
!
!	OPEN(1,file=FILE1)
!	OPEN(1003,file=FILE1003)
!	READ(1,*) CH
!	
!	LIFE=0.0
!	DO K=1,TOTALDATA
!		!２時点間の劣化状態，点検期間，特性ベクトル
!		READ(1,*) IK,JK,ZK,(XK(M,1), M=1,MMAX)
!		DO J=2,JMAX-1;DO M=1,MMAX
!			XK(M,J)=XK(M,1)
!		END DO; END DO
!		CALL THETAVALUE(IK,JK,LK,RK,XK,BETA)
!		DO J=1,JMAX-1
!			LIFE(J)=LIFE(J)+1.0/THETA(J)
!		END DO
!	END DO
!	LIFE=LIFE/TOTALDATA
!	PRINT *, "平均期待寿命"
!	PRINT *, (LIFE(J),J=1,JMAX-1)
!	WRITE(1003,10) (LIFE(J),J=1,JMAX-1)
!  10	FORMAT(1X,6(F10.4,","),F10.4)
!
!	CLOSE(1)
!	RETURN;	END SUBROUTINE
*************************************************************
*	行列作成（要素がすべて0の行or列を削除）
*************************************************************
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

**************************************************
*	LU分解を用いた連立方程式の解法 Ax=b	 N次元
**************************************************
      SUBROUTINE MATRIX_1(N,A,b,X)
      REAL*8 A(N,N),W(N),X(N)
	REAL*8,INTENT(IN):: b(N)
      INTEGER IP(N)
      DATA EPS/1.0E-75/

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
      PRINT *, 'この行列特異';STOP
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
*	逆行列(完全ピボット法)
*******************************************************
	SUBROUTINE MATRIX_INVERSE(A,N)
	IMPLICIT NONE
	INTEGER,INTENT(IN):: N
	REAL*8:: A(N,N),AW,MM(N)
	INTEGER:: K,I,J,P,Q
	REAL*8 EPS
	INTEGER LW,MW,L(N),M(N)

	EPS=1.0D-30

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
			PRINT *, "この行列特異" ;STOP
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
	
	RETURN;	END SUBROUTINE