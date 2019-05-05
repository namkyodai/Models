************************************************************************
*	橋梁劣化予測のためのマルコフ推移確率の推計(論文作成時)
************************************************************************
	MODULE VARIABLE
		INTEGER,PARAMETER:: JMAX=5,MMAX=2  !劣化状態、特性数
	!------------------------------------------
		INTEGER,PARAMETER:: TOTALDATA=17384	!使用データ数
	!------------------------------------------
		INTEGER,PARAMETER:: TIMEMAX=100	!繰り返し回数
		INTEGER,PARAMETER:: dimN=(JMAX-1)*MMAX	!連立方程式の数

	!---推定パラメータ数を変更したい場合-------------
		INTEGER,PARAMETER:: D1=0	!削除パラメータ数
		INTEGER,PARAMETER:: redimN=dimN-D1	!修正後パラメータ数
		INTEGER DEL(D1)
		DATA DEL/0/	!削除するパラメータ番号((J-1)*MMAX+M)
		!初期パラメータを操作(BETA(J)=J1*J+B1)
	 REAL*8,PARAMETER:: B1=1.1D-02,J1=3.0D-03	!ベター
	!	REAL*8,PARAMETER:: B1=5.5D-03,J1=2.5D-04	
	!-----------------------------------------
		!データのレコード数	 MAX32902
		CHARACTER(LEN=30),PARAMETER:: FILE1="input\inputdata.csv"
		!データのレコード数	 MAX26748

		CHARACTER(LEN=13),PARAMETER:: FILE1000="output\df.csv"
		CHARACTER(LEN=17),PARAMETER:: FILE1001="output\beta-t.csv"
		CHARACTER(LEN=17),PARAMETER:: FILE1002="output\hazard.csv"
		CHARACTER(LEN=17),PARAMETER:: FILE1003="output\avlife.csv"

		REAL*8:: THETA(JMAX) !ハザード関数Θj
		REAL*8:: THETASA(JMAX,JMAX)	!Θj-Θi
		INTEGER:: MAXTHETA
		CHARACTER CH

		REAL*8 X1(MMAX)
		REAL Z1
		!DATA X1/1.0, 0.226569, 0.043135/	!全データサンプルの平均
		!DATA X1/1.0, 0.237009, 0.049857/	!対象径間の平均(2年間)
		DATA X1/1.0, 0.204445/	!比較分析用 204445
		DATA Z1/1.0/
	END MODULE
***********************************************************************
*	メインプログラム
***********************************************************************
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

*----------------------------
*	最尤推定法
*----------------------------
!	PRINT *, "交通量は？"
!	READ(*,*) X1(2)
!	PRINT *, "床版面積は？"
!	READ(*,*) X1(3)
!-----初期パラメータβ------------------------------
	DO J=1,JMAX-1;	DO M=1,MMAX
		BETA(M,J)=J*J1+B1	!ベター
	!	BETA(M,J)=1.0D-06+1.1D-05*J
	END DO;	END DO
	
	DO I=1,D1
		M=MOD(DEL(I)+(MMAX-1),MMAX)+1;J=(DEL(I)-M)/MMAX+1
		BETA(M,J)=0.0
	END DO
!----------------------------------------------------
	OPTBETA=BETA

	DO 100 TIME=1,TIMEMAX
		!ニュートン・ラプソン法を用いた連立非線形方程式
		!近似解βを求める
		CALL NEWTON_RAPHSON(OPTBETA,G)
	
		EVAL=1	!EVALが1のままならベータが収束している
		DO I=1,JMAX-1;DO M=1,MMAX
			!誤差が大きければベータの再計算
			IF(DABS( (OPTBETA(M,I)-BETA(M,I))/OPTBETA(M,I) ) > EPS1) THEN
				EVAL=0;GOTO 200
			END IF
		END DO; END DO
 200		CONTINUE	!IF文からの脱出

		!βの更新
		BETA=OPTBETA
		IF(EVAL==1)	EXIT	!βが収束なら脱出
		
		IF(MOD(TIME,10)==0) THEN
		PRINT *, TIME ,"回目です"
			DO I=1,JMAX-1
				PRINT *, (BETA(M,I),M=1,MMAX)
			END DO
		END IF
		
		IF(TIME==TIMEMAX) THEN
			PRINT *, "繰り返し計算収束せず";STOP
		END IF
 100	CONTINUE

*-------------------
*	検定
*-------------------
	CALL DEL_MATRIX(G,A)
	A=-A; CALL MATRIX_INVERSE(A,redimN)	!フィッシャー行列の逆行列
!	検定
	OPTOBJ=OBJ(OPTBETA)
	PRINT *, "尤度関数"
	PRINT *, OPTOBJ
	OPEN(1001,file=FILE1001)
	
	!パラメータβ
	PRINT *, "未知パラメータβ********************"
	WRITE(1001,*) "未知パラメータβ***************"
	DO I=1,JMAX-1
		PRINT *, (OPTBETA(M,I),M=1,MMAX)
		WRITE(1001,20) (OPTBETA(M,I),M=1,MMAX)
	END DO 
 
	PRINT *, "t値********"
	WRITE(1001,*) "t値********"
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

	!ハザード関数・推移確率行列・未知パラメータβ
	Z=Z1
	CALL PRINT1(OPTBETA,Z)

	STOP
	END PROGRAM

*******************************************************
*	ニュートン・ラプソン法による連立非線形方程式
*******************************************************
	SUBROUTINE NEWTON_RAPHSON(BETA,HES)
	USE VARIABLE; IMPLICIT NONE

	REAL*8:: BETA(MMAX,JMAX-1)
	REAL*8:: G(dimN,dimN),HES(dimN,dimN),F(dimN),C(redimN)
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
		!２時点間の劣化状態，点検期間，特性ベクトル
		READ(1,*) IK,JK,ZK,(XK(M,1), M=1,MMAX)
		DO J=2,JMAX-1;DO M=1,MMAX
			XK(M,J)=XK(M,1)
		END DO; END DO
		DO I=1,D1
			M=MOD(DEL(I)+(MMAX-1),MMAX)+1;J=(DEL(I)-M)/MMAX+1
			XK(M,J)=0.0
		END DO

		!Θの値
		CALL THETAVALUE(IK,JK,XK,BETA)
		
		!πの１階偏微分　dlnπ_(ij)/dΘ_lを求める
		CALL FIRST_DIFFERENTIAL(IK,JK,ZK,df)
			
		!１階偏微分方程式
		DO L=IK,JK
			IF(L==JMAX) CYCLE
			DO M=1,MMAX
				F((L-1)*MMAX+M)=F((L-1)*MMAX+M)+df(L)*XK(M,L)
 			END DO
		END DO

		!２階偏微分
		CALL HESSIAN(IK,JK,ZK,df2,df)

		!ヘッシアンマトリクス
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

	!GC=Fを解く
	!パラメータ数を部分的に削った場合
	!その分の行，列を削除しておく必要がある
	L=0
	DO 410 I=1,dimN
		DO J=1,D1
			IF(I==DEL(J)) GOTO 410
		END DO
		L=L+1; X(L)=F(I)
 410	CONTINUE

	CALL DEL_MATRIX(G,A)
	HES=G

	CALL MATRIX_1(redimN,A,X,C)

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

********************************************
*	１階偏微分サブルーチン	  df(L)=dlnπ_(IK,JK)/dΘ_L
********************************************
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
		IF(L==JMAX) CYCLE !THETA(JMAX)は定義されていないので微分も不可
		SUM2=0.0	!分母

		DO H=IK,JK
			SUM1=0.0	!分子
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
*	２階偏微分サブルーチン（ヘッセ行列）
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

	!非対角要素
	DO 100 L=IK,JK-1
	DO 200 N=L+1,JK
		IF(N==JMAX) CYCLE	!THETA(JMAX)は定義されていないので微分も不可
		SUM2=0.0	!分母
		
		DO H=IK,JK
			SUM1=0.0	!分子
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

		HES(N,L)=HES(L,N)	!対称行列
 200	CONTINUE
 100	CONTINUE
	
	!対角要素
	DO 300 L=IK,JK
		IF(L==JMAX) CYCLE	!THETA(JMAX)は定義されていないので微分も不可
		SUM2=0.0	!分母
		DO H=IK,JK
			SUM1=0.0	!分子
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

*******************************************
*	Θi(i=IK,JK)の値を求める
*******************************************
	SUBROUTINE THETAVALUE(IK,JK,XK,BETA)
	USE VARIABLE; IMPLICIT NONE
	
	REAL*8,INTENT(IN):: BETA(MMAX,JMAX-1),XK(MMAX,JMAX-1)
	INTEGER,INTENT(IN):: IK,JK
	INTEGER I,J,M
	REAL*8 RESERVE

	!Θの計算
	THETA=0.0
	MAXTHETA=IK;RESERVE=0.0
	DO J=IK,JK
		IF(J==JMAX) CYCLE
		DO M=1,MMAX
			THETA(J)=THETA(J)+BETA(M,J)*XK(M,J)
		END DO
	!オーバーフローを防ぐため最大のΘを調べておく
		IF(RESERVE < THETA(J)) THEN		
			RESERVE=THETA(J);MAXTHETA=J
		END IF
	END DO
	
	!各Θi-Θjの差
	THETASA=0.0
	DO J=IK,JK; DO I=IK,JK
		THETASA(I,J)=THETA(I)-THETA(J)
	END DO;	END DO

	END SUBROUTINE

**********************************
*	関数
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

**************************************************
*	尤度関数
**************************************************
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
		!ファイルからの読み込み
		!２時点間の劣化状態，点検期間，特性ベクトル
		READ(1,*) IK,JK,ZK,(XK(M,1), M=1,MMAX)
		DO J=2,JMAX-1;DO M=1,MMAX
			XK(M,J)=XK(M,1)
		END DO; END DO
		DO I=1,D1
			M=MOD(DEL(I)+(MMAX-1),MMAX)+1;J=(DEL(I)-M)/MMAX+1
			XK(M,J)=0.0
		END DO
		!Θの値
		CALL THETAVALUE(IK,JK,XK,BETA)

		PI=PROB(IK,JK,ZK)
		!尤度関数
		IF(PI<=1.0D-305) THEN	!Θがマイナスになるケースでは
								!πが負になるケースも
			OBJ=OBJ-702.2884
		ELSE
			OBJ=OBJ+DLOG(PI)
		END IF
	END DO
	CLOSE(1)
	RETURN
	END FUNCTION

************************************
*	推移確率行列の要素π_ij(Z)
************************************
	REAL*8 FUNCTION PROB(I,J,Z)
	USE VARIABLE;IMPLICIT NONE
	REAL*8,EXTERNAL:: PROD1
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

*********************************************************
*	プリント用サブルーチン
*********************************************************
	SUBROUTINE PRINT1(BETA,Z)		!ハザード関数，推移確率プリント用
	USE VARIABLE;IMPLICIT NONE
	REAL*8,EXTERNAL:: PROB

	REAL*8 BETA(MMAX,JMAX-1)
	INTEGER I,J,M
	REAL Z
	REAL*8 XK(MMAX,JMAX-1)
	

	DO J=1,JMAX-1;	DO M=1,MMAX
		XK(M,J)=X1(M)
	END DO;END DO

	CALL THETAVALUE(1,JMAX,XK,BETA)	
			
	OPEN(1002,file=FILE1002)
	
	!Θ
	PRINT *, "ハザード関数****************************"
	WRITE(1002,*) "ハザード関数***********"
	PRINT 10, (THETA(J),J=1,JMAX)
	WRITE(1002,10) (THETA(J),J=1,JMAX)
	!π_ij
	PRINT *, "推移確率********************************"
	WRITE(1002,*) "推移確率********************************"
	DO I=1,JMAX
		PRINT 10, (PROB(I,J,Z),J=1,JMAX)
		WRITE(1002,10) (PROB(I,J,Z),J=1,JMAX)
	END DO
 10	FORMAT(1X,6(F10.4,","),F10.4)
	CLOSE(1002)
	
	PRINT *,"期待寿命"
	PRINT 10, (1.0/THETA(J),J=1,JMAX)
	!呼び出すとΘの値が変化してしまう
	CALL ELIFE(BETA)
	END SUBROUTINE
*************************************************************
*	平均レーティング期待寿命
*************************************************************
	SUBROUTINE ELIFE(BETA)
	USE VARIABLE;IMPLICIT NONE
	REAL*8,EXTERNAL:: PROB
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
		!２時点間の劣化状態，点検期間，特性ベクトル
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
	PRINT *, "平均期待寿命"
	PRINT *, (LIFE(J),J=1,JMAX-1)
	WRITE(1003,10) (LIFE(J),J=1,JMAX-1)
  10	FORMAT(1X,6(F10.4,","),F10.4)

	CLOSE(1)
	RETURN;	END SUBROUTINE
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