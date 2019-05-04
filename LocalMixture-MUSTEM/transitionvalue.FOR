************************************************************************
*	マルコフ推移確率の計算
************************************************************************
	MODULE VARIABLE
		INTEGER,PARAMETER:: JMAX=5,MMAX=3  !劣化状態、特性数
	!-----------------------------------------
		REAL*8:: THETA(JMAX) !ハザード関数Θj
		REAL*8:: THETASA(JMAX,JMAX)	!Θj-Θi
		INTEGER:: MAXTHETA
	END MODULE
***********************************************************************
*	メインプログラム
***********************************************************************
	PROGRAM ESTIMATION
	USE VARIABLE; IMPLICIT NONE

	REAL*8,PARAMETER:: EPS1=1.0D-06
	REAL Z
	INTEGER I,J

!---------------------------------------------------------------
*	入力データ
!	ハザード率
	THETA(1)=0.2478;THETA(2)=0.5115;THETA(3)=1.2847;THETA(4)=0.901;
     #	THETA(JMAX)=0.0
!	評価年数
	Z=1.0
!---------------------------------------------------------------
	PRINT *,'期待寿命',1.0/THETA(1)+1.0/THETA(2)+1.0/THETA(3)

	!各Θi-Θjの差
	THETASA=0.0
	DO J=1,JMAX; DO I=1,JMAX
		THETASA(I,J)=THETA(I)-THETA(J)
	END DO;	END DO
	CALL PRINT1(Z)

	STOP
	END PROGRAM

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
	SUBROUTINE PRINT1(Z)		!ハザード関数，推移確率プリント用
	USE VARIABLE;IMPLICIT NONE
	REAL*8,EXTERNAL:: PROB

	REAL*8 BETA(MMAX,JMAX-1)
	INTEGER I,J,M
	REAL Z
	REAL*8 XK(MMAX,JMAX-1)	
			
	OPEN(1002,file="transition.csv")
	
	!π_ij
	PRINT *, "推移確率********************************"
	WRITE(1002,*) "推移確率********************************"
	DO I=1,JMAX
		PRINT 10, (PROB(I,J,Z),J=1,JMAX)
		WRITE(1002,10) (PROB(I,J,Z),J=1,JMAX)
	END DO
 10	FORMAT(1X,4(F10.4,","),F10.4)
	CLOSE(1002)

	END SUBROUTINE
