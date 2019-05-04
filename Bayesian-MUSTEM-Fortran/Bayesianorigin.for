*******************************************************
*     �}���R�t���ڊm���̃x�C�Y����(Gibbs�T���v�����O)
*	(���K���z)
*******************************************************
	MODULE FILE
	CHARACTER(LEN=30),PARAMETER:: INPUT_F1="input\data.csv"
		CHARACTER(LEN=30),PARAMETER:: INPUT_F2="input\noinfo.csv"
		CHARACTER(LEN=30),PARAMETER:: OUTPUT_F1="output\point.csv"
		CHARACTER(LEN=30),PARAMETER:: OUTPUT_F2="output\result.csv"
	END MODULE
!------------------------------------------------------------------------
	MODULE VARIABLE
		INTEGER,PARAMETER:: NDATA=880	!�f�[�^��
		INTEGER,PARAMETER:: JMAX=6	!��ԃ����N��(J-1)
		INTEGER,PARAMETER:: MMAX=3	!�����ϐ��̐�
		REAL*8 BETA(JMAX,MMAX)
		INTEGER,PARAMETER:: D1=2	!�폜�p�����[�^��
		INTEGER DEL(D1)
		DATA DEL/2,6/	!�폜����p�����[�^�ԍ�((J-1)*MMAX+M)
		INTEGER IK(NDATA),JK(NDATA)
		REAL*8 XK(MMAX,NDATA), ZK(NDATA), LAMB(JMAX,NDATA)
		REAL*8 MYU0(MMAX,JMAX),ROH(MMAX,MMAX,JMAX)
		INTEGER PATTERN
		INTEGER,PARAMETER:: LIMIT=1
	END MODULE
**********************************************
*	Gibbs�T���v�����O�ɕK�v�Ȓl
**********************************************
	SUBROUTINE FIRST(IBURN,KMAX,NT,SXK,JM,MM)
	INTEGER,INTENT(OUT):: IBURN	!�}���R�t�`�F�[����������
*	ARS���x�z����p�����[�^
	INTEGER,INTENT(IN):: KMAX,JM,MM	!��������W��
	INTEGER,INTENT(OUT):: NT	!Tk�Fx���W�̍ő吔
	REAL*8,INTENT(OUT):: SXK(JM*MM,KMAX)	!�X�^�[�g�|�C���g
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
	INTEGER,PARAMETER:: NPOINT=10000 !��]�T���v�����O
	REAL*8 SAMP(NPOINT,JMAX*MMAX)	!�T���v�����O�W��
!---------------------------------------------------------
*	�����̎�
	INTEGER IR1,IX,IY,IZ,TANE
	COMMON /IRAND1/IR1
	COMMON /IRAND2/IX,IY,IZ
!---------------------------------------------------------
	CHARACTER(LEN=30):: hour
	INTEGER I

	CALL TIME(hour); PRINT *, hour
	MYU0=0.0; ROH=0.0

*	�����̎�
	CALL CLOCK (IR1,2)
	IX=16251; IY=21157; IZ=3583
	TANE=IR1

*	�M�u�X�T���v�����O�����l
	BETA(1,1)=-1.0	; BETA(1,2)=0.0	;	BETA(1,3)=0.015
	BETA(2,1)=-2.0	; BETA(2,2)=0.0	;	BETA(2,3)=0.016
	BETA(3,1)=-2.1	; BETA(3,2)=0.0	;	BETA(3,3)=0.0
	BETA(4,1)=-1.1	; BETA(4,2)=0.0	;	BETA(4,3)=0.018
	BETA(5,1)=-1.2	; BETA(5,2)=0.0	;	BETA(5,3)=0.0
	BETA(6,1)=-1.5	; BETA(6,2)=0.0	;	BETA(6,3)=0.0

	DO I=1,1000
	PRINT *,"�S���O���(1) or �����I���O���(2) or �����(3)"
		READ *, PATTERN
		SELECT CASE(PATTERN)
			CASE(1:3); EXIT
		END SELECT
	END DO
!	PATTERN=1
*	�f�[�^�t�@�C������̓���
	CALL INPUT()
	CALL MCMC_G(SAMP,NPOINT)
	
!	�������ւ̕ϊ�
	IF(PATTERN/=3) CALL REVERSE(2)

	CALL OUTPUT(TANE,SAMP,NPOINT)
!	CALL BEEPQQ (3000, 500); CALL BEEPQQ (4000, 500)
	CALL TIME(hour); PRINT *, hour
	STOP;	END PROGRAM
**********************************************
*	�֐�H(x)�CH'(x)�̕]��  H(x)=Ln G(x)
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

		DO I=1,D1	! �p�����[�^��Ԃ��珜��
			M = MOD(DEL(I)+(MMAX-1),MMAX) + 1; J = (DEL(I)-M)/MMAX+1
			XKK(M,J)=0.0
		END DO

!----------------------------------------
*	���ʐݒ�i�v�Z�ȗ͉��j
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
		! �ŏ��l�̒T��
		MT = THETA(IKK)
		DO J=IKK+1,JKK
			IF(THETA(J) < MT) MT=THETA(J)
		END DO
		DO L=IKK,JKK	!�ŏ��̃���0�ɂȂ�悤�ɒ���
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
			PRINT *, K, "�l���}�C�i�X�DWHY ? "; STOP
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
*	�C���v�b�g�t�@�C������̓���
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
!		�X�P�[������
	DO K=1,NDATA
		DO M1=2,MMAX; XK(M1,K)=XK(M1,K)/MAXX(M1); END DO
	END DO
	XK=XK*SCALE

	IF(PATTERN==3) RETURN	!�����̏ꍇ

	OPEN(1001,file=INPUT_F2)
		DO J=1,JMAX
			READ(1001,*) CH
			READ(1001,*) CH,(MYU0(M1,J),M1=1,MMAX)
			DO M1=1,MMAX
				READ(1001,*) CH,(ROH(M1,M2,J),M2=1,MMAX)
			END DO
		END DO
	CLOSE(1001)

	PRINT *, "���O���[����]"
	DO J=1,JMAX
		PRINT 30, J,(MYU0(M1,J),M1=1,MMAX)
	END DO
30	FORMAT(1X,I2,6F10.4)

!		�t�s��֕ϊ��E�X�P�[������
	CALL REVERSE(1)

	RETURN
	END SUBROUTINE INPUT
***************************************************
*	�A�E�g�v�b�g
***************************************************
	SUBROUTINE OUTPUT(TANE,SAMP,NPOINT)
	USE VARIABLE; USE FILE; IMPLICIT NONE
	INTEGER,INTENT(IN):: TANE,NPOINT
	REAL*8,INTENT(IN):: SAMP(NPOINT,JMAX*MMAX)	!�T���v�����O�W��
	REAL*8 SAMP1(NPOINT)	!�W��
	REAL*8 MYU(MMAX,JMAX),SS(MMAX,MMAX,JMAX),SM	!���ρC���U
	REAL*8 UBAR(JMAX*MMAX),OBAR(JMAX*MMAX)	!�M�����
	REAL Z(JMAX*MMAX); COMMON /GEWEKE/Z	!GEWEKE�������蓝�v��
	INTEGER I,J,K,M1,M2
	
!--------------------------
*	�T���v�����O�̏o��
!--------------------------
	OPEN(2000,file=OUTPUT_F1)
		WRITE(2000,*) "�����̎�",",",TANE
		WRITE(2000,*) "�T���v����",",",NPOINT
		WRITE(2000,*) "���1",",,,","���2",",,,","���3",",,,","���4"
		WRITE(2000,*) ("B1",",","B2",",","B3",",",J=1,4)
		DO K=1,NPOINT
			WRITE(2000,20) (SAMP(K,J),J=1,JMAX*MMAX)
		END DO
	CLOSE(2000)
 20	FORMAT(1X,30(E15.8,","))

!--------------------------
*	���v�ʂ̏o��
!--------------------------
!	�W�{����
	MYU=0.0
	DO J=1,JMAX; DO M1=1,MMAX
		DO K=1,NPOINT
			MYU(M1,J)=MYU(M1,J)+SAMP(K,(J-1)*MMAX+M1)
		END DO
	END DO;	END DO	
	MYU=MYU/NPOINT

!	�W�{���U
	DO J=1, JMAX; DO M2=1,MMAX; DO M1=1,MMAX
		I = (J-1)*MMAX
		DO K=1,NPOINT
			SS(M1,M2,J)=SS(M1,M2,J)
     &			+(SAMP(K,I+M1)-MYU(M1,J))*(SAMP(K,I+M2)-MYU(M2,J))
		END DO
	END DO; END DO; END DO
	SS=SS/NPOINT

!	�M�����
	DO J=1,JMAX*MMAX
		DO K=1,NPOINT
			SAMP1(K)=SAMP(K,J)
		END DO
		CALL BUBBLE(SAMP1,NPOINT)
		UBAR(J)=SAMP1(NPOINT*5/100)
		OBAR(J)=SAMP1(NPOINT*95/100)
	END DO
 90	FORMAT(1X, F10.2, A12)

!	�t�@�C���o��	
	OPEN(3000,file=OUTPUT_F2)
		SELECT CASE(PATTERN)
			CASE(1); WRITE(3000,*) "�S���O���"
			CASE(2); WRITE(3000,*) "�������O���"
			CASE(3); WRITE(3000,*) "�����"
		END SELECT
	
		IF(PATTERN/=3) THEN
			WRITE(3000,*) "�������"
			WRITE(3000,*) 
     &		",", "���1",",,,","���2",",,,","���3",",,,","���4",",,,"
			WRITE(3000,10) "��������",((MYU0(M1,J),M1=1,MMAX),J=1,JMAX)
			DO M1=1,MMAX
				WRITE(3000,10) "�������U",((ROH(M1,M2,J),M2=1,MMAX),J=1,JMAX)
			END DO
		END IF
				
		WRITE(3000,*) "������"
		WRITE(3000,*) 
     &		",", "���1",",,,","���2",",,,","���3",",,,","���4",",,,"
		WRITE(3000,10) "�W�{����",((MYU(M1,J),M1=1,MMAX),J=1,JMAX)
		DO M1=1,MMAX
			WRITE(3000,10) "�W�{���U",((SS(M1,M2,J),M2=1,MMAX),J=1,JMAX)
		END DO

		WRITE(3000,10) "�M�����",(UBAR(J),J=1,JMAX*MMAX)
		WRITE(3000,10) "�M�����",(OBAR(J),J=1,JMAX*MMAX)
		WRITE(3000,10) "��������",(Z(J),J=1,JMAX*MMAX)
 10	FORMAT(1X,A8,",",20(E15.8,","))
	CLOSE(3000)

	RETURN;	END SUBROUTINE OUTPUT
**********************************************
*	�X�P�[������
**********************************************
	SUBROUTINE REVERSE(IP)
	USE VARIABLE; IMPLICIT NONE
	INTEGER,INTENT(IN):: IP
	REAL*8 A(MMAX,MMAX)
	REAL*8 MAXX(MMAX)	!�����ϐ��̍ő�l
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
	ELSE	! �t����
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
*	�M�u�X�T���v�����O��p����MCMC�@
**************************************************************
	SUBROUTINE MCMC_G(SAMP,NPOINT)
	USE VARIABLE; IMPLICIT NONE
	INTEGER,INTENT(IN):: NPOINT ! ��]�T���v�����O��
	REAL*8,INTENT(OUT):: SAMP(NPOINT,JMAX*MMAX)	!�T���v�����O�W��
	REAL*8 SAMP1(NPOINT)		!SUB�T���v�����O�W��
	INTEGER IBURN	
*	ARS���x�z����p�����[�^
	INTEGER,PARAMETER:: KMAX=2	!��������W��
	INTEGER NT	!Tk�Fx���W�̍ő吔
	REAL*8 SXK(JMAX*MMAX,KMAX)	!�X�^�[�g�|�C���g
!-----------------------------------------
	LOGICAL OB(JMAX,MMAX),UB(JMAX,MMAX)	!�̈�D�̏�������̑���
	REAL*8 OX(JMAX,MMAX),UX(JMAX,MMAX)	!�̈�D�̏������
	REAL,EXTERNAL:: GEWEKEZ
	REAL Z(JMAX*MMAX); COMMON /GEWEKE/Z	!GEWEKE�������蓝�v��
	REAL*8 MAXX(MMAX), SCALE
	COMMON /ARR/MAXX,SCALE
	LOGICAL KENTEI
	CHARACTER(LEN=30):: hour
!-----------------------------------------
	REAL*8 X	!�߂�l
!-----------------------------------------
	INTEGER IPOINT
	INTEGER I,J,M,T,ID, IREV,d
	LOGICAL STAY

*	ARS�����̐ݒ�
	OB=.FALSE.; UB=.FALSE.
	CALL FIRST(IBURN,KMAX,NT,SXK,JMAX,MMAX)
!*********************
*	Burn-in
!*********************
		PRINT *, "Burn-in����"
		DO T=1,IBURN
		I=1
			IF(MOD(T,1000)==0.0) PRINT *, T,"���"
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
		PRINT *, "Sampling����"
		DO T=IPOINT,NPOINT
		I=1
			IF(MOD(T,1000)==0.0) PRINT *, T,"���"
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
*	��������
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
			IF(Z(ID) < 1.96) THEN	!�T���L�Ӑ���
				PRINT 90, ID,Z(ID), "����"
			ELSE
				PRINT 90, ID,Z(ID), "��������"
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
*	�����葱��
!---------------------------------
	CONTAINS
	SUBROUTINE GIBBS(XJ,UBJ,OBJ,UXJ,OXJ)
		LOGICAL,INTENT(IN):: UBJ,OBJ
		REAL*8,INTENT(IN):: UXJ,OXJ
		REAL*8,INTENT(OUT):: XJ
		INTEGER IFAULT,I,K
		REAL*8 XX(KMAX)	!�������W
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
*	�K���I���p�T���v�����O
*******************************************************
	SUBROUTINE REJECT_SAMPLE(XPOINT,NT,NK,UB,OB,UX,OX,XX,NO,IFAULT)
	IMPLICIT NONE
	INTEGER,INTENT(IN):: NO	!�ǂ̏����t���m�����z�𔭐������邩
	INTEGER,INTENT(IN):: NK
	INTEGER IFAULT	!�G���[����
	REAL*8,EXTERNAL:: RANSU1
	REAL*8,PARAMETER:: EPS10=1.0D-16	!��
	REAL*8 XPOINT, XX(NK)
	INTEGER,INTENT(IN):: NT	!x���W�̍ő吔
	LOGICAL,INTENT(IN):: OB, UB	!�̈�D�̏�������̑���
	REAL*8,INTENT(IN):: UX,OX	!�����C���
	REAL*8 HK(NT),DHK(NT),ZK(0:NT),HZK(0:NT),XKK(NT)
	INTEGER,PARAMETER:: IEMAX=700	!EXP�̍ő����
	INTEGER,PARAMETER:: NIT=1000	!�J�Ԃ���
	REAL*8 RX	!�������W
	REAL*8 W	!��l����
	REAL*8 LK,UK,HX,RH,RDH,RZ1,RZ2,RZH1,RZH2
	INTEGER RJ	!�}���ꏊ
	INTEGER DR	!TEST����
	REAL*8 HZMAX	!H�̃X�P�[������
	REAL*8 RESERVE
	INTEGER K
	INTEGER I,IT

!---------------------------
!	Initialization step
!---------------------------
	K=NK
	DO I=1,K
		XKK(I)=XX(I)	!��������W Tk
		XPOINT=XKK(I)
		CALL MAKEH(XPOINT,HK(I),DHK(I),NO)
	END DO

	IF(UB==.FALSE. .AND. DHK(1) < 0.0) THEN
		PRINT *, "ID=",NO,"�X��",DHK(1),"���͈͂������ƍL��"
		IFAULT=1; RETURN
	END IF

	IF(OB==.FALSE. .AND. DHK(K) > 0.0) THEN
		PRINT *, "ID=",NO,"�X��",DHK(K),"�E�͈͂������ƍL��"
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
!	Sk���痐������
		CALL SKSAMPLE(RX,RJ,K,IFAULT)
		IF(IFAULT==6) RETURN
!	��l�����𔭐�
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
				PRINT *, "ID=",NO, "Tk�������ς��ł�"; IFAULT=3; RETURN
			END IF
!	Tk,Zj�������ɕ��בւ�
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

	PRINT *, "�K�v�ȃT���v�����O�����W�܂炸"
	RETURN
!---------------------------------------------------------
*	�����葱��
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
*	Sk����̗�������
!	�⊶�Ȃ���Gilks�̃v���O������傢�ɎQ��
**************************************************
	SUBROUTINE SKSAMPLE(X,JJ,KMAX,IFAULT)
	REAL*8,EXTERNAL:: RANSU2
!	���߂����l
	INTEGER,INTENT(OUT):: JJ	!J:�}���ꏊ
	REAL*8,INTENT(OUT):: X	!X:�T���v�����W
	INTEGER,INTENT(IN):: KMAX
	INTEGER,INTENT(OUT):: IFAULT
	REAL*8 SSUM(KMAX)	!Uk���̍��v�ʐ�
	REAL*8 LSUM
	REAL*8 EH,SIGN	!�l�ۊǗp
	LOGICAL HORIZ
	REAL*8 R	!����
	INTEGER J

!	����Uk���̊e�ʐς����߂�
	CALL UAREA(KMAX,SSUM)
	IF(SSUM(KMAX)<=0.0) THEN
		IFAULT=6; PRINT *, IFAULT,"�ʐς�0�ȉ�"; RETURN
	END IF
	LSUM=DLOG(SSUM(KMAX))
!	���z�m���ɕϊ�
	DO J=1,KMAX-1
		R=SSUM(J)/SSUM(KMAX)
		SSUM(J)=R
	END DO
	SSUM(KMAX)=1.0

!	�i0,1�j�̈�l�����𔭐������C�G���AJ������
	R=RANSU2()
	IF(R==0.0) THEN
		PRINT *, "����2��0��"; STOP
	END IF
	DO J=1,KMAX
		IF(SSUM(J) > R) EXIT	!�G���Aj(Zj-1 < X* < Zj)
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
	SUBROUTINE UAREA(KMAX,SSUM)	! ����Uk���̊e�ʐς����߂�
	INTEGER,INTENT(IN):: KMAX
	REAL*8,INTENT(OUT):: SSUM(KMAX)	!Uk���̍��v�ʐ�
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
*     ��핪�z�ւ̎�������iGEWEKE�̉�������j
*******************************************************
	REAL FUNCTION GEWEKEZ(SAMP,NPOINT)
	IMPLICIT NONE
	REAL*8,EXTERNAL:: TVORA
	INTEGER,INTENT(IN):: NPOINT ! �W�{��
	REAL*8,INTENT(IN):: SAMP(NPOINT)	!�W�{�W��
	INTEGER NP1, NP2
	REAL*8 SAMP1(NPOINT/10), SAMP2(NPOINT/50)
	REAL*8 MYU1,MYU2,S1,S2
	REAL*8 Z	! ���蓝�v��
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
*	�W�{���U�̌v�Z�i���n��j
! �Q�l����
! A Simple, Positive Semi-Definite, Heteroskedasticity
!	and Autocorrelation Consistent Covariance Matrix
*****************************************
	REAL*8 FUNCTION TVORA(SAMPX,MYU,NN)
	IMPLICIT NONE
	INTEGER,INTENT(IN):: NN
	REAL*8,INTENT(IN):: SAMPX(NN),MYU
	INTEGER I,J,M

		M = 10	! �o���h��
		TVORA= OMEGA(0)
		DO J=1, M
			TVORA = TVORA + (1.0 - DBLE(J)/(M+1.0))*2.0*OMEGA(J)
		END DO
	CONTAINS
!--------------------
*	�����葱��
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
*	�^����l�����̔����i��Z�������@�j
***************************************************
	REAL*8 FUNCTION RANSU1()
	INTEGER IR1
	COMMON /IRAND1/IR1
		IR1=IAND(48828125*IR1, 2147483647)
		RANSU1=DBLE(IR1)/DBLE(2147483647)
	END FUNCTION
***************************************************
*	�^����l�����̔����i�}���`��Z�������@�FAS183�j
!	����������
***************************************************
	REAL*8 FUNCTION RANSU2()
!	�����̎�͂������1����30000�͈̔͂ŗ^���Ȃ���΂����Ȃ�
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
*	�o�u���\�[�g
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
*	�t�s�� (���S�s�{�b�g�@)
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
!	(N-K+1)�����s��̒������Βl���ő�̂��̂�T��
		AW=DABS(A(K,K))
		P=K; Q=K
		DO J=K,N; DO I=K,N
			IF(AW < DABS(A(I,J))) THEN
				AW=DABS(A(I,J)); P=I; Q=J
			END IF
		END DO; END DO
		IF(AW < EPS) THEN
			PRINT *, "���̍s�����" ;STOP
		END IF
!	��Βl���ő�̂��̂�����Ɏ����Ă���
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
!	�j�s�ɑ΂��鑀��
		AW=1.0/A(K,K);	A(K,K)=1.0
		DO J=1,N; A(K,J)=A(K,J)*AW;	END DO
!	�j�s�ȊO�ւ̑���
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