*******************************************************
*	A Poisson Hidden Markov Model
*	Version 1.0 Made by Nam Lethanh (namkyodai@gmail.com)
*	This program takes a full advantage of program Bayesian for Multi-state Markov Hazard Model
*	which was developed by Tsuda. This is for estimation of Beta and hidden variable so as to estimate the Markov transition probability
*	input of model is information related to visual inspection and also the output of potholes model obtained by running R package.
*************************************************************
	Module FILE
!	Module FILE: for number of data input/output files of the program
	! INPUT
	
	!Condition state of base course (i=1,...,M)
	character(len=30),parameter:: input_base="input/data.csv" 
	! initial value of beta, which is hazard rate of deterioration of base course
	character(len=30),parameter:: input_beta_b="input/infobeta.csv" 
	! OUTPUT


	!Base course
	character(len=30),parameter:: output_b1="output/samp_b.csv" ! 
	character(len=30),parameter:: output_b2="output/psamp_b.csv" ! 
	character(len=30),parameter:: output_b3="output/qsamp_b.csv" ! 
	character(len=30),parameter:: output_b4="output/fmi_b.csv" ! 
	character(len=30),parameter:: output_b5="output/pi_b.csv" ! transition probability
	character(len=30),parameter:: output_b6="output/result_b.csv" ! overal results
	character(len=30),parameter:: output_b7="output/beta_b.csv" ! unknown parameters
	character(len=30),parameter:: output_b8="output/av_theta_b.csv" ! average hazar rate
	character(len=30),parameter:: output_b9="output/theta_b.csv" !hazard rate
	character(len=30),parameter:: output_b10="output/sksamp_b.csv" ! new value of sk 
	character(len=30),parameter:: output_b11="output/likeli_b.csv" ! 
	
	End module
*******************************************************************************************
*------------Dimension of data and variables in the program---------------------
	Module dim

	integer, parameter:: jmax_b=5 ! condition state of base course will be jmax_b+1
	integer, parameter:: mmax_b=2 ! No. of covariates
	integer, parameter:: jmjm_b=jmax_b*mmax_b !
	integer, parameter:: ndata_b= 17000 ! 0 !1738 !8 !4 ! Obtained numbers of
	integer, parameter:: tmax=2 ! Total number of inspections
	integer, parameter:: npoint_b=10000 !00 !0 !No. of samples in simulation
	integer, parameter:: iburn_b=3000 !00 !0 !burn-in samples of base course

!-----------------------------------------------------------------------------------------
	integer, parameter:: timemax=1
	integer, parameter:: smethod=2
	integer, parameter:: qpattern=2
	logical, parameter:: preject=.true.
	integer,parameter:: infotype=3
	integer, parameter:: limit=0
	End module
***********************************************************************************
!----------------------------------------------------------------------------------
	Module Variable
		use DIM
!---------------------------------

	integer pattern
	logical kenzendo

	!Base course
	real*8 beta_b(jmax_b,mmax_b)
	data beta_b(1,1) /0.0385/
	data beta_b(2,1) /0.1/
	data beta_b(3,1) /0.055/
	data beta_b(4,1) /0.59/
	data beta_b(1,2) /0.4/
	data beta_b(2,2) /0.11/
	data beta_b(3,2) /0.01/
	data beta_b(3,2) /0.13/
	data beta_b(4,2) /0.026/
	! Numbers of parameters to be removed (If not satisfying statistical tests)
	integer,parameter:: d1_b =2
	! Numbers of parameters to be removed, for lkehood ration test)
	integer,parameter:: d2_b=2
	integer del_b1(d1_b) ! for removing unknown parameter
	integer del_b2(d2_b) ! for likelihood test

	real*8, parameter:: myubar_b=0.0
	real*8, parameter:: sigmabar_b=0.08

	data del_b1/9,10/
	data del_b2/9,10/
	real*8:: likeli_b(d2_b), sumlikeli_b(d2_b)
	real*8 myu0_b(mmax_b, jmax_b+1), roh_b(mmax_b,mmax_b,jmax_b+1)
	real*8 lamb_b(jmax_b, tmax, ndata_b)
	real*8 p_b(jmax_b+1,jmax_b+1), q_b(jmax_b+1,jmax_b+1)
	real*8 zeta_b(mmax_b, jmax_b), sigma_b(mmax_b,mmax_b, jmax_b)
	real*8 myu_b(jmax_b+1,jmax_b+1), nyu_b(jmax_b+1,jmax_b+1)
	real*8 avx_b(mmax_b)
	! Average explanatory variable of surface course
!	data avx_b(1)/0.3684/
	data avx_b(2)/0.3000/
	integer sk_b(tmax, ndata_b), uk_b(tmax, ndata_b)
	end module
!-----------------------------------
*****************************************************************
	module bound
	use dim
!--------------For base course--------------------------	
	logical ob_b(jmjm_b), ub_b(jmjm_b)
	data ob_b/jmjm_b*.False./
	data ub_b/jmjm_b*.False./
	real*8 ox_b(jmjm_b), ux_b(jmjm_b)
	real*8 htheta_b(jmjm_b)
	real*8 svs_b(jmjm_b)
	DATA SVS_b(1:4)/3.5, 3.5, 3.5, 3.5/
	DATA SVS_b(5:8)/3.5, 3.5, 3.5, 3.5/
	integer, parameter:: kmax_b=2
	integer nt_b
	real*8 sxk_b(jmax_b*mmax_b, kmax_b) ! starting point


	End module
****************************************************************
	module database
	use dim
		
	real*8 xk_b(mmax_b,mmax_b,ndata_b), zk_b(mmax_b,ndata_b)
	integer ik_b(ndata_b), jk_b(ndata_b), mk_b(mmax_b,ndata_b)!mk_b is condition states
	end module
!---------------------------------------------------
	module mtheta_b
	use dim
	real*8 theta_b(jmax_b+1), tsa_b(jmax_b+1,jmax_b+1),psum_b(jmax_b+1)
	real*8 zkk_b
	integer ikk_b, jkk_b
	integer resi_b, resj_b
	real*8 resz_b, resx_b(mmax_b,jmax_b)
	end module
	
*********************************************************************************
*********************************************************************************
						! MAIN PROGRAM STARTS HERE
*********************************************************************************
*********************************************************************************
	Program bayes_markov
	use file; use variable; use file; use mtheta_b
	use bound; use database; implicit none
	character(len=30):: hour
	integer ir1, ix, iy, iz, tane_b, tane_s, d
	common /irand1/ ir1
	common /irand2/ ix, iy, iz
	integer i,j,m,k,id_b, id_s, ii
	!! ii is for base course
	real res_b, res_s
!	real theta_b(jmax_b+1)

!---------------For base course -------------------------
	real*8 samp_b(0:npoint_b,jmjm_b) !
	real*8 psamp_b(npoint_b,jmax_b+1,jmax_b+1) !
	real*8 qsamp_b(npoint_b,jmax_b+1,jmax_b+1) !
	real*8 fmi_b(npoint_b,jmax_b+1,jmax_b+1) !
!-----------------------------------------------------------------

	
********************************************************************************
*     START FIRSTLY WITH BASE COURSE TO ESTIMATE BETA_b
!-----------~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~-------------------------
*	Initial value of MH
	SAMP_b=0.0
	DO J=1,JMAX_b; DO M=1,MMAX_b
		ID_b = (J-1)*MMAX_b+M
		SAMP_b(0,ID_b)=BETA_b(J,M)
		HTHETA_b(ID_b)=BETA_b(J,M)
	END DO; END DO
!	SAMP=4.0*SAMP/5.0
		KENZENDO=.TRUE.

	RES_b=0.0
	DO I=1,JMAX_b
	RES_b=RES_b+I
	END DO
	DO I=1,JMAX_b; DO K=1,JMAX_b
!This part is important as it defines how P and Q will be, For Poisson model, the values of P and Q must be related
		
		P_b(I,K)=(JMAX_b+1-I)/RES_b; Q_b(I,K)=(JMAX_b+1-I)/RES_b
	END DO; END DO

	MYU_b=1.0; NYU_b=1.0

!------------------------------------------------
	CALL TIME(hour); PRINT *, hour
	MYU0_b=0.0; ROH_b=0.0
*	 Random seed
	CALL CLOCK (IR1,2)
	
	IX=16251; IY=21157; IZ=3583
	TANE_b=IR1
*	Data input from
	call input_b()
	SK_b=MK_b  ! Assumption for initial of hidden vector is actually equal to

	! observed condition states
	UK_b=MK_b ! Observed condition states.
	DO 80 D=1,D2_b
	CALL MCMC_MH_b(SAMP_b,PSAMP_b,QSAMP_b,FMI_b,D)
80	END DO	
		LIKELi_b=SUMLIKeli_b/(NPOINT_b-IBURN_b)

	open(234,file=output_b11)
	do d=1, d2_b
	write(234,15) d, likeli_b(d), likeli_b(1)-likeli_b(d)
	end do
	close(234)

!	Conversion to initial information
	IF(INFOTYPE/=3) CALL REVERSE_b(2)
	Call output_b(tane_b,samp_b,psamp_b,qsamp_b,fmi_b)
	CALL TIME(hour); PRINT *, hour

!	The value of hazard function, which is estimated for base course will be used
!	as the input for calculating hazard function of surface course.


!	open(555, file=output_b8)
!	read(555,*) theta_b(i)=
!	do ii=1, jmax_b
!	print*, theta_b(2)
	stop
***************************************************************************


15	FORMAT(1X,I5,",",(F20.8,","),F12.8)

!	Conversion to initial information
	
!	END DO ! end do for thete_b to be end

1	FORMAT(1X,2(F10.4,","),F10.4)
	end program

*******************************************************************
	!	INPUT SUBROUTINE
	! Reading data from input folder!
******************************************************************
!		Input for base course
	subroutine input_b()
	use database; use variable; use file; implicit none
	character ch
	integer i,j,k,m,m1,m2, t
	real*8 maxx_b(mmax_b),xk2_b(mmax_b,ndata_b)
	!---------------------------------------------------
!	Open inspection data of base course
	open(2000,file=input_base)
	read(2000,*) ch
		do k=1, ndata_b
			read(2000,*) (mk_b(t,k),t=1,tmax),(zk_b(t,k),t=1,tmax),
     &			(xk2_b(m1,k),m1=1,mmax_b)
		
			do t=1, tmax; do m=1,mmax_b
				xk_b(m,t,k)=xk2_b(m,k)
			end do; end do
		end do
	Close(2000)
	if (infotype==3) return  ! if there is no information	
	open(2001, file=input_beta_b)
		do j=1,jmax_b
			read(2001,*) ch
			read(2001,*) ch, (myu0_b(m1,j), m1=1,mmax_b)
			do m1=1, mmax_b
				read (2001,*) ch, (roh_b(m1,m2,j), m2=1,mmax_b)
			end do
		end do
	close(2001)
	call reverse_b(1)
	return
	End subroutine

*****************************************************************************
!		OUTPUT SUBROUTINE
*****************************************************************************
	! For base course
!	subroutine output_b()

	Subroutine output_b(tane_b,samp_b,psamp_b,qsamp_b,fmi_b)
	use variable; use bound; use file; use mtheta_b; implicit none
	real, external:: prob_b, gewekez
	integer, intent(in):: tane_b
	real*8, intent(in):: samp_b(0:npoint_b,jmjm_b)
	real*8, intent(in):: psamp_b(npoint_b,jmax_b+1, jmax_b+1)
	real*8, intent(in):: qsamp_b(npoint_b,jmax_b+1, jmax_b+1)
	real*8, intent(in):: fmi_b(npoint_b,jmax_b+1, jmax_b+1)
	real*8 samp1_b(npoint_b)
	real*8 myu2_b(mmax_b,jmax_b), ss_b(mmax_b,jmax_b,jmax_b), sm_b
	real*8 ubar_b(jmax_b*mmax_b), obar_b(jmax_b*mmax_b)
	real z_b(jmax_b*mmax_b); common /geweke/z_b
	logical rej_b
	integer irej_b(jmjm_b), ntim_b(jmjm_b)
	common /acept/irej_b, ntim_b
	integer i,j,k,m1,m2,m
	real*8 pi_b(jmax_b+1, jmax_b+1), zz_b, dpi_b

	open(200,file=output_b1)
	write(200,*) "Random seed base course",",", tane_b
	write(200,*) "No. of Samples base course ",",", npoint_b
	WRITE(200,*) "All trial frequency base course"
	WRITE(200,21)  (NTIM_b(I),I=1,JMJM_b)
	WRITE(200,*) "Refusal frequency base course"
	WRITE(200,21) (IREJ_b(I),I=1,JMJM_b)
	WRITE(200,*) "state1",",,,","state2",",,,","state3",",,,","state4"
	WRITE(200,*) ("B1",",","B2",",","B3",",",J=1,4)
	DO K=1,NPOINT_b
		WRITE(200,20) (SAMP_b(K,J),J=1,JMAX_b*MMAX_b)
	END DO
	CLOSE(200)
	!psamp_b
 	OPEN(201,file=output_b2)
		WRITE(201,*) "P","M","L"
		DO K=1,NPOINT_b
		WRITE(201,20) (PSAMP_b(K,J,1),J=1,JMAX_b+1),
     &		(PSAMP_b(K,J,2),J=1,JMAX_b+1)
     & ,(PSAMP_b(K,J,3),J=1,JMAX_b+1),(PSAMP_b(K,J,4),J=1,JMAX_b+1),
     &(PSAMP_b(K,J,5),J=1,JMAX_b+1)
		END DO
	CLOSE(201)

	!Qsamp_b
	OPEN(202,file=output_b3)
		WRITE(202,*) "Q","L","I"
		DO K=1,NPOINT_b
			WRITE(202,20) (QSAMP_b(K,J,1),J=1,JMAX_b+1),
     &			(QSAMP_b(K,J,2),J=1,JMAX_b+1)
     & ,(QSAMP_b(K,J,3),J=1,JMAX_b+1),(QSAMP_b(K,J,4),J=1,JMAX_b+1),
     &(QSAMP_b(K,J,5),J=1,JMAX_b+1)
		END DO
	CLOSE(202)
	OPEN(203,file=output_b4)
		WRITE(203,*) "FMI_b"
		DO K=1,NPOINT_b
			WRITE(203,20) (FMI_b(K,J,1),J=1,JMAX_b+1),
     &			(FMI_b(K,J,2),J=1,JMAX_b+1)
     & ,(FMI_b(K,J,3),J=1,JMAX_b+1),(FMI_b(K,J,4),J=1,JMAX_b+1),
     &(FMI_b(K,J,5),J=1,JMAX_b+1)
		END DO
	CLOSE(203)

20	FORMAT(1X,30(E15.8,","))
21	FORMAT(1X,30(I7,","))
!--------------------------
*	 Output of statistic
!--------------------------
!	Sample average
	MYU2_b(1,:)= 0.0 ! 1.0
	MYU2_b(2,:)=0.0 !0.199
	DO J=1,JMAX_b; DO M1=1,MMAX_b
		DO K=1,NPOINT_b
			MYU2_b(M1,J)=MYU2_b(M1,J)+SAMP_b(K,(J-1)*MMAX_b+M1)
		END DO
	END DO;	END DO	
	MYU2_b=MYU2_b/NPOINT_b

!	Mean value calculation of É∆
	THETA_b=0.0
		DO J=1,JMAX_b
		DO M=1,MMAX_b
			THETA_b(J)=THETA_b(J) + MYU2_b(M,J)*AVX_b(M)
		END DO
!		THETA(J) = DEXP(THETA(J))
		END DO

	OPEN(205,file=output_b8)
	DO J=1,JMAX_b
	WRITE(205,20) THETA_b(J)
	END DO
	CLOSE(205)

	!PI
	ZZ_b=1.0
		DO I=1,JMAX_b+1; DO J=1,JMAX_b+1
	DPI_b=PROB_b(I,J,ZZ_b,THETA_b)
	PI_b(I,J)=DPI_b
		END DO; END DO
 	OPEN(206,file=output_b5)
		WRITE(206,*) "PI_b"
	DO I=1,JMAX_b+1
			WRITE(206,20) (PI_b(I,J),J=1,JMAX_b+1)
	END DO
	CLOSE(206)

!  Sample variance
	DO J=1, JMAX_b; DO M2=1,MMAX_b; DO M1=1,MMAX_b
		I = (J-1)*MMAX_b
		DO K=1,NPOINT_b
			SS_b(M1,M2,J)=SS_b(M1,M2,J)
     &	+(SAMP_b(K,I+M1)-MYU2_b(M1,J))*(SAMP_b(K,I+M2)-MYU2_b(M2,J))
		END DO
	END DO; END DO; END DO
	SS_b=SS_b/NPOINT_b

!	 Confidence interval
	DO J=1,JMAX_b*MMAX_b
		DO K=1,NPOINT_b; SAMP1_b(K)=SAMP_b(K,J);	
	END DO
		CALL BUBBLE(SAMP1_b,NPOINT_b) ! Rearrange the obtained data in order by 
	                               ! using BUBBLE sorting algorithm
		UBAR_b(J)=SAMP1_b(NPOINT_b*5/100)
		OBAR_b(J)=SAMP1_b(NPOINT_b*95/100)
		Z_b(J)= GEWEKEZ(SAMP1_b,NPOINT_b)
		END DO
 90	FORMAT(1X, F10.2, A12)

!	 File output
	OPEN(207,file=output_b6)
		SELECT CASE(INFOTYPE)
			CASE(1); WRITE(207,*) "Prior all information"
			CASE(2); WRITE(207,*) "Partial prior information"
			CASE(3); WRITE(207,*) "No Information"
		END SELECT
		SELECT CASE(SMETHOD)
		CASE(1); WRITE(207,*) "Sampling methods",",","Metropolis"
		CASE(2); WRITE(207,*) "Sampling methods",",","Independent MH"
		CASE(3); WRITE(207,*) "Sampling methods",",","ARMH"
		END SELECT
	
		SELECT CASE(QPATTERN)
	CASE(1); WRITE(207,*) "Prodistribution",",","Taylor"
	CASE(2); WRITE(207,*) "Prodistribution",",","Randomwalk"
	CASE(3); WRITE(207,*) "Proposeddistribution",",","Regularproper"
			CASE DEFAULT; WRITE(207,*) "Unmaking"
		END SELECT

		IF(INFOTYPE/=3) THEN
			WRITE(207,*) "Initial information"
			WRITE(207,*) 
     &		",", "state1",",,,","state2",",,,","state3",",,,","state4",",,,"
	   WRITE(207,10) "Initial average",((MYU0_b(M1,J),
     &	   M1=1,MMAX_b),J=1,JMAX_b)
	DO M1=1,MMAX_b
	WRITE(207,10) "Inidispersion",
     &	((ROH_b(M1,M2,J),M2=1,MMAX_b),J=1,JMAX_b)
	END DO
		END IF
				
		WRITE(207,*) "Posterior information"
		WRITE(207,*) 
     &",", "state1",",,,","state2",",,,","state3",",,,","state4",",,,"
		WRITE(207,10) "Sample average",((MYU2_b(M1,J),
     &		M1=1,MMAX_b),J=1,JMAX_b)
		DO M1=1,MMAX_b
	WRITE(207,10) "Sample variance",((SS_b(M1,M2,J),
     &	M2=1,MMAX_b),J=1,JMAX_b)
		END DO

	WRITE(207,10) "Confidence interval",(UBAR_b(J),J=1,JMAX_b*MMAX_b)
	! Upper level for testing (equivalent to 5%)
	WRITE(207,10) "Confidence interval",(OBAR_b(J),J=1,JMAX_b*MMAX_b)
	! Lower level for testing (Equivalent to 95% )
	WRITE(207,10) "Convergence criterion",(Z_b(J),J=1,JMAX_b*MMAX_b) 
	! Z(J) is result of Geweke test
 10	FORMAT(1X,A8,",",20(E15.8,","))
	CLOSE(207)
	OPEN(208,file=output_b7)
	do j=1, jmax_b
	WRITE(208,15) (MYU2_b(M1,J),M1=1,MMAX_b)
	END DO
	CLOSE(208)

 15	FORMAT(2(F15.8,","))


	RETURN;	END SUBROUTINE OUTPUT_b
********************************************************************


****************************  CODING PART OF  ********************************
**************************     MAIN PROGRAM ********************************
*******************************************************************************
**************************************************************
*	/ Law MH MCMC Algorithm
**************************************************************
	SUBROUTINE MCMC_MH_b(SAMP_b,PSAMP_b,QSAMP_b,FMI_b,DD_b)
	USE VARIABLE; USE BOUND; IMPLICIT NONE
	REAL,EXTERNAL:: GEWEKEZ
	REAL*8,INTENT(OUT):: SAMP_b(0:NPOINT_b,JMJM_b) !samplingsetÉ¿
	REAL*8,INTENT(OUT):: PSAMP_b(NPOINT_b,JMAX_b+1,JMAX_b+1) 
	REAL*8,INTENT(OUT):: QSAMP_b(NPOINT_b,JMAX_b+1,JMAX_b+1) 
	REAL*8,INTENT(OUT):: FMI_b(NPOINT_b,JMAX_b+1,JMAX_b+1)   
	INTEGER,INTENT(IN):: DD_b	! Likelihood ratio test for
	REAL*8 MYU2(MMAX_b,JMAX_b+1),SS(MMAX_b,MMAX_b,JMAX_b+1),SM	! Mean/SS
	REAL*8 SAMP1_b(NPOINT_b)		!SUB sampling set
	REAL Z(JMJM_b); COMMON /GEWEKE/Z
	INTEGER I,J,M,T,ID,IPOINT_b,IREV,M1,M2,K,T2,II
	LOGICAL KENTEI

!*********************
*	Burn-in
!*********************

	DO 500 T2=1,1	!
	! To converge to the mean and variance
	DO J=1,JMAX_b; DO M=1,MMAX_b
		ID = (J-1)*MMAX_b+M
		SAMP_b(0,ID)=BETA_b(J,M)
		HTHETA_b(ID)=BETA_b(J,M)
	END DO; END DO
	IPOINT_b=1
	CALL sp_b(1,SAMP_b,PSAMP_b,QSAMP_b,FMI_b,IPOINT_b,DD_b)
	DO J=1,JMAX_b
	SVS_b(JMAX_b+J)=4.0/T2
	END DO
500	END DO

	DO ID=1,JMJM_b
	 SAMP_b(0,ID) = SAMP_b(IBURN_b,ID)
	END DO
	
	

!*********************
*	Sampling
!*********************
	IPOINT_b=1; IREV=0

 100	CONTINUE
	CALL sp_b(2,SAMP_b,PSAMP_b,QSAMP_b,FMI_b,IPOINT_b,DD_b)

!*********************
*	Hypothesis testing
!*********************
	IF(DD_b==100)THEN	!----Skip the test force----------
	KENTEI=.TRUE.
	I=1
		DO J=1,JMAX_b; DO M=1,MMAX_b
			ID = (J-1)*MMAX_b+M
			DO T=1,NPOINT_b
				SAMP1_b(T)=SAMP_b(T,ID)
			END DO
!	Parameters for delete (review required)
!			IF(DEL(I)==ID) THEN
!				I=I+1; CYCLE
!			END IF
			Z(ID)= GEWEKEZ(SAMP1_b,NPOINT_b)
		
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
		IPOINT_b=NPOINT_b*3/10
		DO J=1,JMAX_b; DO M=1,MMAX_b
			ID = (J-1)*MMAX_b+M
			DO T=1,IPOINT_b
				SAMP_b(T,ID) = SAMP_b(NPOINT_b-IPOINT_b+T,ID)
			END DO
		END DO; END DO
		IPOINT_b=IPOINT_b+1
		GOTO 100
	END IF
	END IF
	RETURN; END SUBROUTINE MCMC_MH_b
***************************************************


***********************************************************************
*	MH algorithm
***************************************************
	SUBROUTINE sp_b(ITYPE_b,SAMP_b,PSAMP_b,QSAMP_b,FMI_b,IPOINT_b,DD_b)
	USE VARIABLE; USE FILE; IMPLICIT NONE
	INTEGER,INTENT(IN):: ITYPE_b
	REAL*8,INTENT(OUT):: SAMP_b(0:NPOINT_b,JMJM_b)	 !
	REAL*8,INTENT(OUT):: PSAMP_b(NPOINT_b,JMAX_b+1,JMAX_b+1) !
	REAL*8,INTENT(OUT):: QSAMP_b(NPOINT_b,JMAX_b+1,JMAX_b+1) !
	REAL*8,INTENT(OUT):: FMI_b(NPOINT_b,JMAX_b+1,JMAX_b+1)   !
	INTEGER,INTENT(IN):: DD_b	! For the likelihood ratio test
	REAL*8, EXTERNAL:: OBJ_b
	! Set sampling
	INTEGER,INTENT(IN):: IPOINT_b
	REAL*8 X	! Returns
	LOGICAL REJ
	INTEGER IREJ(JMJM_b),NTIM(JMJM_b); COMMON /ACEP/IREJ,NTIM
	CHARACTER(LEN=30):: hour
	INTEGER T,I,ID,J,M,B,TT,L,TNUM,II,C,TTT,K,D,S
	REAL*8 RESERVE_b,FMID_b,W_b
	REAL*8 AVTHETA_b(JMAX_b+1)
		INTEGER NEWU_b(TMAX,NDATA_b)!	 Subroutine forÉø
		REAL*8 NEWNYU_b(JMAX_b,JMAX_b),NEWMYU(JMAX_b,JMAX_b)
		INTEGER NEWS_b(mmax_b,NDATA_b)	!	

	IREJ=0; NTIM=0
	CALL TIME(hour); PRINT *, hour
		SELECT CASE(ITYPE_b)
			CASE(1); PRINT *, "Burn-in period"
	TNUM=IBURN_b
			CASE(2); PRINT *, "Sampling period"
	TNUM=NPOINT_b
		END SELECT

!	FMI=0.0
	DO 821 TT=IPOINT_b, TNUM

			IF(MOD(TT,100)==0.0) THEN
			PRINT*, SK_b
			END IF
			IF(MOD(TT,2)==0.0) THEN
			PRINT *, TT,"Round"
				DO J=1,JMAX_b
					PRINT*, (BETA_b(J,M),M=1,MMAX_b)
				END DO
				DO J=1,JMAX_b
					AVTHETA_b(J)=(BETA_b(J,1)+AVX_b(2)*BETA_b(J,2))
				END DO
				PRINT*,"AVG-THETA"
				PRINT*,AVTHETA_b
			END IF

			DO M=1,JMAX_b+1; DO J=1,JMAX_b+1
				PSAMP_b(TT,M,J)=P_b(M,J)
				QSAMP_b(TT,M,J)=Q_b(M,J)
			END DO; END DO
	!FMI samples (S can be derived in the process occurs)
			DO I=1,JMAX_b+1; DO M=1,JMAX_b+1
			FMID_b=0.0
			DO L=1,JMAX_b+1
			FMID_b=FMID_b+PSAMP_b(TT,M,L)*QSAMP_b(TT,L,I)
			END DO
			FMI_b(TT,M,I)=FMID_b
			END DO; END DO

	! Bat dau tinh toan gia tri Beta bang MCMC su dung Metropolis-Hasting
	CALL NORMAL_MH_b(DD_b)
!	Kiem tra gia tri cua Beta???
	Print*, DD_b

!	Stop

!Remove step parameter----------------------
		IF(D1_b/=0) THEN
		DO S=1,JMAX_b; DO M=1,MMAX_b
		DO D=1,D1_b
		ID = (S-1)*MMAX_b+M
		IF(MMAX_b*(S-1)+M==DEL_b1(D))THEN
		BETA_b(S,M)=0.0
		END IF
		IF(MMAX_b*(S-1)+M==DEL_b2(DD_b)) THEN
		BETA_b(S,M)=0.0
		END IF
				SAMP_b(TT,ID)=BETA_b(S,M)
		END DO;END DO
		END DO
		ELSE
		DO S=1,JMAX_b; DO M=1,MMAX_b
		ID = (S-1)*MMAX_b+M
		IF(MMAX_b*(S-1)+M==DEL_b2(DD_b)) THEN
		BETA_b(S,M)=0.0
		END IF
				SAMP_b(TT,ID)=BETA_b(S,M)
		END DO;END DO
		END IF
!-------------------------------------------------
			CALL GIBBS_b(NEWS_b)
			SK_b=NEWS_b
			W_b=OBJ_b(BETA_b)
			IF(TT>=IPOINT_b) SUMLIKELi_b(DD_b)=SUMLIKELi_b(DD_b)+W_b

821	END DO

	open(123, file=output_b10)
	do k=1, ndata_b
		write(123,22) (sk_b(t,k),t=1,tmax)
	end do
	close(123)
	
	
22	FORMAT(1X,5(I5,","))
	

	CONTAINS


**********************************************************************
!---------------------------------
*S_T^K  Internal procedures for S_T ^ K
!---------------------------------
	SUBROUTINE GIBBS_b(NEWS_b)
	USE VARIABLE; 	USE DATABASE ; IMPLICIT NONE
	INTEGER,INTENT(OUT):: NEWS_b(mmax_b,NDATA_b)
	REAL*8,EXTERNAL:: PROB_b
	REAL*8,EXTERNAL:: RANSU1
	REAL*8 FJM_b(JMAX_b+1,JMAX_b+1),THETA_b(JMAX_b+1)
	REAL*8 WJT_b(JMAX_b+1,mmax_b)
	REAL*8 RES,MAXP,PSTKI,W_b,A,B,C,D,E
	INTEGER I,L,M,TK,J,IDUM,TKL,K
	DO 300 K=1,NDATA_b
	RES=0.0 ; THETA_b=0.0
	DO I=1,JMAX_b; DO M=1,JMAX_b

	RES=0.0
		DO L=1,I
			RES=RES+P_b(M,L)*Q_b(L,I)
		END DO
		FJM_b(I,M)=RES
	END DO; END DO

	WJT_b=0.0
	DO TK=1,mmax_b
	THETA_b=0.0 
	DO J=1,JMAX_b
	IF(J==JMAX_b+1) THEN
	THETA_b(J)=0.0
	ELSE
			DO M = 1, MMAX_b
				THETA_b(J)=THETA_b(J) + BETA_b(J,M)*XK_b(M,TK,K)
			END DO
	END IF
	END DO
	DO J=1,JMAX_b
		theta_b(jmax_b+1)=0.0
		IF(TK==1) THEN
!!!Check out!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			WJT_b(J,TK)=PROB_b(1,J,ZK_b(1,K),THETA_b)*
     &			PROB_b(J,SK_b(2,K),ZK_b(1,K),THETA_b)

		ELSE IF(2 <= TK .AND. TK <= mmax_b-1) THEN
			WJT_b(J,TK)=PROB_b(SK_b(TK-1,K),J,ZK_b(TK,K),THETA_b)*
     & PROB_b(J,SK_b(TK+1,K),ZK_b(TK,K),THETA_b)
		ELSE
			WJT_b(J,TK)=PROB_b(SK_b(mmax_b-1,K),J,ZK_b(TK,K),THETA_b)
		END IF
	END DO

	END DO
	DO TK=1,mmax_b
	RES=0.0
! Stock room to reconsider (?)
		IF(TK==1)THEN	! TK = 1 is special
		DO J=1,SK_b(TK+1,K)
			RES=RES+WJT_b(J,TK)*FJM_b(J,MK_b(TK,K))
		END DO
		ELSE
		DO J=SK_b(TK-1,K),SK_b(TK,K)
			RES=RES+WJT_b(J,TK)*FJM_b(J,MK_b(TK,K))
		END DO
		END IF

	PSTKI=0.0
	W_b=RANSU1()

	IF(TK==1)THEN
	DO I=1,SK_b(TK+1,K)
	PSTKI=PSTKI+WJT_b(I,TK)*FJM_b(I,MK_b(TK,K))/RES
	IF(W_b < PSTKI) THEN
	TKL=I ;EXIT
	END IF
	END DO

	ELSE
	DO I=SK_b(TK-1,K),JMAX_b
	PSTKI=PSTKI+WJT_b(I,TK)*FJM_b(I,MK_b(TK,K))/RES
	IF(W_b < PSTKI) THEN
	TKL=I ;EXIT
	END IF
	END DO
	END IF

*******************Eliminate the force error************************************
	TKL=MK_b(TK,K)
*****************************************************************************

		NEWS_b(TK,K)=TKL
		SK_b(TK,K)=TKL

	END DO
300	END DO
	END SUBROUTINE GIBBS_b
!---------------------------------------------------------

	END SUBROUTINE sp_b

*************************************************************************
!	Subroutine MCMC_MH using Metro-polish hasting algorithm
*******************************************************
	SUBROUTINE NORMAL_MH_b(DD_b)
	USE VARIABLE; USE FILE; USE BOUND; USE DATABASE;IMPLICIT NONE
	REAL*8, EXTERNAL:: OBJ_b,RANSU1,RSEIKI
	INTEGER,INTENT(IN):: DD_b	! For the likelihood ratio test
	REAL*8:: BETAA_b(JMAX_b,MMAX_b)
	REAL*8:: DBETA_b(JMAX_b,MMAX_b)
	REAL*8:: BETAD_b(JMAX_b,MMAX_b)
	REAL*8:: BETAE_b(JMAX_b,MMAX_b)
	REAL*8:: POS_b(JMAX_b,MMAX_b),W_b,H_b,LR_b(JMAX_b,MMAX_b)
	REAL*8:: OBJE_b(D2_b)
	REAL*8:: RNOR_b	!regular random
	REAL*8:: R_b
	INTEGER:: S_b,TIME,M_b,D_b,T
	REAL*8,PARAMETER:: EPS=1.0D-15
	OBJE_b=0.0
	
	! Initial value of Beta has been defined in MODULE VARIABLE

	DO 88 T=1,1
!	---Remove step parameter-------------------
	IF(D1_b/=0) THEN
	DO S_b=1,JMAX_b; DO M_b=1,MMAX_b
	DO D_b=1,D1_b

	IF(MMAX_b*(S_b-1)+M_b==DEL_b1(D_b))THEN
	BETA_b(S_b,M_b)=0.0
	END IF
	IF(MMAX_b*(S_b-1)+M_b==DEL_b2(DD_b)) THEN
	BETA_b(S_b,M_b)=0.0
	END IF
	END DO;END DO
	 END DO
	ELSE
	DO S_b=1,JMAX_b; DO M_b=1,MMAX_b
	IF(MMAX_b*(S_b-1)+M_b==DEL_b2(DD_b)) THEN
	BETA_b(S_b,M_b)=0.0
	END IF
	END DO;END DO
	END IF
	
	! Initial value of Beta has been defined in MODULE VARIABLE

!-------------------------------------------------
	BETAD_b=BETA_b

	DO 22 TIME=1,TIMEMAX

	DO M_b=1,MMAX_b; DO S_b=1,JMAX_b

	IF(MMAX_b*(S_b-1)+M_b==DEL_b2(DD_b)) GO TO 55

	DO D_b=1,D1_b
	IF(MMAX_b*(S_b-1)+M_b==DEL_b1(D_b)) GO TO 55
	END DO
	R_b=RANSU1()
	RNOR_b=RSEIKI(MYUBAR_b,SIGMABAR_b)

			BETAA_b=BETAD_b
			POS_b=BETAD_b
			BETAA_b(S_b,M_b)=BETAA_b(S_b,M_b)+RNOR_b

	H_b=OBJ_b(BETAD_b)
			IF(OBJ_b(BETAA_b)-H_b > 0.0) THEN
				POS_b(S_b,M_b)=BETAA_b(S_b,M_b)
			ELSE IF(OBJ_b(BETAA_b)-H_b > DLOG(R_b)) THEN
				POS_b(S_b,M_b)=BETAA_b(S_b,M_b)
			ELSE
			END IF
			BETAD_b=POS_b

55	END DO
	END DO

			BETA_b=BETAD_b

20	FORMAT(1X,3(F12.8,","),F12.8)
	
	W_b=OBJ_b(BETAD_b)
	PRINT*,"Log likelihood",W_b

22	END DO
88	END DO
      DBETA_b=BETA_b
!	Kiem tra gia tri moi cua Beta.
!	do k=1, jmax
		Print*, "Gia tri moi cua Beta"
		Print*, DBETA_b
!	end do
!	stop
99	CONTINUE

	END SUBROUTINE NORMAL_MH_b
****************************************************

*******************************************************************************************************************
	subroutine thetavalue_b(tik,tjk,xkk_b,beta_b,thetasa_b,mintheta_b)
	use DIM; use mtheta_b; implicit none
	real*8, intent(in):: xkk_b(mmax_b,jmax_b), beta_b(jmax_b,mmax_b)
	real*8, intent(out):: thetasa_b(jmax_b+1,jmax_b+1)
	integer, intent(in):: tik, tjk
	integer, intent(out):: mintheta_b
	integer i,j,m
	real*8 reserve_b
	theta_b=0.0
		do j=tik, tjk
			if(j==jmax_b+1) then
	theta_b(j)=0.0
	else
		do m=1, mmax_b
		theta_b(j)=theta_b(j)+beta_b(j,m)*xkk_b(m,j)
		end do
	end if
	end do
	thetasa_b=0.0
	do j=tik, tjk; do i=tik, tjk
		thetasa_b(i,j)=theta_b(i)-theta_b(j)
	end do
	end do
	! minimum to prevent overflow
	mintheta_b=tik; reserve_b=dexp(600.0D0)
	do j=tik, tjk
		if (j==jmax_b+1) cycle
		if(reserve_b>theta_b(j)) then
		reserve_b=theta_b(j); mintheta_b=j
		end if
	end do
	end subroutine
**************************************************************************************************************
*	Log-likelihood function / For complete posterior distribution/ For the base course
****************************************************
	real*8 function obj_b(dbeta_b)
	use variable; use file; use bound; use database
	use mtheta_b;implicit none
	real*8, intent(in):: dbeta_b(jmax_b,mmax_b)
	real*8, external:: prob_b
	real*8 thetasa_b(jmax_b+1,jmax_b+1)
	real*8 xkk_b(mmax_b,jmax_b)
	integer i,j,l,r,k,m,mintheta_b
	integer ek,fk,gk ! Former condition states
	integer v,x,y,tk,nk
	real*8 zk1,zk2,z,w,u
	real*8 pi,qq,aa
	obj_b=0.0
	v=0
	l=0
	do 101 nk=1,ndata_b ! Loading data
	do tk=1, tmax
		!! Between the two as a weak state, inspection intervals, characteristic vector
!	 Sampling or from anywhere
	if (kenzendo==.TRUE.) then
		if (tk==1) then
			ik_b(nk)=1; jk_b=sk_b(tk,nk)
			ikk_b=1; jkk_b=sk_b(tk,nk)
		else
	ik_b(nk)=sk_b(tk-1,nk); jk_b(nk)=sk_b(tk,nk)
	!! sk condition state of each extraction
	ikk_b=sk_b(tk-1,nk); jkk_b=sk_b(tk,nk)
	end if
	else
		if (tk==1) then
			ik_b(nk)=1; jk_b(nk)=mk_b(tk,nk)
			ikk_b=1; jkk_b=mk_b(tk,nk)
		else
	ik_b(nk)=mk_b(tk-1,nk); jk_b(nk)=mk_b(tk,nk) !
	! mk is extracted from the condition states of each time
	ikk_b=mk_b(tk-1,nk); jkk_b=mk_b(tk,nk) 
	end if
	end if
!	
	do j=1, jmax_b; do m=1, mmax_b
		xkk_b(m,j)=xk_b(m,tk,nk)
	end do; end do
		do i=1,d1_b ! removed from parameter space
			m=mod(del_b1(i)+(mmax_b-1), mmax_b)+1
			j=(del_b1(i)-m)/mmax_b+1
		xkk_b(m,j)=0.0
		end do
	! É∆ values
	call thetavalue_b(ikk_b,jkk_b,xkk_b,dbeta_b,thetasa_b,mintheta_b)
		pi=prob_b(ikk_b,jkk_b,zk_b(tk,nk),theta_b)
		!! Likelihood function
		if(pi<=1.0D-305) then
	obj_b=obj_b-702.0
	else if(pi>1.0000) then
	obj_b=obj_b
	else
		qq = dlog(pi)
		obj_b=obj_b+qq
	end if
	end do
101	end do
	close(1)
	close(1009)
20	format(1x,5(i5,","),f12.8)
	return
	end function
*******************************************
************************************
*	Element ÉŒ_ij(Z) of transition probability matrix
!Function to caculate the fraction 
*	\prod_{e=a,\ne c}^{b} \frac{\theta_k}{\bar{theta}_{e}-\bar{theta}_{c}}
************************************

	REAL*8 FUNCTION PROB_b(I,J,Z,THETA_b)
	USE VARIABLE;IMPLICIT NONE
	INTEGER,INTENT(IN):: I,J
	REAL*8,INTENT(IN):: Z,THETA_b(JMAX_b+1)
	REAL*8,EXTERNAL:: PROD2_b	
	!For the time being in this (later union necessary)
	REAL*8 RESERVE_b
	INTEGER K
	PROB_b=0.0
	IF(Z <= 1.0D-06) RETURN
	RESERVE_b=1.0
	DO K=I,J-1
		RESERVE_b=RESERVE_b*THETA_b(K)
	END DO
	DO K=I,J
		PROB_b=PROB_b+DEXP(-THETA_b(K)*Z)/PROD2_b(I,J,K,THETA_b,JMAX_b)
	END DO
	PROB_b=RESERVE_b*PROB_b
	RETURN; END
!------------------------------------------------------
* Function to caculate the denominator in the fraction 
*	\prod_{e=a,\ne c}^{b} \frac{1}{\bar{theta}_{e}-\bar{theta}_{c}}
!------------------------------------------------------
	REAL*8 FUNCTION PROD2_b(A,B,C,THETA_b,JMAX_b)
	IMPLICIT NONE
	INTEGER,INTENT(IN):: A,B,C,JMAX_b
	REAL*8,INTENT(IN):: THETA_b(JMAX_b+1)
	INTEGER E

	PROD2_b = 1.0
	DO E = A, B
		IF(E /= C) PROD2_b = (THETA_b(E)-THETA_b(C)) * PROD2_b
	END DO
	RETURN;	END
***************************************************
*	Scale-adjustment
**********************************************
	SUBROUTINE REVERSE_b(IP)
	USE VARIABLE; IMPLICIT NONE
	INTEGER,INTENT(IN):: IP
	REAL*8 A_b(MMAX_b,MMAX_b)
	REAL*8 MAXX_b(MMAX_b)	!Maxium numbers of covariates 
	!The maximum value of explanatory variable
	COMMON /ARR/MAXX_b
	INTEGER J,M1,M2
	
	IF(IP==1) THEN
		DO J=1,JMAX_b
		DO M1=2,MMAX_b; MYU0_b(M1,J)=MYU0_b(M1,J)*MAXX_b(M1); END DO
			DO M1=1,MMAX_b; DO M2=1,MMAX_b
				A_b(M2,M1) = ROH_b(M2,M1,J)*MAXX_b(M2)*MAXX_b(M1)
			END DO; END DO
			CALL MATRIX_INVERSE(A_b,MMAX_b)
			DO M1=1,MMAX_b; DO M2=1,MMAX_b
				ROH_b(M2,M1,J) = A_b(M2,M1)
			END DO; END DO
		END DO
	ELSE	! ! Reverse-adjustment
		DO J=1,JMAX_b
			DO M1=2,MMAX_b; MYU0_b(M1,J)=MYU0_b(M1,J)/MAXX_b(M1); END DO
			DO M1=1,MMAX_b; DO M2=1,MMAX_b
				A_b(M2,M1) = ROH_b(M2,M1,J)
			END DO; END DO
			CALL MATRIX_INVERSE(A_b,MMAX_b)
			DO M1=1,MMAX_b; DO M2=1,MMAX_b
				ROH_b(M2,M1,J) = A_b(M2,M1)/MAXX_b(M2)/MAXX_b(M1)
			END DO; END DO
		END DO
	END IF
	
	RETURN; END SUBROUTINE
*******************************************************

*******************************************************



*************************************************************************
************   THIS BELOW SUBROUTINES AND FUNCTIONS *********************
************     are common subroutine and functions **********************
*****************  which can be referred with external sources  *********
*************************************************************************

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
	Print*, rseiki
	!stop
	RETURN;	END FUNCTION
**********************************************

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
**************************************************
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
* Refer to code of FortranÅFAS183-Generation of uniform pseudo-random number
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
