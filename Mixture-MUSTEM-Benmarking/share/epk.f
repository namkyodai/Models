!***********************************************************************
!**   Individual deterioration speed evaluation（Initial value calculation of εk: Unrestricted
!                     **
!***********************************************************************
      implicit real*8(a-h,o-z)
	parameter idim0=1600,jdim0=1600
	dimension ik(idim0,jdim0),jk(idim0,jdim0),idif(idim0)
	dimension zk(idim0,jdim0),xk(idim0,jdim0,10),xk0(10)
	dimension beta(10,10),ep(idim0),ep0(idim0),ep00(idim0)
      character*20 file1,file2,diff0,cdif(idim0)
	dimension ts(30)
	common /idt/ ik,jk,zk,xk,idif
c	common /pr1/ beta,rmu,phi
	common /pr2/ ep00,ep,kno
	common /ini/ jmax,mmax,itotal
c
		!An initial parameter is operated. (BETA(J)=SJ1*J+SB1)
		REAL*8,PARAMETER:: SB1=1.1D-03,SJ1=3.0D-03	
c
	open(5,file="input.dat")
c
	jmax=5
	mmax=3
	itotal=8
	ndim=(jmax-1)*mmax
c
      do 100 i=1,itotal
c
      !Heterogeneity and characteristic vector at state of deterioration for at the 
!	time of two and check period
      	read(5,*) ik0,jk0,zk0,diff0,(xk0(m),m=1,mmax)
c
	!Making of data base of heterogeneity(heterogeneity number and data number)
c
     		if (i.ne.1) then
			iflg=0
			do 200 j=1,kno
				if (iflg.eq.1) then
					exit
      			else if (diff0.eq.cdif(j)) then
      				idif(j)=idif(j)+1
					ik(j,idif(j))=ik0
					jk(j,idif(j))=jk0
					zk(j,idif(j))=zk0
     					do 210 m=1,mmax
     						xk(j,idif(j),m)=xk0(m)
  210					continue
      				iflg=1
				end if
  200			continue
			if (iflg.eq.0) then
				kno=kno+1
				cdif(kno)=diff0
				idif(kno)=1
				ik(kno,1)=ik0
				jk(kno,1)=jk0
				zk(kno,1)=zk0
     				do 220 m=1,mmax
     					xk(kno,1,m)=xk0(m)
  220				continue
      			iflg=1
			end if
c
		else
     			cdif(1)=diff0
     			idif(1)=1
     			ik(1,1)=ik0
     			jk(1,1)=jk0
     			zk(1,1)=zk0
     			do 230 m=1,mmax
     				xk(1,1,m)=xk0(m)
  230			continue
			kno=1
      	end if
c
  100 continue
c      write(*,*) 'Reading end'
	close(5)
c
c     Step 1：Setting of initial value（beta，phi，epk）
c
	do 300 j=1,jmax-1
	do 310 m=1,mmax
c
		BETA(M,J)=J*SJ1+SB1	!ベター
!		BETA(M,J)=1.0D-06+1.1D-05*J
c
!		beta(m,j)=0.01
c		beta(m,j)=j*0.01
!		beta(m,j)=j*2.2d-02+9.0d-03
  310 continue
  300 continue
	rmu=1.5
	phi=0.4
	do 320 j=1,kno
		ep00(j)=1.0
		ep(j)=ep00(j)
  320 continue
c	Output of initial value
	write(*,*) 'Output of initial value'
	do 330 j=1,jmax-1
		write(*,*) (BETA(M,J),m=1,mmax)
  330 continue
	write(*,*) 'rmu=',rmu
	write(*,*) 'phi=',phi
	do 340 j=1,kno
		write(*,*) 'ep00(',j,')',ep00(j)
  340 continue
	
c
c     Step 2：An initial value of the heterogeneity parameter is presumed by using the pattern method. (unrestricted)

c     Step 3：Presumption of limiting attached heterogeneity parameter
c	call epsilon(ep,kno,beta,rmu,phi)
c
c     Step 4：Presumption of characteristic parameter

c
	do 700 j=1,jmax-1
	do 700 m=1,mmax
	    ts((j-1)*mmax+m)=beta(m,j)
  700 continue
	ts((jmax-1)*mmax+1)=rmu
	ts((jmax-1)*mmax+2)=phi
	call Pattern_Method(ts,ndim+2)

	stop
	end
c
!***********************************************************************
!**   Individual deterioration speed evaluation(presumption of εk)                               **
      subroutine epsilon(ep,kno,beta0,rmu0,phi0)
	implicit real*8(a-h,o-z)
	parameter idim0=1600,jdim0=1600
	dimension ik(idim0,jdim0),jk(idim0,jdim0),idif(idim0)
	dimension zk(idim0,jdim0),xk(idim0,jdim0,10),xk0(10)
	dimension beta(10,10),ep(idim0),ep0(idim0),beta0(10,10)
      character*20 file1,file2,diff0,cdif(idim0)
	dimension ts(30)
	common /idt/ ik,jk,zk,xk,idif
	common /pr1/ beta,rmu,phi
	common /ini/ jmax,mmax,itotal
c
c     Step 2：An initial value of the heterogeneity parameter is presumed by using the 
!		pattern method. (unrestricted)
c
	do 300 j=1,jmax-1
		do 310 m=1,mmax
			BETA(M,J)=BETA0(M,J)
  310		continue
  300 continue
	rmu=rmu0
	phi=phi0
c
c	write(*,*) 'Presumption of εk',beta,phi
c
	do 400 k=1,kno
		call Pattern_Methodep(ep(k),k,lth)
		ep0(k)=ep(k)
c	write(*,2000) 'kth==epk', k,'(',idif(k),')/',kno,ep(k),lth
		sumepk0=sumepk0+ep(k)
  400 continue
 2000 format(a,i5,a,i5,a,i5,e12.5,i6)
c
	if(sumepk0.le.kno) then
		ifg=1
	else
		ifg=0
	end if
c
c     Step 3：Presumption of limiting attached heterogeneity parameter

c
	rlambda=0.0
	jfg=2
	kfg=2
	setp=5.0
	do 500 i=1,100
		rlambda=rlambda+setp
		sumepkp=0.0
		sumepkm=0.0
c
		if (jfg.eq.2) then
			do 410 k=1,kno
				ep(k)=ep0(k)
  410			continue
			call optimum(sumepkp,rlambda,ep,kno)
			if (ifg.eq.0) then
				if (sumepkp.gt.sumepk0) then
					jfg=1
				else if (sumepkp.gt.kno) then
					jfg=2
				else
					jfg=3
					exit
				end if
			else
				if (sumepkp.lt.sumepk0) then
					jfg=1
				else if (sumepkp.lt.kno) then
					jfg=2
				else
					jfg=3
					exit
				end if
			end if
		end if
c
		if (kfg.eq.2) then
			do 430 k=1,kno
				ep(k)=ep0(k)
  430			continue
			call optimum(sumepkm,-rlambda,ep,kno)
			if (ifg.eq.0) then
				if (sumepkm.gt.sumepk0) then
					kfg=1
				else if (sumepkm.gt.kno) then
					kfg=2
				else
					kfg=3
					exit
				end if
			else
				if (sumepkm.lt.sumepk0) then
					kfg=1
				else if (sumepkm.lt.kno) then
					kfg=2
				else
					kfg=3
					exit
				end if
			end if
		end if
c
  500 continue
c
	if (jfg.eq.3) then
		fiba=rlambda-setp
		fibb=rlambda
	else if (kfg.eq.3) then
		fiba=-rlambda-setp
		fibb=-rlambda
	else
		write(*,*) 'The retrieval range is undecided. '
		stop
	end if
c
	k=20
	call Fibonacci(fiba,fibb,k,ep,kno)
	rlambda=(fiba+fibb)/2.0
		do 450 k=1,kno
			ep(k)=ep0(k)
  450		continue
	call optimum(sumepk,rlambda,ep,kno)
		do 460 k=1,kno
	write(*,2100) 'kth==epk', k,'(',idif(k),')/',kno,ep(k),rlambda
  460		continue
 2100 format(a,i5,a,i5,a,i5,e12.5,e15.6)
c
	end subroutine
c
c

call Pattern_Methodep(ts,dim)

!***********************************************************************
!**   It is personally a subroutine of the search (pattern method). (heterogeneity)             **
!***********************************************************************
	subroutine Pattern_Methodep(ts,kth,lth)
c
      implicit real*8(a-h,o-z)
	parameter idim0=1600,jdim0=1600
	dimension ik(idim0,jdim0),jk(idim0,jdim0),idif(idim0)
	dimension zk(idim0,jdim0),xk(idim0,jdim0,10),xk0(10)
	dimension beta(10,10),ep(idim0),ep00(idim0)
      character file1*20,file2*20,diff0*20,cdif(idim0)*20
	common /idt/ ik,jk,zk,xk,idif
	common /pr1/ beta,rmu,phi
	common /pr2/ ep00,ep,kno
	common /ini/ jmax,mmax,itotal
c
	dimension maxobj(2),t(2,1),del(2,1)
	real*8 ts,b,dels,t0,t1,t2,t3
	integer i,j,k
c
	open(1001,file='tansakue.txt')
!****Setting of width of search****
	maxobj(1)=-10000
	maxobj(2)=1
	dels=0.1d0
!*******************
!	eps=0.0000001
	eps=1.0e-10
	b=ts
c
	do 200 j=1,10000
!		write(*,*) 'epk',j
c
!		print*,OBJ(ts)
!		write(*,*) 'dels=',dels
		if (dels.gt.eps) then
c
			t(1,1)=ts
			del=0
			del(2,1)=dels
c
			t1=t(1,1)
			t2=t(1,1)+del(2,1)
			t3=t(1,1)-del(2,1)
			delsorg=dels
	        do 100 i=1,1000
	        if (t3.le.0) then
	            dels=0.5*dels
				del(2,1)=dels
				t3=t(1,1)-del(2,1)
!				write(*,*) 'Dels decrease t3',i,t3,dels
			else
!				write(*,*) 'Dels decrease EXIT',i,t3,dels
				exit
			end if
  100         continue
			call likelihoodep(op1,t1,kth)
			call likelihoodep(op2,t2,kth)
			call likelihoodep(op3,t3,kth)
c			write(*,*) 'epk-op1-op3',j,op1,t1,op2,t2,op3,t3
c			write(*,'(a,i5,6e12.3)') 
c     #               'epk-op1-op3',j,op1,t1,op2,t2,op3,t3
c
			if(op2.gt.op1.and.op2.gt.op3)then
				del(2,1)=delsorg
				dels=delsorg
				t(2,1)=t(1,1)+del(2,1)
			else if(op3.gt.op1.and.op3.gt.op2)then
				t(2,1)=t(1,1)-del(2,1)
			else if(op1.ge.op3.and.op1.ge.op2)then
				dels=delsorg
				t(2,1)=t(1,1)
			end if
c
			call likelihoodep(objept,t(2,1),kth)
			call likelihoodep(objepb,b,kth)
   			if(objept.gt.objepb)then
				ts=2.0*t(2,1)-b
				b=t(2,1)
			else
				ts=b
				dels=1.0/2.0*dels
			end if
c
		else 
			exit
		end if
c
c	write(*,*) 'ep(k)-ts(k)',kth,ep(kth),ts
c
		call likelihoodep(objepb,b,kth)
		if(OBJepb.gt.MAXOBJ(1))then
			MAXOBJ(1)=OBJepb
			MAXOBJ(2)=b
		end if
c
!		write(*,*) '<b-objepb>',b,OBJepb
!	print*,"MAXOBJ=",MAXOBJ(1)
!	print*,"optimal=",(MAXOBJ(i),i=2,dim+1)
c
		if(j.eq.50.or.j.eq.100.or.j.eq.500.or.j.eq.1000.or.j.eq.
     #             2500.or.j.eq.5000.or.j.eq.8000.or.j.eq.10000) then
		write(1001,*) MAXOBJ(1)
		write(1001,*) j,dels
		write(1001,*) MAXOBJ(2)
		end if
c
  200 continue
      lth=j
c
	write(1001,*) "MAXOBJ=",MAXOBJ(1)
	write(1001,*) j
	write(1001,*) "optimal=",MAXOBJ(2)
c
	end subroutine
c
c
**************************************************
*	Log likelihood function(heterogeneity)
**************************************************
	subroutine likelihoodep(objep,epk,kth)
!      real*8 function objep(aa,kth) 
      implicit real*8(a-h,o-z)
	parameter idim0=1600,jdim0=1600
	dimension ik(idim0,jdim0),jk(idim0,jdim0),idif(idim0)
	dimension zk(idim0,jdim0),xk(idim0,jdim0,10),xk0(10)
	dimension beta(10,10),ep(idim0),ep00(idim0)
      character file1*20,file2*20,diff0*20,cdif(idim0)*20
	dimension theta(10),thetasa(10,10),reserve(10)
	common /idt/ ik,jk,zk,xk,idif
	common /pr1/ beta,rmu,phi
	common /pr2/ ep00,ep,kno
	common /ini/ jmax,mmax,itotal
c
	ep(kth)=epk
	objep=0.0d0
c
	do 100 ith=1,idif(kth)
c
	!Calculation of Θ
		do 190 j=1,jmax
			theta(j)=0.0
			reserve(j)=0.0
  190		continue
		maxtheta=ik(kth,ith)
		do 200 j=ik(kth,ith),jk(kth,ith)
			if (j.eq.jmax) cycle
			do 210 m=1,mmax
				theta(j)=theta(j)+beta(m,j)*xk(kth,ith,m)
  210			continue
!	!To prevent the overflow, the maximum Θ is examined. 
			if (reserve(j) < theta(j)) then
				reserve(j)=theta(j)
				maxtheta=j
			end if
  200		continue	
	!Difference of each Θi-Θj

		do 240 i=1,jmax
		do 240 j=1,jmax
			thetasa(i,j)=0.0
  240		continue
  		do 250 i=ik(kth,ith),jk(kth,ith)
		do 250 j=ik(kth,ith),jk(kth,ith)
			thetasa(i,j)=theta(i)-theta(j)
!	write(*,'(a35,2i3,3f10.5)') 'thetasa(i,j)=theta(i)-theta(j)',i,j,
!     #            thetasa(i,j),theta(i),theta(j)
  250		continue
c
		sum1=0.0
		do 300 l=ik(kth,ith),jk(kth,ith)
			product=1.0
			do 310 m=ik(kth,ith),l-1
				product=product*theta(m)/thetasa(m,l)
  310			continue
			do 320 m=l,jk(kth,ith)-1
				product=product*theta(m)/thetasa(m+1,l)
  320			continue
			sum1=sum1+product*exp(-theta(l)*zk(kth,ith)*ep(kth))
  300		continue
		if (sum1.gt.0) then
			objep=objep+dlog(sum1)
		else
			sum1=1.0D-305 !The case that becomes negative also :
			! in the case where Θ is subtracted. 

			objep=objep+dlog(sum1)
		end if
  100	continue
c
c	For the gamma distribution
c	objep=objep+phi*dlog(phi)-dlog(dgamma(phi))
c     #                  +(phi-1)*dlog(ep(kth))-phi*ep(kth)
c
c	For lognormal distribution
	objep=objep-dlog(phi)-dlog(ep(kth))
     #         -((dlog(ep(kth))-dlog(rmu)+(phi**2)/2.0)**2)/(2*phi**2)	
c
c	For normal distribution (Gaussian distribution)

c	objep=objep-dlog(phi)-0.5*dlog(2*3.141592654)
c     #                  -((ep(kth)-1)**2)/(2*phi**2)	
c
c	For reverse-Gaussian distribution (Wald distribution)
c	objep=objep-0.5*dlog(phi)-0.5*dlog(2*3.141592654)
c     #           -1.5*dlog(ep(kth))-(phi*(ep(kth)-1)**2)/(2*ep(kth))	
c
	return
	end
!	end function
c
c

!******************************************************
!**			Subroutine of Fibonacci search			 **
!******************************************************

	subroutine Fibonacci(a,b,K,ep,kno)		!a,b is an edge point，
	!K is repeatedly a frequency. 
      implicit real*8(a-h,o-z)

	real(8) F(100),ep(1600),epa(1600),epb(1600)
	real(8) x1,x2,a,b,MOBJa,MOBJb
	integer i,k

!****Fibonacci sequence making****
	F(1)=1
	F(2)=1
	do 100 i=1,K-2
		F(i+2)=F(i)+F(i+1)
  100 continue
!****************************

	x1=F(K-2)/F(K)*(b-a)+a
	x2=F(K-1)/F(K)*(b-a)+a

	do 200 i=2,K-2

	do 150 kk=1,kno
		epa(kk)=ep(kk)
		epb(kk)=ep(kk)
  150 continue

		call optimum(MOBJa,a,epa,kno)
		call optimum(MOBJb,b,epb,kno)
		MOBJa=abs(MOBJa-kno)
		MOBJb=abs(MOBJb-kno)
c	write(*,*) 'MOBJa',MOBJa,a,kno
c	write(*,*) 'MOBJb',MOBJb,b,kno
		if(MOBJa.ge.MOBJb)then
			a=x1
			b=b
			x1=x2
			x2=F(K-i)/F(K+1-i)*(b-a)+a
		else
			a=a
			b=x2
			x2=x1
			x1=F(K-1-i)/F(K+1-i)*(b-a)+a
		end if
  200 continue

c	print*,a,x1,x2,b

	end subroutine

c
**************************************************
*	（Heterogeneity，λ）
**************************************************
	subroutine optimum(obj,rlambda,ep,kno)
!      real*8 function objep(aa,kth) 
      implicit real*8(a-h,o-z)
	parameter idim0=1600,jdim0=1600
	dimension ik(idim0,jdim0),jk(idim0,jdim0),idif(idim0)
	dimension zk(idim0,jdim0),xk(idim0,jdim0,10),xk0(10)
	dimension beta(10,10),ep(idim0)
      character file1*20,file2*20,diff0*20,cdif(idim0)*20
	dimension theta(10),thetasa(10,10),reserve(10)
	common /idt/ ik,jk,zk,xk,idif
	common /pr1/ beta,rmu,phi
c	common /pr2/ ep,kno
	common /ini/ jmax,mmax,itotal
c
c	Settling judgment value
	eps=1.0e-8
c
c	Agree to the lambda. 
	do 200 kth=1,kno
		dels=0.1d0	!Initial value of width of carving of εk
		do 300 i=1,10000
			if (dels.gt.eps) then
c	Decision of width of carving of step
				epk1=ep(kth)
				epk2=ep(kth)-dels
				epk3=ep(kth)+dels
				call likelihoodep2(obj1,epk1,kth)
				call likelihoodep2(obj2,epk2,kth)
				call likelihoodep2(obj3,epk3,kth)
				obj1=obj1-rlambda
				obj2=obj2-rlambda
				obj3=obj3-rlambda
c	write(*,'(2i6,8e12.4)') kth,i,dels,ep(kth),obj2,obj1,obj3,
c     #                        epk2,epk1,epk3
c
				if (abs(obj1).le.eps) then
c	write(*,*) 'Epk is an optimal value. '
					exit
				else if (obj1.ge.0.and.obj2.ge.0.and.obj3.ge.0) then
					if (obj1.ge.obj2.or.obj1.ge.obj3) then
						if (obj2.le.obj3) then
							ep(kth)=ep(kth)-dels
c	write(*,*) 'Obj2 is near 0. ，Plus side',ep(kth),dels
						else
							ep(kth)=ep(kth)+dels
c	write(*,*) 'Obj3 is near 0. ，Plus side',ep(kth),dels
						end if
					else
						dels=0.5*dels
c	write(*,*) 'Obj1 is near 0. ，Plus side',ep(kth),dels
						cycle
					end if 
				else if (obj1.le.0.and.obj2.le.0.and.obj3.le.0) then
					if (obj1.le.obj2.or.obj1.le.obj3) then
						if (obj2.ge.obj3) then
							ep(kth)=ep(kth)-dels
c	write(*,*) 'Obj2 is near 0. ，Minus side',ep(kth),dels
						else
							ep(kth)=ep(kth)+dels
c	write(*,*) 'Obj3 is near 0.，Minus side',ep(kth),dels
						end if
					else
						dels=0.5*dels
c	write(*,*) 'Obj1 is near 0. ，Minus side',ep(kth),dels
						cycle
					end if 
				else
					dels=0.5*dels
c	write(*,*) 'The sign is a failure. ',ep(kth),dels
					cycle
				end if 
			else
				exit
			end if
  300		continue
  200 continue
c
c	Total of εk and check on size of K
	obj=0.0
	do 500 kth=1,kno
		obj=obj+ep(kth)
  500 continue
c	write(*,'(i5,2f15.5,i5,f15.5)') l,kno-sumep,sumep,kno,rlambda
  100 continue
c
c	do 600 i=1,kno
c		write(*,*) 'epk(',i,')',ep(i)
c  600 continue
c
c
	return
	end
c
c
c
**************************************************
*One dimension search(heterogeneity)	
*　　The first floor partial differential subroutine  ∂lnL(θΞ:ε)/∂ε^k
**************************************************
	subroutine likelihoodep2(objep,epk,kth)
!      real*8 function objep(aa,kth) 
      implicit real*8(a-h,o-z)
	parameter idim0=1600,jdim0=1600
	dimension ik(idim0,jdim0),jk(idim0,jdim0),idif(idim0)
	dimension zk(idim0,jdim0),xk(idim0,jdim0,10),xk0(10)
	dimension beta(10,10),ep(idim0)
      character file1*20,file2*20,diff0*20,cdif(idim0)*20
	dimension theta(10),thetasa(10,10),reserve(10)
	common /idt/ ik,jk,zk,xk,idif
	common /pr1/ beta,rmu,phi
c	common /pr2/ ep,kno
	common /ini/ jmax,mmax,itotal
c
	ep(kth)=epk
	objep=0.0d0
c
	do 100 ith=1,idif(kth)
c
	!Calculation of Θ
		do 190 j=1,jmax
			theta(j)=0.0
			reserve(j)=0.0
  190		continue
		maxtheta=ik(kth,ith)
		do 200 j=ik(kth,ith),jk(kth,ith)
			if (j.eq.jmax) cycle
			do 210 m=1,mmax
				theta(j)=theta(j)+beta(m,j)*xk(kth,ith,m)
  210			continue
!	!To prevent the overflow, the maximum Θ is examined. 
			if (reserve(j) < theta(j)) then
				reserve(j)=theta(j)
				maxtheta=j
			end if
  200		continue	
	!各Θi-Θjの差
		do 240 i=1,jmax
		do 240 j=1,jmax
			thetasa(i,j)=0.0
  240		continue
  		do 250 i=ik(kth,ith),jk(kth,ith)
		do 250 j=ik(kth,ith),jk(kth,ith)
			thetasa(i,j)=theta(i)-theta(j)
  250		continue
c
		sum1=0.0
		sum2=0.0
		do 300 l=ik(kth,ith),jk(kth,ith)
			product=1.0
			do 310 m=ik(kth,ith),l-1
				product=product*theta(m)/thetasa(m,l)
  310			continue
			do 320 m=l,jk(kth,ith)-1
				product=product*theta(m)/thetasa(m+1,l)
  320			continue
			sum1=sum1+product*exp(-theta(l)*zk(kth,ith)*ep(kth))
			sum2=sum2+product*(-theta(l)*zk(kth,ith))*
     #                   exp(-theta(l)*zk(kth,ith)*ep(kth))
  300		continue
		objep=objep+(sum2/sum1)
  100	continue
c
c	For the gamma distribution
c	objep=objep+(phi-1)/ep(kth)-phi
c
c	For lognormal distribution
	objep=objep-1.0/ep(kth)
     #     -(dlog(ep(kth))-dlog(rmu)+(phi**2)/2.0)/(ep(kth)*phi**2)	
c
c	For normal distribution (Gaussian distribution)

c	objep=objep-dlog(phi)-0.5*dlog(2*3.141592654)
c     #                  -((ep(kth)-1)**2)/(2*phi**2)	
c
c	For reverse-Gaussian distribution (Wald distribution)
c	objep=objep-0.5*dlog(phi)-0.5*dlog(2*3.141592654)
c     #           -1.5*dlog(ep(kth))-(phi*(ep(kth)-1)**2)/(2*ep(kth))	
c
c
	return
	end
!	end function
c
call Pattern_Method(ts,dim)

!***********************************************************************
!**   It is personally a subroutine of the search (pattern method). (whole)                **
!***********************************************************************
	subroutine Pattern_Method(ts,dim)
c
      implicit real*8(a-h,o-z)
	parameter idim0=1600,jdim0=1600
	dimension ik(idim0,jdim0),jk(idim0,jdim0),idif(idim0)
	dimension zk(idim0,jdim0),xk(idim0,jdim0,10),xk0(10)
	dimension beta(10,10),ep(idim0),ep00(idim0)
      character file1*20,file2*20,diff0*20,cdif(idim0)*20
	integer dim
	common /idt/ ik,jk,zk,xk,idif
c	common /pr1/ beta,rmu,phi
	common /pr2/ ep00,ep,kno
	common /ini/ jmax,mmax,itotal
c
	dimension maxobj(dim+1)
	dimension ts(dim),b(dim),dels(dim),t0(dim),t1(dim),t2(dim),t3(dim)
	dimension t(dim+1,dim),del(dim+1,dim)
	integer i,j,k
c
	open(1002,file='tansaku.txt')
!****Setting of width of search****
	maxobj(1)=-10000
      do 100 i=1,dim
		maxobj(i+1)=1
		dels(i)=3.2e-03
c		dels(i)=2.8e-03
  100 continue
c 
      dels(dim-1)=0.4
      dels(dim)=0.05
c
!*******************
!	eps=0.0000001
	eps=1.0e-10
	do 110 i=1,dim
		b(i)=ts(i)
  110 continue
c	t=1.0
c
	do 200 j=1,10000
		write(*,*) 'ALL',j
c
!		print*,OBJ(ts)
		write(*,'(a,20e15.5)') 'dels',dels
		if (minval(dels).gt.eps) then
c
			do 210 i=1,dim
				t(1,i)=ts(i)
  210			continue
			del=0
			do 220 i=1,dim
				del(i+1,i)=dels(i)
  220			continue
c
			do 300 i=1,dim
				do 230 k=1,dim
c					write(*,*) 'Data replacement',t(i,k)
					t1(k)=t(i,k)
					t2(k)=t(i,k)+del(i+1,k)
					t3(k)=t(i,k)-del(i+1,k)
					if (k.eq.dim) then
c					if (k.eq.dim.and.i.eq.dim) then
						delsorg=dels(k)
				        do 235 l=1,1000
							if (t3(k).le.0) then
								dels(k)=0.5*dels(k)
								del(i+1,k)=dels(k)
								t3(k)=t(i,k)-del(i+1,k)
c								write(*,'(a,i5,4e15.5)') 
c     #			'Dels decrease t3',l,t3(k),t(i,k),del(i+1,k),dels(k)
							else
c				write(*,*) 'Dels decrease EXIT',l,t3(k),dels(k)
								exit
							end if
  235						continue
						if (k.eq.dim) then
							if (t3(k).le.0.0) then
								t3(k)=1.0D-5
!								t3(k)=1.0D-305
							end if
						end if
					end if
  230				continue

c				write(*,*) 'call likelihood(op1,t1)'
c				write(*,'(a,20e15.5)') 't1',t1
				call likelihood(op1,t1)
c				write(*,*) 'call likelihood(op2,t2)'
c				write(*,'(a,20e15.5)') 't2',t2
				call likelihood(op2,t2)
c					write(*,*) 'call likelihood(op3,t3)'
c				write(*,'(a,20e15.5)') 't3',t3
				call likelihood(op3,t3)
				write(*,'(a,2i5,3e15.7)') 'ALL',j,i,op1,op2,op3
c
				if(op2.gt.op1.and.op2.gt.op3)then
					if (i.eq.dim) then
						del(i+1,dim)=delsorg
						dels(dim)=delsorg
					end if
					do 240 k=1,dim
						t(i+1,k)=t(i,k)+del(i+1,k)
c						write(*,*) 'Data replacement 240',t(i,k),t(i+1,k)
  240					continue
				else if(op3.gt.op1.and.op3.gt.op2)then
					do 250 k=1,dim
						t(i+1,k)=t(i,k)-del(i+1,k)
						if (k.eq.dim) then
							if (t(i+1,k).le.0.0) then
								t(i+1,k)=1.0D-5
!								t(i+1,k)=1.0D-305
							end if
						end if
c					write(*,*) 'Data replacement 250',t(i,k),t(i+1,k)
  250					continue
c				else if(op1.ge.op3.and.op1.ge.op2)then
				else
					if (i.eq.dim) then
						dels(dim)=delsorg
					end if
					do 260 k=1,dim
						t(i+1,k)=t(i,k)
c						write(*,*) 'Data replacement 260',t(i,k),t(i+1,k)
  260					continue
				end if
  300			continue
c
			do 310 k=1,dim
				t0(k)=t(dim+1,k)
c				write(*,*) 'Data replacement 310',t0(k)
  310			continue
c
			call likelihood(objt,t0)
			call likelihood(objb,b)
			write(*,*) 'objt=',objt,'objb=',objb
			do 340 jj=1,kno
				write(*,*) 'ep(',jj,')',ep(jj)
  340			continue
			if(objt.gt.objb)then
				do 270 k=1,dim
					ts(k)=2.0*t(dim+1,k)-b(k)
					if (k.eq.dim) then
						if (ts(k).le.0.0) then
							ts(k)=1.0D-5
!							ts(k)=1.0D-305
						end if
					end if
					b(k)=t(dim+1,k)
				write(*,*) 'Data replacement 270',k,ts(k),b(k)
  270				continue
			else
				do 280 k=1,dim
					ts(k)=b(k)
					dels(k)=1.0/2.0*dels(k)
				write(*,*) 'Data replacement 280',k,ts(k),dels(k)
  280				continue
			end if
c
		else
			exit
		end if
c
c
		call likelihood(objb,b)
		if(OBJb.gt.MAXOBJ(1))then
			MAXOBJ(1)=OBJb
			do i=1,dim
				MAXOBJ(i+1)=b(i)
			end do
		end if
c
c		print*,b,OBJ(b)
!		print*,"MAXOBJ=",MAXOBJ(1)
!		print*,"optimal=",(MAXOBJ(i),i=2,dim+1)
c
		if(j==50.or.j==100.or.j==500.or.j==1000.or.j==2500.or.
     #						j==5000.or.j==8000.or.j==10000) then
			write(1002,*),MAXOBJ(1)
			write(1002,*),j,dels
			write(1002,*),(MAXOBJ(i),i=2,dim+1)
		end if
c
  200 continue
c
	write(1002,*)"MAXOBJ=",MAXOBJ(1)
	write(1002,*)j
	write(1002,*)"optimal=",(MAXOBJ(i),i=2,dim+1)
c
	end subroutine
c
c
**************************************************
*	Log likelihood function(whole)
**************************************************
	subroutine likelihood(obj,ts)
!	real*8 function obj(ts) 
      implicit real*8(a-h,o-z)
	parameter idim0=1600,jdim0=1600
	dimension ik(idim0,jdim0),jk(idim0,jdim0),idif(idim0)
	dimension zk(idim0,jdim0),xk(idim0,jdim0,10),xk0(10)
	dimension beta(10,10),ep(idim0),ep00(idim0)
      character file1*20,file2*20,diff0*20,cdif(idim0)*20
	dimension theta(10),thetasa(10,10),reserve(10)
	dimension ts(jmax*mmax+1)
	common /idt/ ik,jk,zk,xk,idif
c	common /pr1/ beta,rmu,phi
	common /pr2/ ep00,ep,kno
	common /ini/ jmax,mmax,itotal
c
	do 400 i=1,jmax-1
	do 400 j=1,mmax
		beta(j,i)=ts((i-1)*mmax+j)
  400 continue
      rmu=ts((jmax-1)*mmax+1)
      phi=ts((jmax-1)*mmax+2)
c
	obj=0.0
c
	do 120 j=1,kno
		ep(j)=ep00(j)
  120 continue
	call epsilon(ep,kno,beta,rmu,phi)
c
	do 50 kth=1,kno
		do 100 ith=1,idif(kth)
c
	!Calculation of Θ
			do 190 j=1,jmax
				theta(j)=0.0
				reserve(j)=0.0
  190			continue
			maxtheta=ik(kth,ith)
			do 200 j=ik(kth,ith),jk(kth,ith)
				if (j.eq.jmax) cycle
				do 210 m=1,mmax
					theta(j)=theta(j)+beta(m,j)*xk(kth,ith,m)
  210				continue
!	!To prevent the overflow, the maximum Θ is examined. 
				if (reserve(j) < theta(j)) then
					reserve(j)=theta(j)
					maxtheta=j
				end if
  200			continue	
	!Difference of each Θi-Θj
			do 240 i=1,jmax
			do 240 j=1,jmax
				thetasa(i,j)=0.0
  240			continue
  			do 250 i=ik(kth,ith),jk(kth,ith)
			do 250 j=ik(kth,ith),jk(kth,ith)
				thetasa(i,j)=theta(i)-theta(j)
  250			continue
c
			sum1=0.0
			do 300 l=ik(kth,ith),jk(kth,ith)
				product=1.0
				do 310 m=ik(kth,ith),l-1
					product=product*theta(m)/thetasa(m,l)
  310				continue
				do 320 m=l,jk(kth,ith)-1
					product=product*theta(m)/thetasa(m+1,l)
  320				continue
				sum1=sum1+product*exp(-theta(l)*zk(kth,ith)*ep(kth))
  300			continue
			if (sum1.gt.0) then
				obj=obj+dlog(sum1)
			else
				sum1=1.0D-305 !The case that becomes negative also : 
				!in the case where Θ is subtracted. 
				obj=obj+dlog(sum1)
			end if
c

  100		continue
c	write(*,*) 'The entire degree calculation<φ>',phi,dgamma(phi),ep(kth)
c	For the gamma distribution
c		obj=obj+phi*dlog(phi)-dlog(dgamma(phi))
c     #                    +(phi-1)*dlog(ep(kth))-phi*ep(kth)
c
c	For lognormal distribution
		obj=obj-dlog(phi)-dlog(ep(kth))
     #     -((dlog(ep(kth))-dlog(rmu)+(phi**2)/2.0)**2)/(2*phi**2)	
c
c	For normal distribution (Gaussian distribution)

c		obj=obj-dlog(phi)-0.5*dlog(2*3.141592654)
c     #                  -((ep(kth)-1)**2)/(2*phi**2)	
c
c	For reverse-Gaussian distribution (Wald distribution)

c		obj=obj-0.5*dlog(phi)-0.5*dlog(2*3.141592654)
c     #           -1.5*dlog(ep(kth))-(phi*(ep(kth)-1)**2)/(2*ep(kth))	
c
   50 continue
	write(*,*) 'The entire degree calculation',obj
c
	return
	end
!	end function
c
c
cc
! Gamma function in double precision
!
      function dgamma(x)
      implicit real*8 (a - h, o - z)
      parameter (
     &    p0 = 0.999999999999999990d+00, 
     &    p1 = -0.422784335098466784d+00, 
     &    p2 = -0.233093736421782878d+00, 
     &    p3 = 0.191091101387638410d+00, 
     &    p4 = -0.024552490005641278d+00, 
     &    p5 = -0.017645244547851414d+00, 
     &    p6 = 0.008023273027855346d+00)
      parameter (
     &    p7 = -0.000804329819255744d+00, 
     &    p8 = -0.000360837876648255d+00, 
     &    p9 = 0.000145596568617526d+00, 
     &    p10 = -0.000017545539395205d+00, 
     &    p11 = -0.000002591225267689d+00, 
     &    p12 = 0.000001337767384067d+00, 
     &    p13 = -0.000000199542863674d+00)
      n = nint(x - 2)
      w = x - (n + 2)
      y = ((((((((((((p13 * w + p12) * w + p11) * w + p10) * 
     &    w + p9) * w + p8) * w + p7) * w + p6) * w + p5) * 
     &    w + p4) * w + p3) * w + p2) * w + p1) * w + p0
      if (n .gt. 0) then
          w = x - 1
          do k = 2, n
              w = w * (x - k)
          end do
      else
          w = 1
          do k = 0, -n - 1
              y = y * (x + k)
          end do
      end if
      dgamma = w / y
      end