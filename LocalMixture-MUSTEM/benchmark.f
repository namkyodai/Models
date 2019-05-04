!***********************************************************************
!**   benchmark                                                       **
!**   Individual deterioration speed evaluation（Calculation of εk: Unrestricted）**
!**   Lk is not multiplied.                                            **
!**   Deletion of parameter in which trouble is found                  **
!***********************************************************************
      implicit real*8(a-h,o-z)
	parameter idim0=1600,jdim0=1600
	dimension ik(idim0,jdim0),jk(idim0,jdim0),idif(idim0)
	dimension zk(idim0,jdim0),xk(idim0,jdim0,10),xk0(10)
	dimension beta(10,10),ep(idim0),ep0(idim0),ep00(idim0)
      character*20 file1,file2,diff0,cdif(idim0)!diff0 is noheterogeneity Number
	dimension ts(30)
	integer timemax   !It is repeatedly a frequency.
	common /idt/ ik,jk,zk,xk,idif,kno
	common /pr2/ ep00,ep
	common /ini/ jmax,mmax,itotal
	common /tmx/ timemax
!	parameter file1="input.dat"
!	parameter file2="output.dat"
!	parameter (b1=9.0d-03)
!	parameter (j1=2.2d-02)
c
		!An initial parameter is operated. (BETA(J)=SJ1*J+SB1)
		REAL*8,PARAMETER:: SB1=1.1D-03,SJ1=3.0D-03	
c
      jmax=5 !Number of ratings
      mmax=3 !Number of characteristic vectors
      itotal=1237 !Data number
      timemax=200 !It is a calculation repeatedly frequency.
c
	open(10,file="input.dat")
c
	ndim=(jmax-1)*mmax
c
      do 100 i=1,itotal
c
      !Heterogeneity and characteristic vector at state of deterioration for at the 
!time of two and check period
      	read(10,*) ik0,jk0,zk0,diff0,(xk0(m),m=1,mmax)
c
	!Making of data base of heterogeneity（Heterogeneity number and data number
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
  100 continue
**********************************************************
c      write(*,*) 'Reading end'
	close(10)
	open(11,file="ini_beta.txt")
	open(12,file="out.txt")
	do 300 j=1,jmax-1
	do 310 m=1,mmax
		read(11,*) beta(m,j)
  310 continue
  300 continue
	read(11,*) eta
c	Output of initial value
	write(*,*) 'Output of initial value'
	write(12,*) 'Output of initial value'
	do 330 j=1,jmax-1
		write(*,*) (BETA(M,J),m=1,mmax)
		write(12,*) (BETA(M,J),m=1,mmax)
  330 continue
	write(*,*) 'eta=',eta
	write(12,*) 'eta=',eta
c     Step 1: Presumption of characteristic parameter（β，φ）
	do 700 j=1,jmax-1
	do 700 m=1,mmax
	    ts((j-1)*mmax+m)=beta(m,j)
  700 continue
	ts((jmax-1)*mmax+1)=eta
	call Pattern_Method(ts,ndim+1)
	!stop
c     Step 2: The heterogeneity parameter is presumed ..the use of the pattern method... (unrestricted)
	do 320 j=1,kno
		ep00(j)=1.0
		ep(j)=ep00(j)
  320 continue
c	Output of initial value
	do 340 j=1,kno
		write(*,*) 'ep00(',j,')',ep00(j)
		write(12,*) 'ep00(',j,')',ep00(j)
  340 continue
	do 400 k=1,kno
		call Pattern_Methodep(ep(k),k)
		ep0(k)=ep(k)
	write(*,2000) 'kth==epk', k,'(',idif(k),')/',kno,ep(k)
	write(12,2000) 'kth==epk', k,'(',idif(k),')/',kno,ep(k)
	open(110,file="epsilon.csv")
	write(110,*) ep(k)

  400 continue
 2000 format(a,i5,a,i5,a,i5,e12.5,i6)
	close(11)
	close(12)
	stop
	end
!***********************************************************************
!**   It is personally a subroutine of the search (pattern method). (whole)**
!***********************************************************************
	subroutine Pattern_Method(ts,dim)
c
      implicit real*8(a-h,o-z)
	parameter idim0=1600,jdim0=1600
	dimension ik(idim0,jdim0),jk(idim0,jdim0),idif(idim0)
	dimension zk(idim0,jdim0),xk(idim0,jdim0,10),xk0(10)
	dimension beta(10,10)
      character file1*20,file2*20,diff0*20,cdif(idim0)*20
	integer dim
	integer timemax   !It is repeatedly a frequency. 
	common /idt/ ik,jk,zk,xk,idif,kno
	common /ini/ jmax,mmax,itotal
	common /tmx/ timemax
c
	dimension maxobj(dim+1)
	dimension ts(dim),b(dim),dels(dim),t0(dim),t1(dim),t2(dim),t3(dim)
	dimension t(dim+1,dim),del(dim+1,dim)
	integer i,j,k
c
	open(1002,file='tansaku3.txt')
!****Setting of width of search****
	maxobj(1)=-10000
      do 100 i=1,dim
		maxobj(i+1)=1
		read(11,*) dels(i)
  100 continue
c
!*****Settling judgment value*****
	eps=1.0e-10
c
	do 110 i=1,dim
		b(i)=ts(i)
  110 continue
c
	do 200 j=1,timemax
		write(*,*) 'ALL',j
		write(12,*) 'ALL',j
c
		write(*,'(a,20e15.5)') 'dels',dels
		write(12,'(a,20e15.5)') 'dels',dels
		if (minval(dels,mask=dels.gt.0.0).gt.eps) then
c
			do 210 i=1,dim
				t(1,i)=ts(i)
  210			continue
			del=0.0
			do 220 i=1,dim
				del(i+1,i)=dels(i)
  220			continue
c
			do 300 i=1,dim
				do 230 k=1,dim
					t1(k)=t(i,k)
					t2(k)=t(i,k)+del(i+1,k)
					t3(k)=t(i,k)-del(i+1,k)
					delsorg=dels(k)
					if (k.eq.dim) then
				        do 235 l=1,1000
							if (t3(k).le.0) then
								dels(k)=0.5*dels(k)
								del(i+1,k)=dels(k)
								t3(k)=t(i,k)-del(i+1,k)
							else
								exit
							end if
  235						continue
						if (t3(k).le.0.0) then
							t3(k)=1.0D-5
							del(i+1,k)=t(i,k)-t3(k)
							dels(k)=del(i+1,k)
						end if
					end if
  230				continue

				call likelihood(op1,t1)
				call likelihood(op2,t2)
				call likelihood(op3,t3)
				write(*,'(a,2i5,3e20.12,3e18.10)')
     #                      'ALL',j,i,op1,op2,op3,t1(i),t2(i),t3(i)
				write(12,'(a,2i5,3e20.12,3e18.10)')
     #                      'ALL',j,i,op1,op2,op3,t1(i),t2(i),t3(i)
c
				if(op2.gt.op1.and.op2.gt.op3)then
					if (i.eq.dim) then
						del(i+1,dim)=delsorg
						dels(dim)=delsorg
					end if
					do 240 k=1,dim
						t(i+1,k)=t(i,k)+del(i+1,k)
  240					continue
				else if(op3.gt.op1.and.op3.gt.op2)then
					do 250 k=1,dim
						t(i+1,k)=t(i,k)-del(i+1,k)
						if (k.eq.dim) then
							if (t(i+1,k).le.0.0) then
								t(i+1,k)=1.0D-5
							end if
						end if
  250					continue
				else
					if (i.eq.dim) then
						dels(dim)=delsorg
					end if
					do 260 k=1,dim
						t(i+1,k)=t(i,k)
  260					continue
				end if
  300			continue
c
			do 310 k=1,dim
				t0(k)=t(dim+1,k)
  310			continue
c
			call likelihood(objt,t0)
			call likelihood(objb,b)
			write(*,*) 'objt=',objt,'objb=',objb
			write(12,*) 'objt=',objt,'objb=',objb
c
			if(objt.gt.objb)then
				do 270 k=1,dim
					ts(k)=2.0*t(dim+1,k)-b(k)
					if (k.eq.dim) then
						if (ts(k).le.0.0) then
							ts(k)=1.0D-5
						end if
					end if
					b(k)=t(dim+1,k)
				write(*,*) 'Data replacement 270',k,ts(k),b(k)
				write(12,*) 'Data replacement 270',k,ts(k),b(k)
  270				continue
			else
				do 280 k=1,dim
					ts(k)=b(k)
					dels(k)=1.0/2.0*dels(k)
				write(*,*) 'Data replacement 280',k,ts(k),dels(k)
				write(12,*) 'Data replacement 280',k,ts(k),dels(k)
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
      implicit real*8(a-h,o-z)
	parameter idim0=1600,jdim0=1600
	dimension ik(idim0,jdim0),jk(idim0,jdim0),idif(idim0)
	dimension zk(idim0,jdim0),xk(idim0,jdim0,10),xk0(10)
	dimension beta(10,10)
      character file1*20,file2*20,diff0*20,cdif(idim0)*20
	dimension theta(10),thetasa(10,10),reserve(10)
	dimension ts(jmax*mmax+1)
	common /idt/ ik,jk,zk,xk,idif,kno
	common /pr1/ beta,eta, phi
	common /ini/ jmax,mmax,itotal
c
	open(1004,file='exp2.txt')
c
	do 400 i=1,jmax-1
	do 400 j=1,mmax
		beta(j,i)=ts((i-1)*mmax+j)
  400 continue
c
      eta=ts((jmax-1)*mmax+1)
	obj=0.0
c
	do 50 kth=1,kno
		sum2=0.0
		do 100 ith=1,idif(kth)
c
c	!Calculation of Θ
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
				theta(j)=exp(theta(j))
c!	!To prevent the overflow, the maximum Θ is examined. 
				if (reserve(j) < theta(j)) then
					reserve(j)=theta(j)
					maxtheta=j
				end if
  200			continue	
c	!Difference of each Θi-Θj
			do 240 i=1,jmax
			do 240 j=1,jmax
				thetasa(i,j)=0.0
  240			continue
  			do 250 i=ik(kth,ith),jk(kth,ith)
			do 250 j=ik(kth,ith),jk(kth,ith)
				thetasa(i,j)=theta(i)-theta(j)
  250			continue
			sum1=0.0
			do 300 l=ik(kth,ith),jk(kth,ith)
				product=1.0
				do 310 m=ik(kth,ith),l-1
					product=product*theta(m)/thetasa(m,l)
  310				continue
				do 320 m=l,jk(kth,ith)-1
					product=product*theta(m)/thetasa(m+1,l)
  320				continue
		
				sum1=sum1+product
     #     *exp(-theta(l)*zk(kth,ith))
     #*(1+((1/phi)*theta(l)*zk(kth,ith)**2)*0.5)
		
	

c
  300			continue
			if (sum1.gt.0) then
				sum2=sum2+dlog(sum1)
			else
				sum1=1.0D-305 !The case that becomes negative also : 
!				in the case where Θ is subtracted. 
				sum2=sum2+dlog(sum1)
			end if
  100		continue
		obj=obj+sum2	
   50 continue
	return
	end
c
!***********************************************************************
!**   It is personally a subroutine of the search (pattern method). (heterogeneity)  **
!***********************************************************************
	subroutine Pattern_Methodep(ts,kth)
c
      implicit real*8(a-h,o-z)
	parameter idim0=1600,jdim0=1600
	dimension ik(idim0,jdim0),jk(idim0,jdim0),idif(idim0)
	dimension zk(idim0,jdim0),xk(idim0,jdim0,10),xk0(10)
	dimension beta(10,10)
	integer timemax   !It is repeatedly a frequency. 
      character file1*20,file2*20,diff0*20,cdif(idim0)*20
	common /idt/ ik,jk,zk,xk,idif,kno
	common /pr1/ beta,eta, phi
	common /ini/ jmax,mmax,itotal
	common /tmx/ timemax
c
	dimension maxobj(2),t(2,1),del(2,1)
	real*8 ts,b,dels,t0,t1,t2,t3
	integer i,j,k
c
	open(1001,file='tansakue3.txt')
!****Setting of width of search****
	maxobj(1)=-10000
	maxobj(2)=1
	dels=0.1d0
!*******************
	eps=1.0e-10
	b=ts
c
	do 200 j=1,timemax
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
			else
				exit
			end if
  100         continue
			call likelihoodep(op1,t1,kth)
			call likelihoodep(op2,t2,kth)
			call likelihoodep(op3,t3,kth)
c
			if(op2.gt.op1.and.op2.gt.op3)then
				del(2,1)=delsorg
				dels=delsorg
				t(2,1)=t(1,1)+del(2,1)
			else if(op3.gt.op1.and.op3.gt.op2)then
				t(2,1)=t(1,1)-del(2,1)
			else
				dels=delsorg
				t(2,1)=t(1,1)
			end if
c
			call likelihoodep(objept,t(2,1),kth)
			call likelihoodep(objepb,b,kth)
c
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
		call likelihoodep(objepb,b,kth)
		if(OBJepb.gt.MAXOBJ(1))then
			MAXOBJ(1)=OBJepb
			MAXOBJ(2)=b
		end if
c
		if(j.eq.50.or.j.eq.100.or.j.eq.500.or.j.eq.1000.or.j.eq.
     #             2500.or.j.eq.5000.or.j.eq.8000.or.j.eq.10000) then
		write(1001,*) MAXOBJ(1)
		write(1001,*) j,dels
		write(1001,*) MAXOBJ(2)
		end if
c
  200 continue
c
	write(1001,*) "MAXOBJ=",MAXOBJ(1)
	write(1001,*) j
	write(1001,*) "optimal=",MAXOBJ(2)
c
	end subroutine
c
**************************************************
*	Log likelihood function(heterogeneity)
**************************************************
	subroutine likelihoodep(objep,epk,kth)
      implicit real*8(a-h,o-z)
	parameter idim0=1600,jdim0=1600
	dimension ik(idim0,jdim0),jk(idim0,jdim0),idif(idim0)
	dimension zk(idim0,jdim0),xk(idim0,jdim0,10),xk0(10)
	dimension beta(10,10),ep(idim0),ep00(idim0)
      character file1*20,file2*20,diff0*20,cdif(idim0)*20
	dimension theta(10),thetasa(10,10),reserve(10)
	common /idt/ ik,jk,zk,xk,idif,kno
	common /pr1/ beta,eta, phi
	common /pr2/ ep00,ep
	common /ini/ jmax,mmax,itotal
c
	open(1003,file='exp1.txt')
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
				theta(j)=exp(theta(j))
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
c
		if (sum1.gt.0) then
			objep=objep+dlog(sum1)
		else
			sum1=1.0D-305 !The case that becomes negative also : 
!			in the case where Θ is subtracted. 
			objep=objep+dlog(sum1)
		end if
  100	continue
c
	objep=objep+((phi-1)*dlog(ep(kth))-phi*ep(kth))
	

c
	return
	end
***********************************************************************
