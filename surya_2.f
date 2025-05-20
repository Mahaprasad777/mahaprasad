C                  The solar dynamo code SURYA
C
C
C   INTRODUCTION :
C     This is a 2-dimensional code for solving the kinematic dynamo
C     equation in the solar convection zone with a meridional
C     circulation.  The potential user of this code should first
C     read the accompanying 'Guide', which explains the code and
C     describes how it can be run. This code is divided in several
C     parts for the benifit of the users. These different parts
C     are explained in the 'Guide'.  Somebody familiar with the
C     solar dynamo problem should be able to run this code after
C     reading the 'Guide'.  This 'Guide' divides the potential
C     users in three levels and provides appropriate instructions
c     for users at these three different levels.
C              
C     This code was developed over the years at Indian Institute
C     of Science, Bangalore, by Arnab Rai Choudhuri and his
C     successive PhD students - Mausumi Dikpati, Dibyendu Nandy,
C     Piyali Chatterjee.
C
        implicit real*8(a-h,o-z)
        parameter (nmax=129,lmax=96)
        parameter (n_c=10230)
        common u(nmax,nmax),t,it
        common/lmx/l
        common /plsincos/pl(nmax,lmax),sn(nmax)
        double precision a(nmax,nmax),b(nmax,nmax),c(nmax,nmax),
     & d(nmax,nmax),fun(nmax),uint(nmax),
     &  e(nmax,nmax),f(nmax,nmax),vp(2*nmax,2*nmax),vq(2*nmax,2*nmax)
     &  ,a1(nmax),b1(nmax),c1(nmax),r(nmax),phi(nmax,nmax),
     &  phib(nmax,nmax),oldu(nmax,nmax),ss1(nmax),uub(nmax),
     &  uu(nmax),al(nmax,nmax),dom(nmax,nmax),ub(nmax,nmax),ra(nmax)
        double precision ab(nmax,nmax),bb(nmax,nmax),cb(nmax,nmax),
     &  db(nmax,nmax),eb(nmax,nmax),fb(nmax,nmax),eta(nmax,nmax),
     &  dror(nmax,nmax),drot(nmax,nmax),vp1(2*nmax,2*nmax),
     &  vq1(2*nmax,2*nmax),psi(2*nmax,2*nmax),etab(nmax,nmax),
     &  ss1_p(nmax),etab2(2*nmax,2*nmax),dvp(2*nmax,2*nmax),
     &  vpb(2*nmax,2*nmax),deta(2*nmax,2*nmax),uboldn,ubolds,
     &  uin(nmax,nmax),ubave(nmax),gamm_r(n_c),qs(nmax)
        external ss,erf
        
        n=nmax-1
	  pi=4.0d0*atan(1.0d0)
C
C----------------------------------------------------------------------
C	
C	PART I. This is the ONLY part of the code in which Level I users
C	may want to make some changes.  Read Sect. 3 of the 'Guide' to
C	learn about the variables specified in this part.  
C
        irelax=1
        !tmax=0.0259d0
! correlation time = 1 months
        tmax=0.726d0*9.0769d0/20.3897d0
! correlation time is different. period: 9.0769 yr
! 100% fluctuation in meri and 200% fluc in alpha

        t=0.0d0

        !s1=0.07d0
        ett=0.05d0/5.0d0
        ets=2.0d0/5.0d0
        dt=0.0002d0
C
C----------------------------------------------------------------------
C
C	PART II.  This is the part where the alpha coefficient, 
C	diffusivity, differential rotation and meridional circulation 
C	are specified. Level II Users wishing to make changes in 
C	this part should read Sect. 4 of the 'Guide'.
C       
        pm=6.96d0
	pb=0.55d0*pm
        pw=2.5d0*pm
        qm=pi
        dp=(pm-pb)/float(n)
        dq=-qm/float(n)
        ita=int((0.7d0-0.55d0)*pm/dp)+1

        n1=int((0.7d0-0.025d0-0.55d0)*pm/dp)+2
        n12=int((0.7d0+0.025d0-0.55d0)*pm/dp)+1
!        nr=int((0.71d0-0.55d0)*pm/dp)
!        mthn=int(((180.0d0-65.0d0)*pi/180.0d0)/abs(dq))
!        mths=int(((180.0d0-115.0d0)*pi/180.0d0)/abs(dq))
        open(unit=11,file='rand_al')

        do  kk=1,n_c
          read(11,*) gamm_r(kk)
        enddo

        do kk=1,n_c
        s1=0.07*gamm_r(kk)
        !s1=0.046875*gamm_r(kk)
C	
C	Here begin the do loops to calculate profiles of alpha, 
C	diffusivity and differential rotation at all the grid 
C	points. The differential rotation is written in the file 
C	'diffrot.dat'.
C
!        open(27,file='diffrot.dat')
!        open(29,file='eta.dat')
!        open(26,file='alpha.dat')   
	do i = 1, nmax
	do j = 1, nmax
	    ra(i)=pb+float(i-1)/float(n)*(pm-pb) 
            p=pb+float(i-1)/float(n)*(pm-pb)
            q=qm-float(j-1)/float(n)*qm
	    co = dcos(q)
C ALPHA PROFILE
           if(q.le.pi/2.0d0) then
          al(i,j)=s1*0.25d0*(1.0+erf((p/pm-0.95d0)/0.05d0))
     &  *(1.0d0-erf((p/pm-1.0d0)/0.01d0))*dsin(q)*co
     &  /(1.0d0+dexp(30.0d0*(pi/4.0-q)))
           else
          al(i,j)=s1*0.25d0*(1.0+erf((p/pm-0.95d0)/0.05d0))
     &  *(1.0d0-erf((p/pm-1.0d0)/0.01d0))*dsin(q)*co
     &  /(1.0d0+dexp(30.0d0*(q-3.0d0*pi/4.0d0)))
           endif
!           if(q.le.pi/2.0d0) then
!          al(i,j) =0.59230444657d0*0.2d0*dsin(6.0d0*(q-pi/2.0))
!     &  *dexp(-10.0d0*(q-pi/4.0)**2)*
!     &  0.25d0*(1.0d0+erf((p/pm-0.705d0)/0.01d0))*
!     &  (1.0d0-erf((p/pm-0.725d0)/0.01d0))
!     &  /(dexp(20.0d0*(q-pi/3.0d0))+1.0d0)
!           else
!          al(i,j) =0.59230444657d0*0.2d0*dsin(6.0d0*(q-pi/2.0))
!     &  *dexp(-10.0d0*(3.0d0*pi/4.0-q)**2)*
!     &  0.25d0*(1.0d0+erf((p/pm-0.705d0)/0.01d0))*
!     &  (1.0d0-erf((p/pm-0.725d0)/0.01d0))
!     &  /(dexp(20.0d0*(2.0d0*pi/3.0d0-q))+1.0d0)
!           endif

C DIFFUSIVITY PROFILES
          eta(i,j) = 0.00050d0 + ((ett/2.0d0)*
     &    (1.0d0+erf((p-0.7d0*pm)/(0.02d0*pm))))+((ets/2.0d0)*
     &       (1.0d0 +erf((p-0.90d0*pm)/(0.02d0*pm))))
          etab(i,j) = eta(i,j)
C SOLAR DIFFERENTIAL ROTATION PROFILE 
! here the overall dom is devided by 10 to match unit
          dom(i,j) =0.20d0*pi*(432.8d0 + 0.50d0*(1.0d0 +
     &       erf(2.0d0*(p-0.7d0*pm)/(0.05d0*pm)))*(460.7d0-
     &       62.69d0*co*co - 67.13d0*co*co*co*co-432.8d0))         
!        write(27,*)q,p,dom(i,j)
!        write(26,*) p/pm, q, al(i,j)
        end do
!         write(26,27) p, al(i,193)
  27      format(2(f13.5,1x))
!         write(29,27) p, eta(i,1)
        end do
!        close(26)
!        close(27)
!        close(29)
C 
C	The do loops are closed.  We also need the diffusivity at
C	mid-points in the grid for some calculations.  This is 
C	obtained and stored now.
C
	do i=1,2*n+1
	do j=1,2*n+1
	    p=pb+float(i-1)*(pm-pb)/float(2*n)
        etab2(i,j) = 0.00050d0 + ((ett/2.0d0)*
     &    (1.0d0+erf((p-0.7d0*pm)/(0.02d0*pm))))+((ets/2.0d0)*
     &       (1.0d0 +erf((p-0.90d0*pm)/(0.02d0*pm))))
	end do
	end do
C	
C	Derivatives of differential rotation are obtained and stored
C	now for future use.
C
      do i=2,n
	do j=2,n
	    dror(i,j)=(dom(i+1,j)-dom(i-1,j))/(2.0d0*dp)  
	    drot(i,j)=(dom(i,j+1)-dom(i,j-1))/(2.0d0*dq)
      end do
      end do
C
C MERIDIONAL CIRCULATION. Note that this is calculated both at the 
C   grid-points and mid-points. First the stream function psi(i,j)is 
C   calculated. It is written in the file 'psi.dat'.	
C
C	Now components of velocity are calculated at grid-points and 
C	mid-points from the stream function psi(i,j).
C
        v0= 10.0d0
!        open(87,file='v_theta.dat')
        p0=0.620d0*pm
        pp=p0
        gi0=pm/p0-1.0d0
        cm=0.5d0 !here I define cm for m
        cp=0.25d0 !here I define cp for p
        cc1=(2.0d0*cm+1.0d0)*(cm+cp)/((cm+1.0d0)*cp)
        cc2=(2.0d0*cm+cp+1.0d0)*cm/((cm+1.0d0)*cp)

        do i=2,2*n
        do j=2,2*n
           p=pb+float(i-1)*(pm-pb)/float(2*n)
           q=qm-float(j-1)*qm/float(2*n)
           gi=pm/p-1.0d0
        vp1(i,j)=v0*(10.0d0/4.1359d0)*((pm/p)**2)*(-1.0d0/(cm+1.0d0)+
     &  (cc1*(gi/gi0)**cm/(2.0d0*cm+1.0d0))-
     &  (cc2*(gi/gi0)**(cm+cp)/(2.0d0*cm+cp+1.0d0)))*
     &  gi*(2.0d0*dcos(q)**2-dsin(q)**2)

        vq1(i,j)=v0*(10.0d0/4.1359d0)*(pm/p)**3*(-1.0d0+
     &  cc1*(gi/gi0)**cm-cc2*(gi/gi0)**(cm+cp))*dsin(q)*dcos(q)
           if(p.le.pp)then 
           vp1(i,j)=0.0d0
           vq1(i,j)=0.0d0
           endif
           deta(i,j)=(etab2(i+1,j)-etab2(i-1,j))/(dp)
           vpb(i,j)=(vp1(i,j)-deta(i,j))*dt/(2.d0*p*dp)
           vp(i,j)=vp1(i,j)*dt/(2.d0*p*dp)
           vq(i,j)=vq1(i,j)*dt/(2.d0*p*dsin(q)*dq)
!          write(87,*)p/pm,q*180.0/pi, vp1(i,j), vq1(i,j)
        end do
!          write(87,*)p/pm,vp1(i,368)
!          write(87,*)p/pm,vq1(i,383)
        end do
!        close(87)
C
c------------------------------------------------------------
C
C	PART II-A.  Legendre polynomials are calculated here with the
C	help of a subroutine and stored for future use.
C
           do l=1,96
            cl=1.0d0/(1.0d0+float(l)/float(l+1)*(pm/pw)**(2*l+1))
            ss1_p(l)=float(2*l+1)/(float(l)*pm)*(cl-
     &        float(l)/float(l+1)*(1.0d0-cl))
        do j=1,n+1
        qm=4.0d0*atan(1.0d0)
        q=qm+float(j-1)*dq
        x=dcos(q)
        pl(j,l) = plgndr(l,1,x)
        sn(j) = dsin(q)
        end do
           end do
C
C------------------------------------------------------------------
C
C	PART III. This is the part where initial values of the 
C	variables u(i,j) and ub(i,j) are chosen, before beginning 
C	the time advancement.

 42    format(5(f13.7,1x))
       if(kk.eq.1)then
c bidya
	if(irelax.eq.0) then
	  do i = 1,n+1
	  do j = 1,n+1
            p=pb+float(i-1)/float(n)*(pm-pb)
            q=qm-float(j-1)/float(n)*qm
	  u(i,j) = 0.0d0
	  ub(i,j) = 1.0*dsin(2.0d0*q)*dsin(pi
     &      *((p-pb)/(pm-pb)))
	  end do
	  end do
	  go to 10
	end if

        open(12,file='init.dat',status='unknown')
        do i=1,n+1
          p=pb+(pm-pb)*float(i-1)/float(n)
         do j=1,n+1
          read(12,42) q,p,uin(i,j),ub(i,j)
           if(j.eq.n+1) then
            u(i,j)=0.0d0
           else
c            if (p.gt.0.8*pm) then
c            u(i,j)=fac*uin(i,j)/(p*dsin(q))
c            else
            u(i,j)=uin(i,j)/(p*dsin(q))
c            end if
           end if
         end do
        end do
        close(12)
c        open(14,file='u.dat',status='unknown')
c        do i=1,n+1
c         do j=1,n+1
c          read(14,42) q,p,u(i,j),junk1,junk2
c           if(j.eq.n+1) then
c            u(i,j)=0.0d0
c           else
c            u(i,j)=u(i,j)/(p*dsin(q))
c           end if
c         end do
c        end do
c        close(14)

 44     format(4(f13.5,1x))
        endif
c by bidya
C
C----------------------------------------------------------------
C
C	PART IV. For advancing in time, we need some matrix 
C	elements arising out of the difference schemes used. We 
C	calculate and store the matrix elements now.  Look at 
C	Sect. 5 of the 'Guide' for a discussion about these 
C	matrix elements.  Only Level III Users need to bother
C	about the matrix elements.
C
10	do i=2,2*n
	do j=2,2*n	          
	  	dvp(i,j)=(vp1(i+1,j)-vp1(i-1,j))/(dp)
	end do
	end do  
       do i=2,n
         do j=2,n
           p=pb+float(i-1)/float(n)*(pm-pb)
           q=qm-float(j-1)/float(n)*qm
        a(i,j)=-(eta(i,j)*dt/(2.0d0*(p*dq)**2)-eta(i,j)*dt/
     & (4.0d0*dtan(q)*p*p*dq)+
     & vq(2*i-1,2*j-1)*dsin(q-dq/
     & 2.0d0)*(1.0d0+vq(2*i-1,2*j-2)*dsin(q-dq))/2.0d0)
        b(i,j)=-(-eta(i,j)*dt/((p*dq)**2)-
     & 0.0*eta(i,j)*dt/(2.0d0*dtan(q)*p*p*dq)-
     & eta(i,j)*dt/(4.0d0*(p*dsin(q))**2)-vq(
     & 2*i-1,2*j-1)*(dsin(q+dq/2.0d0)*(1.0d0+vq(2*i-1,2*j)*dsin(q))-
     & dsin(q-dq/2.0d0)*(1.0d0-vq(2*i-1,2*j-2)*dsin(q)))/2.0d0)
        c(i,j)=-(eta(i,j)*dt/(2.0d0*(p*dq)**2)+
     & eta(i,j)*dt/(4.0d0*dtan(q)*p*p*
     & dq)-vq(2*i-1,2*j-1)*dsin(q+dq/2.0d0)*(1.0d0-vq(2*i-1,2*j)*dsin(q+
     & dq))/2.0d0)
        d(i,j)=-(eta(i,j)*dt/(2.0d0*dp**2)-eta(i,j)*dt/(2.0*p*dp)+
     &  vp(2*i-1,2*j-1)*(p-dp/2.0d0)*
     & (1.0d0+vp(2*i-2,2*j-1)*(p-dp))/2.0d0-0.0d0*eta(i,j)*dt/(p*dp))
        e(i,j)=-(-eta(i,j)*dt/(dp**2)-
     & 0.0d0*eta(i,j)*dt/(p*dp)-eta(i,j)*dt/(4.0d0*
     & (p*dsin(q))**2)-vp(2*i-1,2*j-1)*(dp+p*((p+dp/2.0d0)*vp(2*i,2*j-
     & 1)+(p-dp/2.0d0)*vp(2*i-2,2*j-1)))/2.0d0)
        f(i,j)=-(eta(i,j)*dt/(2.0d0*dp**2)+
     & eta(i,j)*dt/(2.0*p*dp)-vp(2*i-1,2*j-1)*
     & (p+dp/2.0d0)*(1.0d0-vp(2*i,2*j-1)*(p+dp))/2.0d0)
        end do
        end do
C
         do i=2,n
         do j=2,n
           p=pb+float(i-1)/float(n)*(pm-pb)
           q=qm-float(j-1)/float(n)*qm
        ab(i,j)=-(etab(i,j)*dt/(2.0d0*(p*dq)**2)-etab(i,j)*dt/
     & (4.0d0*dtan(q)*p*p*dq)+
     & vq(2*i-1,2*j-2)*dsin(q-dq/
     & 2.0d0)*(1.0d0+vq(2*i-1,2*j-3)*dsin(q-dq))/2.0d0)
        bb(i,j)=-(-etab(i,j)*dt/((p*dq)**2)-
     & etab(i,j)*dt/(4.0d0*(p*dsin(q))**2)-(vq(
     & 2*i-1,2*j)*dsin(q+dq/2.0d0)*(1.0d0+vq(2*i-1,2*j-1)*dsin(q))
     &-vq(2*i-1,2*j-2)*dsin(
     & q-dq/2.0d0)*(1.0d0-vq(2*i-1,2*j-1)*dsin(q)))/2.0d0)
        cb(i,j)=-(etab(i,j)*dt/(2.0d0*(p*dq)**2)+
     & etab(i,j)*dt/(4.0d0*dtan(q)*p*p*
     & dq)-vq(2*i-1,2*j)*dsin(q+dq/2.0d0)*(1.0d0-vq(2*i-1,2*j+1)*dsin(q+
     & dq))/2.0d0)
        db(i,j)=-(etab(i,j)*dt/(2.0d0*dp**2)-etab(i,j)*dt/(2.d0*p*dp)+
     & vpb(2*i-1,2*j-1)*(p-dp/2.0d0)*(
     & 1.0d0 +vpb(2*i-2,2*j-1)*(p-dp))/(2.0d0))
        eb(i,j)=-(-etab(i,j)*dt/(dp**2)-
     & 0.0d0*etab(i,j)*dt/(p*dp)-etab(i,j)*dt/(4.0d0*
     & (p*dsin(q))**2)-dvp(2*i-1,2*j-1)*dt/2.0d0-vpb(2*i-1,2*j-1)*(dp+
     & vpb(2*i,2*j-1)*p*(p+dp/2.0d0)+(p-dp/2.0d0)*p*
     & vpb(2*i-2,2*j-1))/2.0d0)
        fb(i,j)=-(etab(i,j)*dt/(2.0d0*dp**2)+
     & etab(i,j)*dt/(2.d0*p*dp)-vpb(2*i-1,2*j-1)*
     & (p+dp/2.0d0)*(1.0d0 -vpb(2*i,2*j-1)*(p+dp))/2.0d0)
          end do
          end do
C
C-------------------------------------------------------------------------
C
C	PART V. Now we come to the central part of the programme where the
C	time advancement takes place. 'k' is the counter to keep tab on
C	the time.  In each step of the do loop 'do k=1,kend', the magnetic
C	fields are advanced through one time step 'dt'.

        kend = tmax/dt

!        t=0.0d0

        do k=1,kend
c        if (t.gt..9461d0) then
c         al0=0.0d0
c        end if
!        ubolds=ub(90,70)
!        uboldn=ub(90,190)
C
C PART V-A. CALCULATING TIME ADVANCED MAGNETIC FIELDS AT INTERIOR GRID 
C POINTS. Sect. 5 of the 'Guide' explains how magnetic fields are 
C advanced by solving tridiagonal matrices with the subroutine 'tridag'
C
         do j=2,n
           aver=0.0d0
          do i=n1,n12
            aver=aver+ub(i,j)
          end do
            ubave(j)=aver/float(n12-n1+1)
         end do

          do i=2,n
           do j=2,n
            phi(i,j)=-d(i,j)*u(i-1,j)+(-e(i,j)+1.0d0)*u(i,j)
     &        -f(i,j)*u(i+1,j)+(al(i,j)*ubave(j)*dt/2.0d0)*0.25d0*
     &      (1.0d0+ erf((ubave(j))**2 - 1.0d0))*
     &	    (1.0d0 - erf((ubave(j))**2 - 49.0d0))
           end do
          end do

          do i=2,n
           do j=2,n
            p=pb+float(i-1)/float(n)*(pm-pb)
            q=qm-float(j-1)/float(n)*qm
            br = (u(i,j+1)*dsin(q+dq)-u(i,j-1)*dsin(q-dq))*dt/(4.0d0*dq)
            bt=(u(i-1,j)*(p-dp)-u(i+1,j)*(p+dp))*dsin(q)*dt/(4.0d0*p*dp)
            phib(i,j)=-db(i,j)*ub(i-1,j)+(-eb(i,j)+1.0d0)*ub(i,j)
     &        -fb(i,j)*ub(i+1,j)+br*dror(i,j)+bt*drot(i,j)
           end do
          end do

          do i=2,n
           do j=2,n
            a1(j-1)=a(i,j)
            b1(j-1)=b(i,j)+1.0d0
            c1(j-1)=c(i,j)
            r(j-1)=phi(i,j)
           end do
            r(1)=phi(i,2)-a1(1)*u(i,1)
            r(n-1)=phi(i,n)-c1(n-1)*u(i,n+1)
            call tridag(a1,b1,c1,r,uu,n-1)
           do j=2,n
	    u(i,j) = uu(j-1)
            phi(i,j)=-phi(i,j)+2.0d0*uu(j-1)
           end do
          end do

          do i=2,n
           do j=2,n
            a1(j-1)=ab(i,j)
            b1(j-1)=bb(i,j)+1.0d0
            c1(j-1)=cb(i,j)
            r(j-1)=phib(i,j)
           end do
            r(1)=phib(i,2)-a1(1)*ub(i,1)
            r(n-1)=phib(i,n)-c1(n-1)*ub(i,n+1)
            call tridag(a1,b1,c1,r,uub,n-1)
           do j=2,n
	    ub(i,j) = uub(j-1)
            phib(i,j)=-phib(i,j)+2.0d0*uub(j-1)
!            phi(i,j)=phi(i,j)+al(i,j)*ub(i,j)*dt/2.0d0
           end do
          end do

         do j=2,n
           aver=0.0d0
          do i=n1,n12
            aver=aver+ub(i,j)
          end do
            ubave(j)=aver/float(n12-n1+1)
         end do

          do i=2,n
           do j=2,n
            phi(i,j)=phi(i,j)+(al(i,j)*ubave(j)*dt/2.0d0)*0.25d0*
     &      (1.0d0+ erf((ubave(j))**2 - 1.0d0))*
     &	    (1.0d0 - erf((ubave(j))**2 - 49.0d0))
           end do
          end do

         do i=2,n
          do j=2,n
            p=pb+float(i-1)/float(n)*(pm-pb)
            q=qm-float(j-1)/float(n)*qm
            br = (u(i,j+1)*dsin(q+dq)-u(i,j-1)*dsin(q-dq))*dt/(4.0d0*dq)
            bt=(u(i-1,j)*(p-dp)-u(i+1,j)*(p+dp))*dsin(q)*dt/(4.0d0*p*dp)
            phib(i,j)=phib(i,j)
     &        +br*dror(i,j)+bt*drot(i,j)
           end do
          end do

          do j=2,n
           do i=2,n
            a1(i-1)=d(i,j)
            b1(i-1)=e(i,j)+1.0d0
            c1(i-1)=f(i,j)
            r(i-1)=phi(i,j)
           end do
            r(1)=phi(2,j)-a1(1)*u(1,j)
            r(n-1)=phi(n,j)-c1(n-1)*u(n+1,j)
            call tridag(a1,b1,c1,r,uu,n-1)
           do i=2,n
            u(i,j)=uu(i-1)
           end do
          end do

          do j=2,n
           do i=2,n
            a1(i-1)=db(i,j)
            b1(i-1)=eb(i,j)+1.0d0
            c1(i-1)=fb(i,j)
            r(i-1)=phib(i,j)
           end do
            r(1)=phib(2,j)-a1(1)*ub(1,j)
            r(n-1)=phib(n,j)-c1(n-1)*ub(n+1,j)
            call tridag(a1,b1,c1,r,uub,n-1)
	   do i=2,n
            ub(i,j)=uub(i-1)
	   end do
          end do

C PART V-B. BOUNDARY CONDITIONS. These are discussed in Sect. 2.1 of .
C the 'Guide'. 
C Boundary condition at bottom
          do j=1,n+1
           u(1,j)=0.0d0
           ub(1,j)=0.0d0
          end do
C Boundary condition at top surface
          do nup=1,11
           do l=1,96
            ss1(l)=ss1_p(l)*ss(n,dq)
           end do
           do j=1,n+1
            call coeff(n,ss1,dq,j,cf)
            u(n+1,j)=u(n,j)+cf*dp
           end do
          end do
	  do j=1,n+1
	    ub(n+1,j)=0.0d0
	  end do
C Boundary conditions at the poles
          do i=2,n+1
           u(i,n+1)=0.0d0
           ub(i,n+1)=0.0d0
           u(i,1)=0.0d0
           ub(i,1)=0.0d0
          end do
C
C PART V-C.  INCORPORATING THE MAGNETIC BUOYANCY by radial transport
C of toroidal field exceeding a specified critical value:
!	   tau = .0088d0
!	   nmb = tau/dt
!           kmb = k/nmb
!           kch = nmb*kmb
!           if(k.eq.kch)then
!	     n2 = n/2.
!	     nt = n - 6.		
!           do i=1,n2
!           do j=2,n
!            if (ra(i).gt.0.71*pm)then
!            if (ABS(ub(i,j)).GT.0.8d0)then
!             jer=j
!             ier=i
!             qjer=qm-float(jer-1)*qm/float(n)	
!             ber=ub(ier,jer)
!             bold=ub(nt,jer)
!             ub(nt,jer)=((ra(ier)/ra(nt))*0.5d0*ber) + bold
!             ub(ier,jer)=0.5d0*ber
!              open(25,file='ber.dat',status='unknown',access='append')
!              write(25,47) t,qjer,ber
!  47          format(3(d15.8,1x))
!              close(25)
!            end if
!            end if
!           end do
!           end do
!           end if
c      if (k/500*500.eq.k) then
c      open(10,file='snapshot.dat',status='unknown',access='append')
c      do i=1,n+1
c          do j=1,n+1
c			p=pb+float(i-1)/float(n)*(pm-pb)
c			q=qm-float(j-1)/float(n)*qm
c			write(10,42) t,q,p,p*dsin(q)*u(i,j),ub(i,j)
c          end do
c      end do
c      close(10)
c      end if

          t=t+dt
C
C PART V-D. WRITING IN THE FILES 'rad.dat' and 'butbot.dat' after 
C certain intervals of time
C
        if(t.gt.0.0d0) then
	  p20 = pi/20.
	  istep = t/p20
	  tdiff = t - float(istep)*p20
 	    if(tdiff.ge.dt)go to 55
  	  open(17,
     &file='rad.dat',status='unknown',access='append') 
  	  open(18,
     &file='butbot.dat',status='unknown',access='append') 
c          open(19,
c    &file='apot.dat',status='unknown',access='append')
 	 do j = 2, n
           aver=0.0d0
          do i=n1,n12
            aver=aver+ub(i,j)
          end do
            ubave(j)=aver/float(n12-n1+1)
            q=qm-float(j-1)/float(n)*qm
 	    br = -(u(n-1,j-1)*dsin(q-dq) - u(n-1,j+1)
     & *dsin(q+dq))/(2.0d0*dq*dsin(q))/pm
 	    write(17,37) t, q, br
 	    write(18,37) t, q, ubave(j)
c 	    write(19,37) t, q, u(122,j)
 	  end do
  37      format(3(f13.7,1x))
 	  close(17)
 	  close(18)
c 	  close(19)

          do j=115,128
           qs(j)=qm-float(j)/float(n)*qm
          enddo
            flux=0.0d0

          do j=115,128
            q=qm-float(j-1)/float(n)*qm
            brs =-(u(n-1,j-1)*dsin(q-dq) - u(n-1,j+1)
     & *dsin(q+dq))/(2.0d0*dq*dsin(q))/pm
       flux=flux+2.d0*pi*(pm**2)*(dcos(0.5d0*(qs(j)+qs(j+1)))-
     &  dcos(0.5d0*(qs(j-1)+qs(j))))*brs
          enddo

          open(20,file='rad_field',status='unknown',access='append')
          write(20,37) t, flux
          close(20)

  55    end if	
C
C WRITING IN THE FILE 'run.dat' after every 40 time steps
C
        if(k/400*400.eq.k)then
          open(95,file='series_b',status='unknown',access='append')
         fac=1.0d8/(3600.*24*365)
          write(95,49)(t*fac),abs(ub(ita,74))
 49    format(2(f13.7,1x))
          close(95)
        end if
!         open(95,file='bnbs.dat',status='unknown',access='append')
!          write(95,90) ub(nr, mthn), ub(nr, mths)
!  90      format(2(d15.9,1x))
!         close(95)
!        end if
!	if(k/40*40.eq.k)then
c 	  open(95,
c     &file='run.dat',status='unknown',access='append')
c	  write(95,90)t,ub(90,74*2),ub(90,70),ub(90,95*2),
c     &   u(240,122*2),u(240,12)
!  90      format(6(d15.9,1x))
c	  close(95)
!	end if
c        if (ub(90,70)*ubolds.lt.0.0d0) then
c        open(96,file='szero.dat',status='unknown',access='append')
c        write(96,*) t
c        close(96)
c        end if
c        if(ub(90,190)*uboldn.lt.0.0d0) then
c       open(97,file='nzero.dat',status='unknown',access='append')
c        write(97,*) t
c        close(97)
c        end if
C here I calculate the polar field outside the solar surface.
!      open(10,file='draw.dat',status='unknown',access='append')
!         fac=1.0d8/(3600.*24.*365.)
!          aw=0.0d0
!         do l=1,96,2
!       ssp=pl((n+2)/2,l)*ss(n,dq)
!       clc=1.0d0/(1.0d0+float(l)/float(l+1)*(pm/pw)**(2*l+1))
!       facl=(float(2*l+1)/float(l*(l+1)))*((pm/pw)**(l+1))*
!     & (clc+((1.d0-clc)*(pw/pm)**(2*l+1)))
!         aw=aw+facl*ssp
!         end do
!           write(10,49) t*fac,aw
! 49    format(2(f13.7,1x))
!           close(10)
C
C    WRITING TOROIDAL FLUX IN TACHOCLINE:
C    from 10 to 45deg latitude and r=0.677R_sun to 0.726R_sun
        if(k/400*400.eq.k)then
          fluxn=0.0d0
          fluxp=0.0d0
          fluxt=0.0d0
          do i=36,50
            do j=72,97
               p=pb+float(i-1)/float(n)*(pm-pb)
               q=qm-float(j-1)/float(n)*qm
               p1=pb+float(i)/float(n)*(pm-pb)
               q1=qm-float(j)/float(n)*qm
               if(ub(i,j).lt.0.0d0) then
                  fluxn=fluxn+ub(i,j)*0.5d0*(q-q1)*(p1**2-
     &                 p**2)
               else
                  fluxp=fluxp+ub(i,j)*0.5d0*(q-q1)*(p1**2-
     &                 p**2)
               end if
            end do
          end do
          fluxt=fluxn+fluxp
         open(21,file='torflux_dat',status='unknown',access='append')
         write(21,49) t, fluxt
         close(21)
        end if

        if(k/100000*100000.eq.k)then
         open(13,file='final_temp',status='unknown')
         do i=1,n+1
          do j=1,n+1
           p=pb+float(i-1)/float(n)*(pm-pb)
           q=qm-float(j-1)/float(n)*qm
           write(13,42) q,p,p*dsin(q)*u(i,j),ub(i,j)
          end do
         end do
         close(13)
         endif

         end do
C The mammoth do loop which began at the beginning of PART V ends
C here. 
         end do
c by bidya
C
C---------------------------------------------------------------
C	PART VI. Now the final state is written in the file 
C	'final.dat'
C
         open(12,
     &file='final.dat',status='unknown')
         do i=1,n+1
          do j=1,n+1
           p=pb+float(i-1)/float(n)*(pm-pb)
           q=qm-float(j-1)/float(n)*qm
           write(12,42) q,p,p*dsin(q)*u(i,j),ub(i,j)
          end do
         end do
         close(12)
	
        stop
        end

C
C     MAIN PROGRAM ENDS HERE.
C -----------------------------------------------------
C
C     SUBROUTINE FOR INVERTING TRIDIAGONAL MATRIX :
        subroutine tridag(a,b,c,r,u,n)
        implicit real*8 (a-h,o-z)
        parameter (nmax=129)
        dimension gam(nmax),a(nmax),b(nmax),c(nmax),r(nmax),u(nmax)
        if(b(1).eq.0.d0) stop
        bet=b(1)
        u(1)=r(1)/bet
        do 11 j=2,n
          gam(j)=c(j-1)/bet
          bet=b(j)-a(j)*gam(j)
          if(bet.eq.0.d0) stop
          u(j)=(r(j)-a(j)*u(j-1))/bet
  11    continue
        do 12 j=n-1,1,-1
          u(j)=u(j)-gam(j+1)*u(j+1)
  12    continue
        return
        end


C
C     THREE SUBROUTINES NEEDED FOR UPPER BOUNDARY CONDITION :
C       The subroutine plgndr calculates Legendre polynomials.
        subroutine coeff(n,ss1,dq,j,cf)
        implicit real*8(a-h,o-z)
        external plgndr
        parameter (nmax=129,lmax=96)
        dimension ss1(nmax)
        common u(nmax,nmax),t,it
        common/lmx/l
        common /plsincos/pl(nmax,lmax),sn(nmax)
        cf=0.0d0
        do l=1,96
        cf=cf+ss1(l)*pl(j,l)
        end do
        return
        end



        function ss(n,dq)
        implicit real*8(a-h,o-z)
        external plgndr
        parameter (nmax=129,lmax=96)
        dimension y(nmax)
        common u(nmax,nmax),t,it
        common/lmx/l
        common /plsincos/pl(nmax,lmax),sn(nmax)
        do j=1,n+1
        y(j)=pl(j,l)*u(n+1,j)*sn(j)
        end do
        term1=0.0d0
        do k1=2,n,2
        term1=term1+y(k1)
        end do
        term2=0.0d0
        do k2=3,n-1,2
        term2=term2+y(k2)
        end do
        ss=dq/3.0d0*(y(1)+4.0d0*term1+2.0d0*term2+y(n+1))
        return
        end


        function plgndr(l,m,x)
        implicit real*8(a-h,o-z)
        if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.0d0)stop 'bad arguments'
        pmm=1.0d0
        if(m.gt.0) then
          somx2=dsqrt((1.0d0-x)*(1.0d0+x))
          fact=1.0d0
          do 11 i=1,m
             pmm=-pmm*fact*somx2
             fact=fact+2.0d0
  11      continue
        end if
        if(l.eq.m) then
          plgndr=pmm
        else
          pmmp1=x*(2*m+1)*pmm
          if(l.eq.m+1) then
             plgndr=pmmp1
          else
             do 12 ll=m+2,l
                 pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
                 pmm=pmmp1
                 pmmp1=pll
  12         continue
           plgndr=pll
          end if
        end if
        return
        end

c
      function erf (x)
      implicit real*8 (a-h,o-z)
      data p/0.3275911d0/,a1/0.254829592d0/,a2/-0.284496736d0/,
     &     a3/1.421413741d0/,a4/-1.453152027d0/,a5/1.061405429d0/
      sig=sign(1.0d0,x)
      xabs=dabs(x)
      if (xabs.gt.20) then
c     
c     *** vermeidung von underflow ***
c
       erf=sig
      else
       t=1.0d0/(1.0d0+p*xabs)
       erf=1.0d0-((((a5*t+a4)*t+a3)*t+a2)*t+a1)*t*exp(-xabs*xabs)
       erf=erf*sig
      endif
      return
      end
