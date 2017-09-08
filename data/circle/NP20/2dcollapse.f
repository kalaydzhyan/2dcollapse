#define    ndim    2
c np must be even:
#define    np      20
#define    niter   (np/2)
#define    nstep   80
#define    nv      (2*ndim*np)
#define    nt      102
#define    tension (0.42*0.42)			 
#define    hc      0.197326968
#define    tstep   0.1
#define    ssize   (0.11*1.26)
c String tension (in GeV^2), string width in fm


            program spaghetti

c one should compile this program together with specfun.f

            implicit double precision (a-h,o-z)

c global parameter particle and var number  should be defined here 
c np should be even
c dw,nw the step length and the number of steps 
c in the frequency interval

c  ndim=2 is the new number of dimensions
c ****************** new ********************

c nv is the total number of variables for solver

c             complex rho(nt),correl(ncor),rhoch(nt),correlch(ncor)
              real r,dist,xmin,xst,ekinh,rh              
    
             external SUB
	      external besk1
             external besk0

             common /coor/ xx(np,ndim),vv(np,ndim)
             common /a/   coupl, smass 
             common /seed/ iseed
             common /hyst/ ieq(100),ineq(100),irr(100)
             common /nucleus/ A,b,shift


c iekin will hold a histogramm over kinetic energy
             common /hist/ ipp(100),ipm(100),iekin(100)
             dimension y(nv),f(nv),w(nv,6)

c            randomizer seed in form milliseconds*100 + seconds of
c            current time
             dimension idate_time(8)
             character*10 bb(3)
             call date_and_time(bb(1), bb(2), bb(3), idate_time) 
             iseed = idate_time(8)*100+idate_time(7)

c clear histograms
               do 60 nnn=1,100
               ipp (nnn)=0
               ipm(nnn)=0               
               iekin(nnn)=0
 60            continue
c clear  hyst
c hyst parameters are single precision
                xmin=0.0
                xst=0.1
                nst=40
c                  xmax=xmin+nst*xst

             open(10,file='in.dat', status='unknown')
             open(11,file='outxv.dat', status='unknown')
             open(12,file='outenergy.dat', status='unknown')

c note nstep is now not read but given in param
             read(10,*) smass,coupl,h0,eps,ninit

c the second read if init=2 for AA
            if(ninit.eq.2) then 
            read(10,*)A,b,shift
            end if
             smass = smass / hc
             coupl = coupl 

c we can do friction but for now
               cool=0
               t=0.
              pi=4*atan(1.)


c initialization subroutine
c             ninit=1
             call init(ninit)
             call conserv(t)

c this is the beginnig of the main loop over number of steps

             do 1 istep=0,nstep

c main loop over steps, followed by loop over particles and components

             t1=tstep*istep
             t2=tstep*(istep+1)

c now we load x,v  into one array y

               do 20 i=1,np
               do 21 m=1,ndim 
               y(m + 2*ndim*(i-1))=xx(i,m)
               y(m + ndim + 2*ndim*(i-1))=vv(i,m)
 21             continue
 20             continue

             call DDEQMR(nv,t1,t2,y,h0,eps,sub,w)
                t=t2
c we load xx and vv back

               do 4 i=1,np
               do 5 m=1,ndim 
               xx(i,m)=y(m + 2*ndim*(i-1))
               vv(i,m)=y(m + ndim + 2*ndim*(i-1))
 5          continue
 4          continue

               call conserv(t)

c instanteneous coordinates and velocities are recorded 
c into file 'outxv.dat'
            do 73 ip=1,np
        write(11,116)xx(ip,1),xx(ip,2),vv(ip,1),vv(ip,2),t,ip
 116      format(5f10.5,I3)
 73          continue


 1            continue
c the end of the main cycle


c writing final histograms
             call lev(0.,0.005,40,2,iekin)
             call lev(0.,0.1,40,2,ipp)
             call lev(0.,0.1,40,2,ipm)

cccccccccccccccccccccccccccccccccccccccccccccccccccc
c     HERE were correlators, saved to corr.f       c
cccccccccccccccccccccccccccccccccccccccccccccccccccc
               stop
               end


cccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccc

              subroutine conserv(t)
             implicit double precision (a-h,o-z)
              real ekin,ekinh
              common /hist/ ipp(100),ipm(100),iekin(100)
                common /coor/ xx(np,ndim),vv(np,ndim)
                    ekin=0.
                do 1 i=1,np                
                     ekinh=0.
                do 33 m=1,ndim
                  ekinh=ekinh+vv(i,m)**2/2
 33            continue
c now put into a histogram                  
                  call lens(ekinh,0.,.005,40,iekin)
                  ekin=ekin+ekinh
 1             continue
c now pot. energy 
                 epot=0.
                do 3 i=1,np
                do 35 j=i+1,np
c                if(i.eq.j) goto 3          
                epot=epot+v(i,j)
 35              continue
  3             continue 
                e=ekin+epot
             write(12,*)t,e,ekin,epot
                return
                end

ccccccccccccccccccccccccccccccccccccccccccccccccc

             subroutine init(ninit)
             implicit double precision (a-h,o-z)
             common /coor/ xx(np,ndim),vv(np,ndim)
             common /nucleus/ A,b,shift

             pi = 4.0*atan(1.0)
             if(ninit.eq.1) then
                  bmax=3.
c particle loop
                   do 1 i=1,np
                  randphi=2*pi*rang()
 2                  b=rang()*bmax
c       prof2 is squared profile of the elastic scattering
                  p=prof2(b)
                  ran=rang()
                 if(ran.gt.p)goto 2
c it will repeat it till it finds the point below the curve
                  xx(i,1)=b*cos(randphi)
                  xx(i,2)=b*sin(randphi)
                  vv(i,1)=0.
                  vv(i,2)=0.
 1              continue
             else
ccccccccccccccccccc             
c ANOTHER GEOMETRY, namely AA peripheral is HERE
c            program AAperiph
c 24.3.14
c a geometric model giving wounded nuclei for AA collisions
c using simplified hard spheres
c          parameter( nbin=100 )
c         common /seed/ iseed          
c         open(unit=1,file='inAA.dat',status='unknown')
c          open(unit=2,file='outAA.dat',status='unknown')


c we will get b in the range of centralities

c niter is number of points to be selected
            rnucl=A**.333*1.2
             rnucls=rnucl**2
c sigman is sigma total times density in fm -1
             sigman=1.6
c the niter is here the number of points
           do 10 iter=1,niter

c point x,y is where the interaction happens
           call givepoint(x,y,rnucl,b,sigman)
            call givestrings(x1,y1,x2,y2,x,y,shift,pi)
                  xx(1+2*(iter-1),1)=x1
                  xx(1+2*(iter-1),2)=y1
                   xx(2+2*(iter-1),1)=x2
                  xx(2+2*(iter-1),2)=y2
                  vv(1+2*(iter-1),1)=0.
                  vv(1+2*(iter-1),2)=0.
                  vv(2+2*(iter-1),1)=0.
                  vv(2+2*(iter-1),2)=0.
 10              continue


 
             end if
              return
             end
c--------------------------------------------------------------
c---------------------------------------------------------------------+---
         subroutine givepoint(x,y,rnucl,b,sigman)

c gives a point with a collision probability for two spheres
c in the original program collisions were hard and prob scaled as product
c of two lengths, L1 and L2.
c now it is soft interaction stot , densit

             implicit double precision (a-h,o-z)

 1        continue     
           rnucls=rnucl**2  
           x=2*rnucl*(2*rang()-1.)
           y=rnucl*(2*rang()-1.)
           aL1s=(rnucls-(x-b/2)**2-y**2)
           aL2s=(rnucls-(x+b/2)**2-y**2)
c      write(6,*) x,y,al1s,al2s
c the outL will need some coefficient later
           if(aL1s.lt.0.or.aL2s.lt.0)goto 1
c now the point is in the almond

c sigman is sigma total times density in fm-1

             P=1-exp(-2*sqrt(aL2s)*sigman)
             ran=rang()
             if(P.lt.ran)goto 1
          RETURN
          END
c---------------------------------------------------------------------+---
         subroutine givestrings(x1,y1,x2,y2,x,y,shift,pi)
             implicit double precision (a-h,o-z)
c x,y input, others output, shift is the radius of the displacement
          randphi=2*pi*rang()
           dx=shift*cos(randphi)
           dy=shift*sin(randphi)
           x1=x+dx/2
           x2=x-dx/2
           y1=y+dy/2
           y2=y-dy/2
           return
           end

c---------------------------------------------------------------------+----



cccccccccccccccccccccccccccccccccccccccccccccccc
        function prof2(b)

        implicit double precision (a-h,o-w)
        implicit complex (X)

c       we want to have b in GeV^-1, not fm:
	bGeV = b / hc

	XI = (0, 1)
	XR = (1, 0)

        s=7000**2
        pi = 4.0*atan(1.0)
        c = .167
        cp = .748
        XSs = s**c/(log(s)**cp)+s**c*exp(-XI*Pi*c)/((log(s)-XI*Pi)**cp)
        FF = 1.031925910*bGeV*besk1(0.577*bGeV)
     $  -.2119478974*besk0(0.577*bGeV)
     $  +3.678319512*bGeV*besk1(1.719*bGeV)
     $  -17.76215332*besk0(1.719*bGeV)
     $  +19.36687661*besk0(1.858*bGeV)
	
        prof2 = b*realpart(XR - exp(-XSs*FF))
c        prof2 = prof2**2
       
           return
           end

ccccccccccccccccccccccccccccccccccccccccccccccccccc
              function v(i,j)
            implicit double precision (a-h,o-z)
            real rh,ekinh
 
             common /coor/ xx(np,ndim),vv(np,ndim)
             common /charge/ ich(np) 
              common /hist/ ipp(100),ipm(100),iekin(100)
              dimension xi(ndim),xj(ndim)
                     do 10 m=1,ndim
                  xi(m)=xx(i,m)
                  xj(m)=xx(j,m)
  10              continue
                   v=pot(xi,xj,R)

             return
                end 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                  function pot(x,y,r)
c x and y are two 3-points, r will be returned for hysto
              implicit double precision (a-h,o-z)
               common /a/ coupl, smass
               dimension x(ndim),y(ndim)
                r=sqrt((x(1)-y(1))**2+(x(2)-y(2))**2 + ssize**2)
                pot=-coupl*besk0(smass*r)*2

             return
              end

ccccccccccccccccccccccccccccccccccccccccccccccccccc
              subroutine force(i,m,y,for)
            implicit double precision (a-h,o-z)

            dimension y(nv),xj(ndim),xi(ndim)
            
            for=0.0

                 do 120 mm=1,ndim
                 xi(mm)=y(mm+2*ndim*(i-1))
 120              continue

                
                do 11 j=1,np
c no self force
                  if(i.eq.j) go to 11

                 do 12 mm=1,ndim
                 xj(mm)=y(mm+2*ndim*(j-1))
 12              continue
                  for=for+ff(m,xi,xj)
  11           continue    
 10                continue 
            return
            end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                   function ff(m,x,y)
c it is one component
              implicit double precision (a-h,o-z)
               common /a/ coupl,smass
                  dimension x(ndim),y(ndim)

                rs=(x(1)-y(1))**2+(x(2)-y(2))**2 + ssize**2
                r=sqrt(rs)
                ff=besk1(r*smass)*coupl*smass*2
                ff=ff*(-x(m)+y(m))/r

               return
               end
 
ccccccccccccccccccccccccccccccccccccccccccccccc

           subroutine SUB(t,y,f) 
          implicit double precision (a-h,o-z)
           dimension y(nv),f(nv),w(nv,6)
            do 1 i=1,np

c             do 2 m=1,ndim
c             f(m+2*ndim*(i-1))=0.
c             f(m+ndim+2*ndim*(i-1))=0.
c            call force(i,m,y,for)
c  2          continue

             do 3 m=1,ndim
c these are x-dot=p
            f(m+2*ndim*(i-1))=y(m+ndim+2*ndim*(i-1))
c these are p-dot=force
            call force(i,m,y,for)
            f(m+ndim+2*ndim*(i-1))=for
  3         continue
 1          continue
          return
          end
cccccccccccccccccccccccccccccccccccccccccccccccc
    


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c i took the liberty to modify it and simplify all these ifs
*
* $Id: deqmr64.F,v 1.1.1.1 1996/04/01 15:02:17 mclareni Exp $
*
* $Log: deqmr64.F,v $
* Revision 1.1.1.1  1996/04/01 15:02:17  mclareni
* Mathlib gen
*
*
c#include "gen/pilot.h"
c#if !defined(CERNLIB_DOUBLE)
c      SUBROUTINE RDEQMR(N,XA,XZ,Y,H0,EPS,SUB,W)
c#endif
c#if defined(CERNLIB_DOUBLE)
       SUBROUTINE DDEQMR(N,XA,XZ,Y,H0,EPS,SUB,W)
c I give it double accuracy in my way
         implicit double precision (a-h,o-z)
c#include "gen/imp64.inc"
c#endif
C     Based on a modification of the Runge-Kutta method suggested
C     by Merson. See G.N. Lance, Numerical Methods for High speed
C     Computers, Iliffe & Sons, London 1960, pp. 56-57

c I commented name out
C      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
c#if !defined(CERNLIB_DOUBLE)
c      PARAMETER (NAME = 'RDEQMR')
c#endif
c#if defined(CERNLIB_DOUBLE)
c      PARAMETER (NAME = 'DDEQMR')
c#endif
      LOGICAL LER,LFN
c I modify this statement from f90 into the same as in sub
c       DIMENSION Y(*),W(N,*)
       dimension  y(N),f(N),w(N,6)
      PARAMETER (DELTA = 1D-14)
      PARAMETER (Z1 = 1, R2 = Z1/2, R3 = Z1/3)
      PARAMETER (R4 = 3*Z1/8, R5 = 3*Z1/2, R6 = 9*Z1/2)
      PARAMETER (R7 = 4*Z1/3, R0 = Z1/32)

c#if !defined(CERNLIB_DOUBLE)
c      ENTRY DEQMR(N,XA,XZ,Y,H0,EPS,SUB,W)
c#endif
      IF(N .LT. 1 .OR. XA .EQ. XZ .OR. H0 .EQ. 0) RETURN
      DELTAX=DELTA*ABS(XZ-XA)
      EPS5=5*ABS(EPS)
      EPS0=R0*EPS5
      X=XA
      H1=SIGN(ABS(H0),XZ-XA)
      SGH=SIGN(Z1,H1)

   12 IF(SGH*(X+H1-XZ) .LT. 0) THEN
       HH=H1
       H0=H1
       LFN=.FALSE.
      ELSE
       HH=XZ-X
       IF(ABS(HH) .LT. DELTAX) THEN
        DO 10 I = 1,N
   10   Y(I)=W(I,6)
        RETURN
       END IF
       LFN=.TRUE.
      END IF
      S2=R2*HH
      S3=R3*HH
      S7=R7*HH
      X1=X+HH
      X2=X+S2
      X3=X+S3

      CALL SUB(X,Y,W(1,1))
      DO 1 I = 1,N
      W(I,1)=S3*W(I,1)
    1 W(I,6)=Y(I)+W(I,1)

      CALL SUB(X3,W(1,6),W(1,2))
      DO 2 I = 1,N
      W(I,2)=S3*W(I,2)
    2 W(I,6)=Y(I)+R2*(W(I,1)+W(I,2))

      CALL SUB(X3,W(1,6),W(1,3))
      DO 3 I = 1,N
      W(I,3)=S3*W(I,3)
      W(I,2)=3*W(I,3)
    3 W(I,6)=Y(I)+R4*(W(I,1)+W(I,2))

      CALL SUB(X2,W(1,6),W(1,4))
      DO 4 I = 1,N
      W(I,4)=S7*W(I,4)
    4 W(I,6)=Y(I)+R5*(W(I,1)-W(I,2)+W(I,4))

      CALL SUB(X1,W(1,6),W(1,5))
      DO 5 I = 1,N
      W(I,5)=S3*W(I,5)
    5 W(I,6)=Y(I)+R2*(W(I,1)+W(I,4)+W(I,5))

      DO 8 I = 1,N
      W(I,2)=ABS(W(I,1)-R6*W(I,3)+W(I,4)-R2*W(I,5))
      W(I,1)=ABS(W(I,6))
      IF(W(I,2) .GT. EPS5*W(I,1)) THEN
       H1=R2*HH
       IF(ABS(H1) .LT. DELTAX) THEN
        WRITE(ERRTXT,101) X
c this is my comment since i do not have it
c        CALL MTLPRT(NAME,'D202.1',ERRTXT)
        RETURN
       END IF
       GO TO 12
      END IF
    8 CONTINUE
      LER=.TRUE.
      DO 7 I = 1,N
    7 LER=LER .AND. W(I,2) .LT. EPS0*W(I,1)
      DO 9 I = 1,N
    9 Y(I)=W(I,6)
      IF(LER) THEN
       H0=H1+H1
       H1=HH+HH
      END IF
      IF(LFN) RETURN
      X=X1
      GO TO 12
  101 FORMAT('TOO HIGH ACCURACY REQUIRED NEAR  X = ',1P,D15.8)
      END


c---------------------------------------------------------------------+----

      function rang()
c------------------------------------------------------------------------c
c     various random number generators                                   c
c------------------------------------------------------------------------c
      common /seed/ iseed
c     rang = rnunf()
      rang = ran2(iseed)
c     rang = ran(iseed)

      return
      end

c------------------------------------------------c---------------------------------------------------------------------+---

      function ran2(idum)
             implicit double precision (a-h,o-z)
c----------------------------------------------------------------------c
c     numerical recipes random number generator ran2 (revised version) c
c     copr. 1986-92 numerical recipes software. Reinitialize with idum c
c     negative, then do not later idum between successive calls.       c
c----------------------------------------------------------------------c

      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      real ran2,am,eps,rnmx
      parameter (im1=2147483563,im2=2147483399,am=1./im1,imm1=im1-1,
     *ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,ir2=3791,
     *ntab=32,ndiv=1+imm1/ntab,eps=1.2e-7,rnmx=1.-eps)
      integer idum2,j,k,iv(ntab),iy
      save iv,iy,idum2
      data idum2/123456789/, iv/ntab*0/, iy/0/

      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=ntab+8,1,-1
          k=idum/iq1
          idum=ia1*(idum-k*iq1)-k*ir1
          if (idum.lt.0) idum=idum+im1
          if (j.le.ntab) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/iq1
      idum=ia1*(idum-k*iq1)-k*ir1
      if (idum.lt.0) idum=idum+im1
      k=idum2/iq2
      idum2=ia2*(idum2-k*iq2)-k*ir2
      if (idum2.lt.0) idum2=idum2+im2
      j=1+iy/ndiv
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+imm1
      ran2=min(am*iy,rnmx)

      return
      end

c---------------------------------------------------------------------+---
c---------------------------------------------------------------------+---


c---------------------------------------------------------------------+---
      subroutine lev(xmin,st,n,nsm,ist)
c            implicit double precision (a-h,o-z)
c------------------------------------------------------------------------c
c     plot simple histogram in output file ftn02                         c
c------------------------------------------------------------------------c
c     input: xmin    smallest x value                                    c
c            st      bin width                                           c
c            n       number of bins                                      c
c            nsm     plot symbol                                         c
c            ist(n)  histogram array                                     c
c------------------------------------------------------------------------c

      dimension ist(n),a(6),k(6),t(50)
      data a    /1hx,1h*,1h+,1hc,1he,1h /
      data nhgt /50/

      j=0
      m=0
      s=0.
      d=xmin-st*0.5
      x=d
      do 1 i=1,n
      x=x+st
      s=s+ist(i)*x
      if(ist(i).lt.0)go to 9
      if(i.eq.1.or.i.eq.n)go to 1
      if(ist(i).gt.m)m=ist(i)
    1 j=j+ist(i)
      if(j.lt.1)go to 10
      s=s/j
      m=m/nhgt+1
      do 2 i=1,6
    2 k(i)=m*10*(i-1)
      write(2, 30)k
   30 format(5x,i20,5i10)
      write(2,40)
   40 format (2x,3h---,7(9(1h-),1h+))
      x=d
      d=0.
      do6 i=1,n
      do3 l=1,nhgt
    3 t(l)=a(6)
      nn=ist(i)/m
      if(nn.eq.0) go to 5
         if(nn.gt.nhgt)nn=nhgt
      do 4 l=1,nn
    4 t(l)=a(nsm)
    5 x=x+st
      write(2, 70)i,ist(i),x,t
   70 format(1x,1hi,i4,i7,1pe11.3,1x,50a1,1hi)
    6 d=d+ist(i)*(x-s)**2
      if(j.lt.2)goto 7
      d=sqrt(d/(j-1))
    7 write(2, 40)
      write(2, 30)k
   8  write(2, 90)j,s,d
   90 format(//2x,18hnumber of events =,i8,5x,
     *10haverage = ,e11.4,5x,8hsigma = ,e11.4//)
      return
    9 write(2, 80)i
   80 format(//2x,3hlev,i8,23h channel  less  than  0//)
      return
10    d=0.
      go to 8
      end

c---------------------------------------------------------------------+---
c---------------------------------------------------------------------+---

      subroutine lens(a,amin,st,m,ist)
c------------------------------------------------------------------------c
c     include value a in histogram array ist(n)                          c
c------------------------------------------------------------------------c
c     a      value to be added to histogram array                        c
c     amin   minimum value in histogram                                  c
c     st     bin width                                                   c
c     m      number of bins                                              c
c     ist(n) histogram array                                             c
c------------------------------------------------------------------------c

      dimension ist(150)
      j=(a-amin)/st+1.000001
      if(j.lt.1)j=1
      if(j.gt.m)j=m
      ist(j)=ist(j)+1
      return
      end

c---------------------------------------------------------------------+---
