            program AAperiph
c 24.3.14
c a geometric model giving wounded nuclei for AA collisions
c using simplified hard spheres
          parameter( nbin=100 )
          common /seed/ iseed          
          open(unit=1,file='inAA.dat',status='unknown')
          open(unit=2,file='outAA.dat',status='unknown')

          pi=4*atan(1.0)
          iseed =-56789

c we will get b in the range of centralities
         read(1,*)A,niter,b,shift
c niter is number of points to be selected
c           niter=100000
            rnucl=A**.333*1.2
             rnucls=rnucl**2
c sigman is sigma total times density in fm -1
             sigman=1.6
c the niter is here the number of points
           do 1 iter=1,niter

c point x,y is where the interaction happens
           call givepoint(x,y,rnucl,b,sigman)
            call givestrings(x1,y1,x2,y2,x,y,shift,pi)
            write(2,*) x1,y1,x2,y2,iter

c end of the main loop
 1         continue

            stop
             end
c--------------------------------------------------------------
c---------------------------------------------------------------------+---
         subroutine givepoint(x,y,rnucl,b,sigman)

c gives a point with a collision probability for two spheres
c in the original program collisions were hard and prob scaled as product
c of two lengths, L1 and L2.
c now it is soft interaction stot , densit
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

c------------------------------------------------
c---------------------------------------------------------------------+---
c---------------------------------------------------------------------+---

      function ran2(idum)
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


