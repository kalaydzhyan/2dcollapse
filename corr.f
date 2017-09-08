c HERE WE DO THE DYNAMICAL CORRELATORS
c now the double loop over time
                  do 500 istep=ntcorr+1 ,nstep 
                  do 500 icorr=1,ntcorr
c so one time is istep, the other is istep-icorr+1

c let me start with self-correlation of one particle, van Hove's Gs
                   do 510 i=1,np
                rr=0.
c note that distances to images are NOT included now
                  rr=(xxt(istep,i,1)-xxt(istep-icorr+1,i,1))**2
     2              +(xxt(istep,i,2)-xxt(istep-icorr+1,i,2))**2
                  dist=sqrt(rr)
c now we have to put it into the right place jj in the hystogramm
c here I assume xmin=0. and take care than it is not too large
               if(dist.gt.xmax)goto 510
                jj=dist/xst+1.000001
                gs(icorr,jj)=gs(icorr,jj)+1.
                ngs(icorr)=ngs(icorr)+1
c                 write(6,*)i,icorr,dist,gs(icorr,jj)
 510          continue           
               
c now I will do van Hove's Gd function which has a sum over 2 particles
c and i is never equal to j
                do 501 i=1,np
                do 501 j=1,np
                 if(i.eq.j)goto 501
c let me now add mirror images, im=-1,0,1 in all directions
                  do 502 imx=-1,1
                  do 502 imy=-1,1

                   rr=0.
c note that distances to images are included now
         rr=(xxt(istep,i,1)-xxt(istep-icorr+1,j,1)+al(m)*imx)**2
     2     +(xxt(istep,i,2)-xxt(istep-icorr+1,j,2)+al(m)*imy)**2
                  dist=sqrt(rr)
c now the difficulty is that we have
c now we have to put it into the right place jj in the hystogramm
c here I assume xmin=0. and take care than it is not too large
               if(dist.gt.xmax)goto 501
                jj=dist/xst+1.000001
                gd(icorr,jj)=gd(icorr,jj)+1.
                ngd(icorr)=ngd(icorr)+1

 502            continue
 501            continue

 500         continue

c now the hyst for correlator should be renormalized radially
               anorm=1./4./3.14/dens
               do 560 icorr=1,ntcorr 
                do 561 ii=1,nst
                 rd=xmin+xst*(ii-.5)
                 gs(icorr,ii)=anorm*gs(icorr,ii)/(rd**2*ngs(icorr))
                 gd(icorr,ii)=anorm*gd(icorr,ii)/(rd**2*ngd(icorr))
 561            continue
 560            continue
c let me now write differently, plotting wants columns

                do 563 ii=1,nst
                 rd=xmin+xst*(ii-.5)
               write(11,564)rd,(gs(icorr,ii),icorr=1,ntcorr)
               write(14,564)rd,(nstep*gd(icorr,ii),icorr=1,ntcorr)
  564          format(8f8.4)
  563           continue

C E 
c also making probablity function from hyst array 
               anorm=1./4./3.14/dens
               do 66 ii=1,100
                rd=xmin+xst*(ii-.5)
                eq=anorm*ieq(ii)/rd**2
               aneq=anorm*ineq(ii)/rd**2     
               write(17,*)rd,eq,aneq
 66            continue

               stop

c E: this is skipped for the time being
             
c density-density correlation function f(t)


             do 350 n=0,ncor
             correl(n+1)=(0.0,0.0) 
             correlch(n+1)=(0.0,0.0)

             do 450 istep=n+1,nstep+1
             correl(n+1)=correl(n+1)+rho(istep)*conjg(rho(istep-n))/np
     $       /(nstep+1-n)
             correlch(n+1)=correlch(n+1)
     $       +rhoch(istep)*conjg(rhoch(istep-n))/np/(nstep+1-n)

c! check nstep+1-n


 450         continue
             
             t=n*tstep
             write(17,*)t,real(correl(n+1)),aimag(correl(n+1))

 350         continue
