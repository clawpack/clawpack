c
c
c     =====================================================
      subroutine rpt2(ixy,maxm,meqn,mwaves,mbc,mx,
     &                  ql,qr,aux1,aux2,aux3,
     &                  ilr,asdq,bmasdq,bpasdq)
c     =====================================================
      use geoclaw_module

      implicit none
c
c     # Riemann solver in the transverse direction using an einfeldt
c     Jacobian.

c-----------------------last modified 1/10/05----------------------

      integer ixy,maxm,meqn,mwaves,mbc,mx,ilr

      double precision  ql(meqn,1-mbc:maxm+mbc)
      double precision  qr(meqn,1-mbc:maxm+mbc)
      double precision  asdq(meqn,1-mbc:maxm+mbc)
      double precision  bmasdq(meqn,1-mbc:maxm+mbc)
      double precision  bpasdq(meqn,1-mbc:maxm+mbc)
      double precision  aux1(maux,1-mbc:maxm+mbc)
      double precision  aux2(maux,1-mbc:maxm+mbc)
      double precision  aux3(maux,1-mbc:maxm+mbc)

      double precision  s(3)
      double precision  r(3,3)
      double precision  beta(3)
      double precision  g,tol,abs_tol
      double precision  hl,hr,hul,hur,hvl,hvr,vl,vr,ul,ur,bl,br
      double precision  uhat,vhat,hhat,roe1,roe3,s1,s2,s3,s1l,s3r
      double precision  delf1,delf2,delf3,dxdcd,dxdcu
      double precision  dxdcm,dxdcp,topo1,topo3,eta

      integer i,m,mw,mu,mv
      integer maux
      common /cmmaux/  maux

      g=grav
      tol=drytolerance
      abs_tol=drytolerance

      if (ixy.eq.1) then
	     mu = 2
	     mv = 3
      else
	     mu = 3
	     mv = 2
      endif


      do i=2-mbc,mx+mbc

         hl=qr(1,i-1)
         hr=ql(1,i)
         hul=qr(mu,i-1)
         hur=ql(mu,i)
         hvl=qr(mv,i-1)
         hvr=ql(mv,i)

c===========determine velocity from momentum===========================
       if (hl.lt.abs_tol) then
          hl=0.d0
          ul=0.d0
          vl=0.d0
       else
          ul=hul/hl
          vl=hvl/hl
       endif

       if (hr.lt.abs_tol) then
          hr=0.d0
          ur=0.d0
          vr=0.d0
       else
          ur=hur/hr
          vr=hvr/hr
       endif

       do mw=1,mwaves
          s(mw)=0.d0
          beta(mw)=0.d0
          do m=1,meqn
             r(m,mw)=0.d0
          enddo
       enddo
      dxdcp = 1.d0
      dxdcm = 1.d0

       if (hl.le.drytolerance.and.hr.le.drytolerance) go to 90

*      !check and see if cell that transverse waves are going in is high and dry
       if (ilr.eq.1) then
            eta = qr(1,i-1) + aux2(1,i-1)
            topo1 = aux1(1,i-1)
            topo3 = aux3(1,i-1)
c            s1 = vl-sqrt(g*hl)
c            s3 = vl+sqrt(g*hl)
c            s2 = 0.5d0*(s1+s3)
       else
            eta = ql(1,i) + aux2(1,i)
            topo1 = aux1(1,i)
            topo3 = aux3(1,i)
c            s1 = vr-sqrt(g*hr)
c            s3 = vr+sqrt(g*hr)
c            s2 = 0.5d0*(s1+s3)
       endif
       if (eta.lt.max(topo1,topo3)) go to 90

      if (icoordsys.eq.2) then
         if (ixy.eq.2) then
            dxdcp=(Rearth*pi/180.d0)
            dxdcm = dxdcp
         else
            if (ilr.eq.1) then
               dxdcp = Rearth*pi*cos(aux3(3,i-1))/180.d0
               dxdcm = Rearth*pi*cos(aux1(3,i-1))/180.d0
            else
               dxdcp = Rearth*pi*cos(aux3(3,i))/180.d0
               dxdcm = Rearth*pi*cos(aux1(3,i))/180.d0
            endif
         endif
      endif

c=====Determine some speeds necessary for the Jacobian=================
            vhat=(vr*dsqrt(hr))/(dsqrt(hr)+dsqrt(hl)) +
     &        (vl*dsqrt(hl))/(dsqrt(hr)+dsqrt(hl))

            uhat=(ur*dsqrt(hr))/(dsqrt(hr)+dsqrt(hl)) +
     &        (ul*dsqrt(hl))/(dsqrt(hr)+dsqrt(hl))
            hhat=(hr+hl)/2.d0

            roe1=vhat-dsqrt(g*hhat)
            roe3=vhat+dsqrt(g*hhat)

            s1l=vl-dsqrt(g*hl)
            s3r=vr+dsqrt(g*hr)

            s1=dmin1(roe1,s1l)
            s3=dmax1(roe3,s3r)

            s2=0.5d0*(s1+s3)

           s(1)=s1
           s(2)=s2
           s(3)=s3
c=======================Determine asdq decomposition (beta)============
         delf1=asdq(1,i)
         delf2=asdq(mu,i)
         delf3=asdq(mv, i)

         beta(1) = (s3*delf1/(s3-s1))-(delf3/(s3-s1))
         beta(2) = -s2*delf1 + delf2
         beta(3) = (delf3/(s3-s1))-(s1*delf1/(s3-s1))
c======================End =================================================

c=====================Set-up eigenvectors===================================
         r(1,1) = 1.d0
         r(2,1) = s2
         r(3,1) = s1

         r(1,2) = 0.d0
         r(2,2) = 1.d0
         r(3,2) = 0.d0

         r(1,3) = 1.d0
         r(2,3) = s2
         r(3,3) = s3
c============================================================================
90      continue
c============= compute fluctuations==========================================

            do  m=1,meqn
               bmasdq(m,i)=0.0d0
               bpasdq(m,i)=0.0d0
            enddo
            do  mw=1,3
               if (s(mw).lt.0.d0) then
                 bmasdq(1,i) =bmasdq(1,i) + dxdcm*s(mw)*beta(mw)*r(1,mw)
                 bmasdq(mu,i)=bmasdq(mu,i)+ dxdcm*s(mw)*beta(mw)*r(2,mw)
                 bmasdq(mv,i)=bmasdq(mv,i)+ dxdcm*s(mw)*beta(mw)*r(3,mw)
               elseif (s(mw).gt.0.d0) then
                 bpasdq(1,i) =bpasdq(1,i) + dxdcp*s(mw)*beta(mw)*r(1,mw)
                 bpasdq(mu,i)=bpasdq(mu,i)+ dxdcp*s(mw)*beta(mw)*r(2,mw)
                 bpasdq(mv,i)=bpasdq(mv,i)+ dxdcp*s(mw)*beta(mw)*r(3,mw)
               endif
            enddo
c========================================================================
         enddo
c

c
      return
      end
