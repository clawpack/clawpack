c
c
c     ===============================================================
      subroutine rpt2_ptwise(ixy, ilr, meqn, maux, qLeft, qRight,
     &           auxLeft, auxRight, asdq, bmasdq, bpasdq)
c     ===============================================================
c
c     # Riemann solver in the transverse direction using an einfeldt
c     # Jacobian.
c
c     # Split asdq (= A^* \Delta q, where *=+ if ilr=1 or *=- if ilr=2)
c     # into down-going flux difference bmasdq (= B^- A^* \Delta q)
c     #    and up-going flux difference bpasdq (= B^+ A^* \Delta q)
c
c     # This is a new-style Riemann solver that works pointwise, with
c     # Left state qLeft and right state qRight.  (Note change in notation!!)
c
c     # Often only asdq is needed, not q or aux arrays, but these are
c     # provided in a more general form than in the previous rpt2 format:
c     #   qLeft(0, 1:m)  is q in left state along this grid row,
c     #   qLeft(-1,1:m)  is q in left state along row below,
c     #   qLeft(1, 1:m)  is q in left state along row above,
c     # Simlarly for qRight, auxLeft, and auxRight. 
c
c
      use geoclaw_module
      implicit none
      integer ixy, ilr, meqn, maux
      double precision qLeft(-1:1,meqn), qRight(-1:1,meqn) 
      double precision auxLeft(-1:1,maux), auxRight(-1:1,maux) 
      double precision asdq(meqn), bpasdq(meqn), bmasdq(meqn)
c

c-----------------------last modified 1/10/05----------------------


      double precision  s(3)
      double precision  r(3,3)
      double precision  beta(3)
      double precision  g,tol,abs_tol
      double precision  hl,hr,hul,hur,hvl,hvr,vl,vr,ul,ur,bl,br
      double precision  uhat,vhat,hhat,roe1,roe3,s1,s2,s3,s1l,s3r
      double precision  delf1,delf2,delf3,dxdcd,dxdcu

      integer i,m,mw,mu,mv

      g=grav
      tol=drytolerance
      abs_tol=drytolerance/100.d0

      if (ixy.eq.1) then
	  mu = 2
	  mv = 3
      else
	  mu = 3
	  mv = 2
      endif



         hl=qLeft(0,1)
         hr=qRight(0,1)
         hul=qLeft(0,mu)
         hur=qRight(0,mu)
         hvl=qLeft(0,mv)
         hvr=qRight(0,mv)

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
       
       do mw=1,3
          s(mw)=0.d0
          beta(mw)=0.d0
          do m=1,meqn
             r(m,mw)=0.d0
          enddo
       enddo

       if (hl.le.0.d0.and.hr.le.0.d0) go to 90

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
         delf1=asdq(1)
         delf2=asdq(mu)
         delf3=asdq(mv)

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
               bmasdq(m)=0.0d0
               bpasdq(m)=0.0d0
            enddo
            do  mw=1,3
               if (s(mw).lt.0.d0) then
                     bmasdq(1) =bmasdq(1)  + s(mw)*beta(mw)*r(1,mw)
                     bmasdq(mu)=bmasdq(mu) + s(mw)*beta(mw)*r(2,mw)
                     bmasdq(mv)=bmasdq(mv) + s(mw)*beta(mw)*r(3,mw)
               elseif (s(mw).gt.0.d0) then
                     bpasdq(1) =bpasdq(1)  + s(mw)*beta(mw)*r(1,mw)
                     bpasdq(mu)=bpasdq(mu) + s(mw)*beta(mw)*r(2,mw)
                     bpasdq(mv)=bpasdq(mv) + s(mw)*beta(mw)*r(3,mw)
               endif
            enddo
c
      return
      end
