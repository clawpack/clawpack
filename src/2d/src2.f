c
c
c =========================================================
      subroutine src2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,
     &             q,maux,aux,t,dt)
c =========================================================
      use geoclaw_module
      
      implicit double precision (a-h,o-z)

      dimension   q(meqn,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc)
      dimension aux(maux,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc)
c

c     # incorporates friction using Manning coefficient

      if (coeffmanning.eq.0.d0 .or. frictiondepth.eq.0) return

c     # check for NANs in solution:
c     call check4nans(maxmx,maxmy,meqn,mbc,mx,my,q,t,2)

      g=grav
      coeff = coeffmanning
      tol = 1.d-30  !# to prevent divide by zero in gamma

*     ! friction--------------------------------------------------------
      if (coeffmanning.gt.0.d0.and.frictiondepth.gt.0.d0) then
         do j=1,my
            do i=1,mx

               h=q(1,i,j)
               if (h.le.frictiondepth) then
c                 # apply friction source term only in shallower water
                  hu=q(2,i,j)
                  hv=q(3,i,j)

                  if (h.lt.tol) then
                     q(2,i,j)=0.d0
                     q(3,i,j)=0.d0
                  else
                     gamma= dsqrt(hu**2 + hv**2)*(g*coeff**2)/(h**(7/3))
                     dgamma=1.d0 + dt*gamma
                     q(2,i,j)= q(2,i,j)/dgamma
                     q(3,i,j)= q(3,i,j)/dgamma
                  endif
               endif
            enddo
         enddo
      endif
*     ! ----------------------------------------------------------------

*     ! coriolis--------------------------------------------------------
      if (icoordsys.eq.2.and.icoriolis.eq.1) then
         w = 2.d0*pi/(86400.d0) !angular velocity of earth
         do i=1,mx
            do j=1,my
               ycell = ylower + (j-.5d0)*dy
               cor = 2.d0*w*sin(pi*ycell/180.d0)
               ct = cor*dt
*              !integrate momentum exactly using matrix exponential
*              !forth order term should be sufficient since cor^3 ~= eps
               hu0 = q(2,i,j)
               hv0 = q(3,i,j)
*              !dq/dt = 2w*sin(latitude)*[0 1 ; -1 0] q = Aq
*              !e^Adt = [a11 a12; a21 a22] + I
               a11 = -0.5d0*ct**2 + ct**4/24.d0
               a12 = ct - ct**3/6.0d0
               a21 = -ct + ct**3/6.0d0
               a22 = a11
*              !q = e^Adt * q0
               q(2,i,j) = q(2,i,j) + hu0*a11 + hv0*a12
               q(3,i,j) = q(3,i,j) + hu0*a21 + hv0*a22
            enddo
         enddo
      endif
*     ! ----------------------------------------------------------------
c
      return
      end
