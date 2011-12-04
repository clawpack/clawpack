c
c
c =========================================================
      subroutine src1d(meqn,mbc,mx1d,q1d,maux,aux1d,t,dt)
c =========================================================
      use geoclaw_module

      implicit double precision (a-h,o-z)

      dimension   q1d(meqn, mx1d)
      dimension aux1d(maux, mx1d)
c

c
c     # This routine should be a simplified version of src2
c     # which applies source terms for a 1-d slice of data along the
c     # edge of a grid.  This is called only from qad where the conservative
c     # fix-up is applied and is used to apply source terms over partial
c     # time steps to the coarse grid cell values used in solving Riemann
c     # problems at the interface between coarse and fine grids.


c     # incorporates friction using Manning coefficient
        g=grav
        coeff = coeffmanning
        tol = 1.d-30  !# to prevent divide by zero in gamma

      if (coeffmanning.gt.0.d0.and.frictiondepth.gt.0.d0) then
        do i=1,mx1d
           h=q1d(1,i)
           if (h.lt.frictiondepth) then
c            # apply friction source term only in shallower water
             hu=q1d(2,i)
             hv=q1d(3,i)

             if (h.lt.tol) then
                q1d(2,i)=0.d0
                q1d(3,i)=0.d0
             else
                gamma= dsqrt(hu**2 + hv**2)*(g*coeff**2)/(h**(7/3))
                dgamma=1.d0 + dt*gamma
                q1d(2,i)= q1d(2,i)/dgamma
                q1d(3,i)= q1d(3,i)/dgamma
             endif
           endif
        enddo
      endif

*     ! coriolis--------------------------------------------------------
*     ! coriolis--------------------------------------------------------
      if (icoordsys.eq.2.and.icoriolis.eq.1) then
         w = 2.d0*pi/(86400.d0) !angular velocity of earth
         do i=1,mx1d
            cor = 2.d0*w*sin(aux1d(3,i))
            ct = cor*dt
*           !integrate momentum exactly using matrix exponential
*           !forth order term should be sufficient since cor^3 ~= eps
            hu0 = q1d(2,i)
            hv0 = q1d(3,i)
*           !dq/dt = 2w*sin(latitude)*[0 1 ; -1 0] q = Aq
*           !e^Adt = [a11 a12; a21 a22] + I
            a11 = -0.5d0*ct**2 + ct**4/24.d0
            a12 = ct - ct**3/6.0d0
            a21 = -ct + ct**3/6.0d0
            a22 = a11
*           !q = e^Adt * q0
            q1d(2,i) = q1d(2,i) + hu0*a11 + hv0*a12
            q1d(3,i) = q1d(3,i) + hu0*a21 + hv0*a22
            enddo
         endif
*     ! ----------------------------------------------------------------
c
      return
      end
