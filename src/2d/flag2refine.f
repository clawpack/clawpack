c
c -------------------------------------------------------------------
      subroutine flag2refine(mx,my,mbc,meqn,maux,xlower,ylower,dx,dy,
     &                 t,level,tolsp,q,aux,amrflags,DONTFLAG,DOFLAG)
c -------------------------------------------------------------------

c
c ::::::::::::::::::::: flag2refine ::::::::::::::::::::::::::::::::::
c
c User routine to control flagging of points for refinement.
c
c Specific for GeoClaw for tsunami applications and related problems
c
c
c The logical function allowflag(x,y,t) is called to
c check whether further refinement at this level is allowed in this cell
c at this time.
c
c    q   = grid values including ghost cells (bndry vals at specified
c          time have already been set, so can use ghost cell values too)
c
c  aux   = aux array on this grid patch
c
c amrflags  = array to be flagged with either the value
c             DONTFLAG (no refinement needed)  or
c             DOFLAG   (refinement desired)
c

c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      use geoclaw_module
      use topo_module
      use dtopo_module
      use regions_module
      use qinit_module

      implicit none

      ! Subroutine arguments
      integer, intent(in) :: mx,my,mbc,meqn,maux,level
      real(kind=8), intent(in) :: xlower,ylower,dx,dy,t,tolsp
      
      real(kind=8), intent(in) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      
      ! Flagging
      real(kind=8),intent(inout) :: amrflags(1-mbc:mx+mbc,1-mbc:my+mbc)
      real(kind=8), intent(in) :: DONTFLAG
      real(kind=8), intent(in) :: DOFLAG
      
      logical :: allowflag
      external allowflag
      
      ! Generic locals
      integer :: i,j,m
      integer :: minleveldt,maxleveldt
      logical :: shoreregion,wave,shoreline
      real(kind=8) :: x,x1,x2,xlow,xhi,y,y1,y2,ylow,yhi
      real(kind=8) :: t0dt,tfdt
      real(kind=8) :: surface,shoretol,speed

c     # loop over interior points on this grid:

      do 200 j = 1,my
        y = ylower +  (j-0.5d0)*dy
        y1 = ylower + (j-1)*dy
        y2 = ylower + j*dy
        do 100 i = 1,mx
          x = xlower +  (i-0.5d0)*dx
          x1 = xlower +  (i-1)*dx
          x2 = xlower +  i*dx
c         # (i,j) grid cell is [x1,x2] x [y1,y2].

c         # default for each point is not to flag unless some condition
c         # below is satisfied:

          amrflags(i,j) = DONTFLAG


          do 30 m=1,mtopofiles
c           # check to see if refinement is forced in any topo file region:
            if (level < minleveltopo(m) .and.
     &          t >= tlowtopo(m) .and. t <= thitopo(m)) then
              xlow = xlowtopo(m)
              xhi = xhitopo(m)
              ylow = ylowtopo(m)
              yhi = yhitopo(m)
                 if (x2.gt.xlow.and.x1.lt.xhi.and.
     &               y2.gt.ylow.and.y1.lt.yhi) then
                    amrflags(i,j) = DOFLAG
                    go to 100 !# flagged, so no need to check anything else
                    endif
              endif
   30       continue

           do m=1,num_regions
c           # check to see if refinement is forced in any other region:
            if (level < min_level_region(m) .and.
     &          t >= t_low_region(m) .and. t <= t_hi_region(m)) then
              xlow = x_low_region(m)
              xhi  = x_hi_region(m)
              ylow = y_low_region(m)
              yhi  = y_hi_region(m)
                 if (x2 > xlow .and. x1 < xhi .and.
     &               y2 > ylow .and. y1 < yhi) then
                    amrflags(i,j) = DOFLAG
                    go to 100 !# flagged, so no need to check anything else
                    endif
              endif
            enddo

         do m = 1,num_dtopo
c           # check if we're in the dtopo region and need to refine:
c           # force refinement to level minleveldtopo
            t0dt = t0dtopo(m)
            tfdt = tfdtopo(m)
            minleveldt = minleveldtopo(m)
            if (level.lt.minleveldtopo(m).and.
     &              t.le.tfdtopo(m).and. !t.ge.t0dtopo(m).and.
     &              x2.gt.xlowdtopo(m).and.x1.lt.xhidtopo(m).and.
     &              y2.gt.ylowdtopo(m).and.y1.lt.yhidtopo(m)) then
                amrflags(i,j)=DOFLAG
                go to 100 !# flagged, so no need to check anything else
                endif
         enddo


         if (qinit_type > 0 .and. t == 0.d0) then
c           # check if we're in the region where initial perturbation is
c           # specified and need to force refinement:
            if (level < min_level_qinit.and.
     &              x2 > x_low_qinit .and. x1 < x_hi_qinit.and.
     &              y2 > y_low_qinit .and. y1 < y_hi_qinit) then
                amrflags(i,j)=DOFLAG
                go to 100 !# flagged, so no need to check anything else
                endif
             endif

c        -----------------------------------------------------------------

c        # refinement not forced, so check if it is allowed, and if so,
c        # check if there is a reason to flag this point:

         if (allowflag(x,y,t,level)) then

             surface = q(1,i,j) + aux(1,i,j)
c
c            # RJL: not sure why this next line is done?
c            # Need to fix for arb.  sealevel?
c            surface = dsign(surface,q(1,i,j))

c            # DLG: it was a way to prevent refining on dry land...
c            # probably should be changed if we allow arbitrary sealevel
c            # by adding sealevel to surface or something.

c            # determine region type and refinement params

             shoreregion = dabs(aux(1,i,j)) .lt. depthdeep
             wave = (dabs(surface-eta_init(1)).gt.wavetolerance.and.
     &                q(1,i,j).gt.drytolerance)
c             #DLG: changing following: didn't work so well for non-lat-lon grids
c              shoretol = depthdeep*(dx*dy)
               shoretol = depthdeep
c

             if (wave) then
c               # the surface is not at sea level here

                if (level.lt.maxleveldeep) then
c                   # in deep water we can refine to this level
                    amrflags(i,j)=DOFLAG
                    go to 100 !# flagged, so no need to check anything else
                    endif

c                if (shoreregion) then
c                  shoreline=.false.
c                 # check if any neighboring cell is dry:
c                  do jj=-1,1
c                   do ii=-1,1
c                    shoreline = shoreline.or.q(1,i+ii,j+jj).le.shoretol
c                   enddo
c                  enddo

c                 shoreline=shoreline.and.q(1,i,j).gt.shoretol
                 shoreline = shoreregion

                 if (shoreline.and.q(1,i,j).gt.drytolerance) then
c                    # following comment is regarding commented nested do loop above.
c                    # this cell is wet and a neighbor is dry ==> shoreline
                     amrflags(i,j)=DOFLAG
                     go to 100 !# flagged, so no need to check anything else
                 endif

c               endif
             endif

          endif
          
c        Refine based on momentum or speed of water
c        Had to remove this form the allow flag block as it checks for t > 0
c        and was not allowing refinement before t = 0, we need this as the 
c        storm surge has ramp up time that may need refinement (KTM 2010-8-4)
C           speed = sqrt(q(2,i,j)**2 + q(3,i,j)**2)
C           ! This is only important for layers > 1, otherwise this is already speed
C           if (.not.momentum_refinement) then
C               if (q(1,i,j) > drytolerance) then
C                   speed = speed / q(1,i,j)
C               endif
C           endif
C           do m=1,max_speed_nest
C               if ((speed > speed_refine(m)).and.(level <= m)) then
C                   amrflags(i,j) = DOFLAG
C                   go to 100
C               endif
C           enddo
          

 100     continue  !# end loop on i
 200    continue   !# end loop on j

      return
      end
