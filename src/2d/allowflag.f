c     =========================================
      logical function allowflag(x,y,t,level)
c     =========================================

c     # Indicate whether the grid point at (x,y,t) at this refinement level
c     # is allowed to be flagged for further refinement.
c
c     # Modified for GeoClaw to check whether the point lies in any of
c     # the various regions specified in the data files.
c
c     # This routine is called from routine flag2refine.
c
c     # If Richardson error estimates are used (if tol>0) then this routine
c     # is also called from errf1.

      use geoclaw_module
      use topo_module
      use dtopo_module
      use regions_module
      use qinit_module

      implicit none

C     Function arguments
      real(kind=8), intent(in) :: x,y,t
      integer, intent(in) :: level
      
C     Locals
      integer :: m

      allowflag=.false.

c      following commented by dlg on 10/9/08.
c      my understanding of maxleveldeep might be differnet
c      still shouldn't be allowed if maxlevel allowed in a region is less
c      than maxleveldeep
c      might want to allow high levels of refinement in some deep regions
c      but not others.
c
c      if (level .lt. maxleveldeep) then
c         # allow refinement to next level in deep water
c          allowflag = .true.
c          go to 900  !# no need to check anything else
c          endif

      do m=1,mtopofiles
     	if (level < maxleveltopo(m)) then
      	  if (x > xlowtopo(m) .and. x < xhitopo(m) .and.
     &	      y > ylowtopo(m) .and .y < yhitopo(m) .and.
     &	      t > tlowtopo(m) .and .t < thitopo(m)) then
     		  allowflag=.true.
                  go to 900  !# no need to check anything else
	  endif
	endif
      enddo

      do m=1,num_regions
     	if (level < max_level_region(m)) then
      	  if (x.gt.x_low_region(m) .and. x <  x_hi_region(m).and.
     &	      y.gt.y_low_region(m) .and. y <  y_hi_region(m).and.
     &	      t.ge.t_low_region(m) .and. t <= t_hi_region(m)) then
     		  allowflag=.true.
                  go to 900  !# no need to check anything else
	  endif
	endif
       enddo

      do m=1,num_dtopo
        if (x.gt.xlowdtopo(m).and.x.lt.xhidtopo(m).and.
     &       y.gt.ylowdtopo(m).and.y.lt.yhidtopo(m).and.
     &       t.ge.t0dtopo(m).and.t.le.tfdtopo(m)) then
             if (level.lt.maxleveldtopo(m)) then
                  allowflag=.true.
                  go to 900  !# no need to check anything else
                  endif
        endif
      enddo

      if (t == 0.d0 .and. qinit_type > 0) then
        if (x > x_low_qinit .and. x < x_hi_qinit .and.
     &	    y > y_low_qinit .and. y < y_hi_qinit) then
     		if (level.lt.max_level_qinit) then
                  allowflag=.true.
                  go to 900  !# no need to check anything else
                  endif
        endif
      endif

  900 continue
      return
      end
