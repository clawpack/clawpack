! ========================================
!  logical function allowflag(x,y,t)
!
!  Indicate whether the grid point at (x,y,t) at this refinement level
!  is allowed to be flagged for further refinement.
!
!  Modified for GeoClaw to check whether the point lies in any of
!  the various regions specified in the data files.
!
!  This routine is called from routine flag2refine.
!
!  If Richardson error estimates are used (if tol>0) then this routine
!  is also called from errf1.
!
logical function allowflag(x,y,t,level)

    use amr_module, only: t0
    use geoclaw_module
    use regions_module
    use refinement_module
    use topo_module
    use qinit_module

    implicit none

    ! Function arguments
    real(kind=8), intent(in) :: x,y,t
    integer, intent(in) :: level

    ! Locals
    integer :: m

    allowflag = .false.

!   following commented by dlg on 10/9/08.
!   my understanding of maxleveldeep might be differnet
!   still shouldn't be allowed if maxlevel allowed in a region is less
!   than maxleveldeep
!   might want to allow high levels of refinement in some deep regions
!   but not others.
!
!   if (level .lt. maxleveldeep) then
!      # allow refinement to next level in deep water
!       allowflag = .true.
!       go to 900  !# no need to check anything else
!       endif

    ! Allow flagging everywhere if using test bathymetry
    if (test_topography > 1) then
        allowflag = .true.
        return
    endif
    do m=1,mtopofiles
        if (level < maxleveltopo(m)) then
            if (x > xlowtopo(m) .and. x < xhitopo(m) .and. &
                y > ylowtopo(m) .and. y < yhitopo(m) .and. &
                t >= tlowtopo(m) .and. t < thitopo(m)) then

                allowflag = .true.
                return
            endif
        endif
    enddo
    do m=1,num_regions
        if (level < regions(m)%max_level) then
            if (x > regions(m)%x_low .and. x <  regions(m)%x_hi.and. &
                y > regions(m)%y_low .and. y <  regions(m)%y_hi.and. &
                t >= regions(m)%t_low .and. t <= regions(m)%t_hi) then

                allowflag = .true.
                return
            endif
        endif
    enddo

    do m=1,num_dtopo
        if (x >  xlowdtopo(m) .and. x < xhidtopo(m).and. &
            y >  ylowdtopo(m) .and. y < yhidtopo(m).and. &
            t >= t0dtopo(m)   .and. t <= tfdtopo(m)) then

            if (level < maxleveldtopo(m)) then
                allowflag = .true.
                return
            endif
        endif
    enddo

    if (t == t0 .and. qinit_type > 0) then
        if (x > x_low_qinit .and. x < x_hi_qinit .and. &
            y > y_low_qinit .and. y < y_hi_qinit) then

            if (level < max_level_qinit) then
                allowflag = .true.
                return
            endif
        endif
    endif

end function allowflag
