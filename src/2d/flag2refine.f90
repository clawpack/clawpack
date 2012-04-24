! ::::::::::::::::::::: flag2refine ::::::::::::::::::::::::::::::::::
!
! User routine to control flagging of points for refinement.
!
! Specific for GeoClaw for tsunami applications and related problems
!
!
! The logical function allowflag(x,y,t) is called to
! check whether further refinement at this level is allowed in this cell
! at this time.
!
!    q   = grid values including ghost cells (bndry vals at specified
!          time have already been set, so can use ghost cell values too)
!
!  aux   = aux array on this grid patch
!
! amrflags  = array to be flagged with either the value
!             DONTFLAG (no refinement needed)  or
!             DOFLAG   (refinement desired)
!
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine flag2refine(mx,my,mbc,meqn,maux,xlower,ylower,dx,dy,t,level,tolsp,q,aux,amrflags,DONTFLAG,DOFLAG)

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
    real(kind=8) :: x_c,y_c,x_low,y_low,x_hi,y_hi
    real(kind=8) :: h,speed,surface

    ! Loop over interior points on this grid:
    ! (i,j) grid cell is [x_low,x_hi] x [y_low,y_hi]
    do j=1,my
        y_c= ylower +  (j-0.5d0)*dy
        y_low = ylower + (j-1)*dy
        y_hi = ylower + j*dy
        do i = 1,mx
            x_c = xlower +  (i-0.5d0)*dx
            x_low = xlower +  (i-1)*dx
            x_hi = xlower +  i*dx

            ! Default is to not flag
            amrflags(i,j) = DONTFLAG

            ! Check to see if refinement is forced in any topography file region:
            do m=1,mtopofiles
                if (level < minleveltopo(m) .and. t >= tlowtopo(m) .and. t <= thitopo(m)) then
                    if (  x_hi > xlowtopo(m) .and. x_low < xhitopo(m) .and. &
                          y_hi > ylowtopo(m) .and. y_low < yhitopo(m)) then
                        
                        amrflags(i,j) = DOFLAG
                        cycle
                    endif
                endif
            enddo

            ! Check to see if refinement is forced in any other region:
            do m=1,num_regions
                if (level < min_level_region(m) .and. t >= t_low_region(m) .and. t <= t_hi_region(m)) then
                    if (x_hi > x_low_region(m) .and. x_low < x_hi_region(m) .and. &
                        y_hi > y_low_region(m) .and. y_low < y_hi_region(m)) then
                    
                        amrflags(i,j) = DOFLAG
                        cycle
                    endif
                endif
            enddo

            ! Check if we're in the dtopo region and need to refine:
            ! force refinement to level minleveldtopo
            do m = 1,num_dtopo
                if (level < minleveldtopo(m).and. &
                    t <= tfdtopo(m) .and. & !t.ge.t0dtopo(m).and.
                    x_hi > xlowdtopo(m) .and. x_low < xhidtopo(m).and. &
                    y_hi > ylowdtopo(m) .and. y_low < yhidtopo(m)) then
                    
                    amrflags(i,j) = DOFLAG    
                    cycle
                endif
            enddo

            ! Check if we're in the region where initial perturbation is
            ! specified and need to force refinement:
            ! This assumes that t0 = 0.d0
            if (qinit_type > 0 .and. t == 0.d0) then 
                if (level < min_level_qinit .and. & 
                    x_hi > x_low_qinit .and. x_low < x_hi_qinit .and. &
                    y_hi > y_low_qinit .and. y_low < y_hi_qinit) then
                    
                    amrflags(i,j) = DOFLAG
                    cycle
                endif
            endif

            ! -----------------------------------------------------------------
            ! Refinement not forced, so check if it is allowed, and if so,
            ! check if there is a reason to flag this point:
            if (allowflag(x_c,y_c,t,level)) then
                
                ! Calculate each layer's surface and speed
                h = q(3*(layers-1)+1,i,j) / rho(layers)
                
                if (h > drytolerance) then
                    ! We are in a wet cell, compute surface and speed of this layer
                    surface = h + aux(1,i,j)
                    speed = sqrt(q(3*(layers-1)+2,i,j)**2 + &
                                    q(3*(layers-1)+3,i,j)**2) / q(3*(layers-1)+1,i,j)
                    
                    ! Check to see if the surface is close to the intitial
                    ! value and perhaps refine if it is not                     
                    if (abs(surface - eta_init(layers)) > wavetolerance) then
                        ! Check to see if we are near shore
                        if (h < depthdeep) then
                            amrflags(i,j) = DOFLAG
                            cycle ! Move onto next point
                        ! If we are not, limit the level of refinement
                        else if (level < maxleveldeep) then
                            amrflags(i,j) = DOFLAG
                            cycle ! Move onto next point
                        endif
                    endif
                    
                    ! Check speed tolerances
                    ! TODO: The cycle command unfortunately does not work here
                    ! as it would only kick out of the inner loop, need to add
                    ! more logic to handle this
                    do m=1,max_speed_nest
                        if (speed > speed_tolerance(m) .and. level <= m) then
                            amrflags(i,j) = DOFLAG
                        endif
                    enddo
                endif

                ! Do for the rest of the layers
!                 do m=layers-1,1,-1
!                     h = q(3*(m-1)+1,i,j) / rho(m)
!                     if (h > drytolerance) then
!                         surface(m) = h + surface(m-1)
!                         speed(m) = sqrt(q(3*(m-1)+2,i,j)**2 + &
!                                              q(3*(m-1)+3,i,j)**2) / q(3*(m-1)+1,i,j)
!                         if (h < shore_tolerance) then
!                             cycle
!                         else
!                             
!                         endif
!                         
!                     else
!                         surface(m) = aux(1,i,j)
!                         speed(m) = 0.d0
!                         shore_region(m) = .true.
!                     endif
!                 enddo
!                 
!                 ! Check if each surface is less than the wave tolerance specified
!                 if (q(1,i,j) / rho(1) > drytolerance) then
!                     abs(surface(1) - eta_init(1)) > wave_tolerance(1)
                
                
!c
!c            # RJL: not sure why this next line is done?
!c            # Need to fix for arb.  sealevel?
!c            surface = dsign(surface,q(1,i,j))
!
!c            # DLG: it was a way to prevent refining on dry land...
!c            # probably should be changed if we allow arbitrary sealevel
!c            # by adding sealevel to surface or something.
!
!c            # determine region type and refinement params
!
!!              shoreregion = dabs(aux(1,i,j)) .lt. depthdeep
!              wave = (dabs(surface-eta_init(1)).gt.wavetolerance.and.
!      &                q(1,i,j).gt.drytolerance)
! c             #DLG: changing following: didn't work so well for non-lat-lon grids
! c              shoretol = depthdeep*(dx*dy)
!                shoretol = depthdeep
! c
! 
!              if (wave) then
! c               # the surface is not at sea level here
! 
!                 if (level.lt.maxleveldeep) then
! c                   # in deep water we can refine to this level
!                     amrflags(i,j)=DOFLAG
!                     go to 100 !# flagged, so no need to check anything else
!                     endif
! 
! c                if (shoreregion) then
! c                  shoreline=.false.
! c                 # check if any neighboring cell is dry:
! c                  do jj=-1,1
! c                   do ii=-1,1
! c                    shoreline = shoreline.or.q(1,i+ii,j+jj).le.shoretol
! c                   enddo
! c                  enddo
! 
! c                 shoreline=shoreline.and.q(1,i,j).gt.shoretol
!                  shoreline = shoreregion
! 
!                  if (shoreline.and.q(1,i,j).gt.drytolerance) then
! c                    # following comment is regarding commented nested do loop above.
! c                    # this cell is wet and a neighbor is dry ==> shoreline
!                      amrflags(i,j)=DOFLAG
!                      go to 100 !# flagged, so no need to check anything else
!                  endif
! 
! c               endif
!              endif
! 
!              
! 
            endif  ! Allow flag check
        enddo ! x-loop
    enddo ! y-loop
end subroutine flag2refine
