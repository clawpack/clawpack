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
subroutine flag2refine2(mx,my,mbc,mbuff,meqn,maux,xlower,ylower,dx,dy,t,level, &
                       tolsp,q,aux,amrflags)

    use amr_module, only: mxnest, t0, DONTFLAG, DOFLAG, UNSET
    use geoclaw_module, only: sea_level, rho
    use geoclaw_module, only: spherical_distance, coordinate_system

    use topo_module, only: tlowtopo,thitopo,xlowtopo,xhitopo,ylowtopo,yhitopo
    use topo_module, only: minleveltopo,mtopofiles

    use topo_module, only: tfdtopo,xlowdtopo,xhidtopo,ylowdtopo,yhidtopo
    use topo_module, only: minleveldtopo,num_dtopo

    use qinit_module, only: x_low_qinit,x_hi_qinit,y_low_qinit,y_hi_qinit
    use qinit_module, only: min_level_qinit,qinit_type

    use storm_module, only: storm_specification_type, wind_refine, R_refine
    use storm_module, only: storm_location, wind_forcing, wind_index, wind_refine

    use regions_module, only: num_regions, regions
    use refinement_module, only: speed_tolerance, deep_depth, max_level_deep

    use multilayer_module, only: num_layers, eta_init, dry_tolerance, wave_tolerance

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: mx,my,mbc,meqn,maux,level,mbuff
    real(kind=8), intent(in) :: xlower,ylower,dx,dy,t,tolsp

    real(kind=8), intent(in) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    ! Flagging
    real(kind=8), intent(in out) :: amrflags(1-mbuff:mx+mbuff,1-mbuff:my+mbuff)

    logical :: allowflag
    external allowflag

    ! Generic locals
    integer :: i, j, m, layer
    real(kind=8) :: x_c,y_c,x_low,y_low,x_hi,y_hi
    real(kind=8) :: speed, eta, ds, b

    ! Storm specific variables
    real(kind=8) :: R_eye(2), wind_speed

    ! Don't initialize flags, since they were already
    ! flagged by flagregions2
    ! amrflags = DONTFLAG

    ! Loop over interior points on this grid
    ! (i,j) grid cell is [x_low,x_hi] x [y_low,y_hi], cell center at (x_c,y_c)
    y_loop: do j=1,my
        y_low = ylower + (j - 1) * dy
        y_c = ylower + (j - 0.5d0) * dy
        y_hi = ylower + j * dy

        x_loop: do i = 1,mx
            x_low = xlower + (i - 1) * dx
            x_c = xlower + (i - 0.5d0) * dx
            x_hi = xlower + i * dx

            ! The following conditions are only checked in the horizontal and
            ! override the allowflag routine

            ! ************* Storm Based Refinement ****************
            ! Check to see if we are some specified distance from the eye of
            ! the storm and refine if we are
            if (storm_specification_type /= 0) then
                R_eye = storm_location(t)
                do m=1,size(R_refine,1)
                    if (coordinate_system == 2) then
                        ds = spherical_distance(x_c, y_c, R_eye(1), R_eye(2))
                    else
                        ds = sqrt((x_c - R_eye(1))**2 + (y_c - R_eye(2))**2)
                    end if
                    
                    if ( ds < R_refine(m) .and. level <= m .and. amrflags(i,j) == UNSET) then
                        amrflags(i,j) = DOFLAG
                        cycle x_loop
                    endif
                enddo
                
                ! Refine based on wind speed
                if (wind_forcing) then
                    wind_speed = sqrt(aux(wind_index,i,j)**2 + aux(wind_index+1,i,j)**2)
                    do m=1,size(wind_refine,1)
                        if ((wind_speed > wind_refine(m)) .and. (level <= m) &
                              .and. amrflags(i,j) == UNSET) then
                            amrflags(i,j) = DOFLAG
                            cycle x_loop
                        endif
                    enddo
                endif
            endif
            ! *****************************************************

            ! Check to see if refinement is forced in any topography file region:
            do m=1,mtopofiles
                if (level < minleveltopo(m) .and. t >= tlowtopo(m) .and. t <= thitopo(m)) then
                    if (  x_hi > xlowtopo(m) .and. x_low < xhitopo(m) .and. &
                          y_hi > ylowtopo(m) .and. y_low < yhitopo(m) &
                           .and. amrflags(i,j) == UNSET) then

                        amrflags(i,j) = DOFLAG
                        cycle x_loop
                    endif
                endif
            enddo

            ! Check if we're in the dtopo region and need to refine:
            ! force refinement to level minleveldtopo
            do m = 1,num_dtopo
                if (level < minleveldtopo(m).and. &
                    t <= tfdtopo(m) .and. & !t.ge.t0dtopo(m).and.
                    x_hi > xlowdtopo(m) .and. x_low < xhidtopo(m).and. &
                    y_hi > ylowdtopo(m) .and. y_low < yhidtopo(m) &
                     .and. amrflags(i,j) == UNSET) then

                    amrflags(i,j) = DOFLAG
                    cycle x_loop
                endif
            enddo

            ! Check if we're in the region where initial perturbation is
            ! specified and need to force refinement:
            ! This assumes that t0 = 0.d0, should really be t0 but we do
            ! not have access to that parameter in this routine
            if (qinit_type > 0 .and. t == t0 .and. amrflags(i,j) == UNSET) then
                if (level < min_level_qinit .and. &
                    x_hi > x_low_qinit .and. x_low < x_hi_qinit .and. &
                    y_hi > y_low_qinit .and. y_low < y_hi_qinit) then

                    amrflags(i,j) = DOFLAG
                    cycle x_loop
                endif
            endif

            ! -----------------------------------------------------------------
            ! Refinement not forced, so check if it is allowed and if so,
            ! check if there is a reason to flag this point:
            if (allowflag(x_c,y_c,t,level) .and. amrflags(i,j) == UNSET) then
                b = aux(1,i,j)
                do layer = num_layers, 1, -1
                    if (q(3*layer-2,i,j) / rho(layer) > dry_tolerance(layer)) then

                        eta = q(3*layer-2,i,j) / rho(layer) + b

                        ! Check wave criteria
                        if (abs(eta - eta_init(layer)) > wave_tolerance(layer)) then
                            ! Check to see if we are near shore
                            if (q(3*layer-2,i,j) / rho(layer) < deep_depth) then
                                amrflags(i,j) = DOFLAG
                                cycle x_loop
                            ! Check if we are allowed to flag in deep water
                            ! anyway
                            else if (level < max_level_deep) then
                                amrflags(i,j) = DOFLAG
                                cycle x_loop
                            endif
                        endif

                        ! Check speed criteria, note that it might be useful to
                        ! also have a per layer criteria since this is not
                        ! gradient based
                        speed = sqrt(q(3*layer-1,i,j)**2 + q(3*layer,i,j)**2) / q(3*layer-2,i,j)
                        do m=1,min(size(speed_tolerance),mxnest)
                            if (speed > speed_tolerance(m) .and. level <= m) then
                                amrflags(i,j) = DOFLAG
                                cycle x_loop
                            endif
                        enddo
                        b = b + q(3*layer-2,i,j) / rho(layer)

                    endif
                enddo
            endif

        enddo x_loop
    enddo y_loop
end subroutine flag2refine2
