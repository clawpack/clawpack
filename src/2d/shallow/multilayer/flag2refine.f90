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
subroutine flag2refine(mx,my,mbc,meqn,maux,xlower,ylower,dx,dy,t,level,tolsp, &
                       q,aux,amrflags,DONTFLAG,DOFLAG)

    use amr_module, only: mxnest, t0
    use geoclaw_module, only:dry_tolerance,rho,eta_init,num_layers
    
    use topo_module, only: tlowtopo,thitopo,xlowtopo,xhitopo,ylowtopo,yhitopo
    use topo_module, only: minleveltopo,mtopofiles
    
    use dtopo_module, only: tfdtopo,xlowdtopo,xhidtopo,ylowdtopo,yhidtopo
    use dtopo_module, only: minleveldtopo,num_dtopo
    
    use qinit_module, only: x_low_qinit,x_hi_qinit,y_low_qinit,y_hi_qinit
    use qinit_module, only: min_level_qinit,qinit_type
    
    use regions_module, only: num_regions, regions
    use refinement_module
 
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
    integer :: i,j,m,k,layer_index
    real(kind=8) :: x_c,y_c,x_low,y_low,x_hi,y_hi
    real(kind=8) :: h(num_layers), total_depth, speed, eta, eta_below, ds

    ! Initialize flags
    amrflags = DONTFLAG

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

            ! Check to see if refinement is forced in any topography file region:
            do m=1,mtopofiles
                if (level < minleveltopo(m) .and. t >= tlowtopo(m) .and. t <= thitopo(m)) then
                    if (  x_hi > xlowtopo(m) .and. x_low < xhitopo(m) .and. &
                          y_hi > ylowtopo(m) .and. y_low < yhitopo(m) ) then
                        
                        amrflags(i,j) = DOFLAG
                        cycle x_loop
                    endif
                endif
            enddo

            ! Check to see if refinement is forced in any other region:
            do m=1,num_regions
                if (level < regions(m)%min_level .and. &
                    t >= regions(m)%t_low .and. t <= regions(m)%t_hi) then
                    if (x_hi > regions(m)%x_low .and. x_low < regions(m)%x_hi .and. &
                        y_hi > regions(m)%y_low .and. y_low < regions(m)%y_hi ) then
                    
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
                    y_hi > ylowdtopo(m) .and. y_low < yhidtopo(m)) then
                    
                    amrflags(i,j) = DOFLAG    
                    cycle x_loop
                endif
            enddo

            ! Check if we're in the region where initial perturbation is
            ! specified and need to force refinement:
            ! This assumes that t0 = 0.d0, should really be t0 but we do
            ! not have access to that parameter in this routine
            if (qinit_type > t0 .and. t == 0.d0) then 
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
            if (allowflag(x_c,y_c,t,level)) then
                ! These refinement criteria are checked per layer going backwards
                ! The bottom layer is checked first and eta_below is set to the 
                ! bathymetry
                eta_below = aux(1,i,j)
                
                forall(k=1:num_layers)
                    h(k) = q(3 * (k - 1) + 1,i,j) / rho(k)
                end forall
                total_depth = sum(h)

                do k=num_layers,1,-1
                    layer_index = 3 * (k - 1)
                    
                    ! Extract state
                    if (h(k) > dry_tolerance(k)) then
                        eta = h(k) + eta_below
                    
                        ! Check wave criteria
                        if (abs(eta - eta_init(k)) > wave_tolerance(k)) then
                            ! Check to see if we are near shore
                            if (total_depth < deep_depth) then
                                amrflags(i,j) = DOFLAG
                                cycle x_loop
                            ! Check if we are allowed to flag in deep water
                            ! anyway
                            else if (level < max_level_deep) then
                                amrflags(i,j) = DOFLAG
                                cycle x_loop
                            endif
                        endif
                        eta_below = eta
                    
                        ! Check speed criteria, note that it might be useful to 
                        ! also have a per layer criteria since this is not 
                        ! gradient based
                        speed = sqrt(q(layer_index+2,i,j)**2 &
                                   + q(layer_index+3,i,j)**2) &
                                   / q(layer_index+1,i,j)
                        do m=1,min(size(speed_tolerance),mxnest)
                            if (speed > speed_tolerance(m) .and. level <= m) then
                                amrflags(i,j) = DOFLAG
                                cycle x_loop
                            endif
                        enddo
                    endif
                enddo
            endif
            
        enddo x_loop
    enddo y_loop
end subroutine flag2refine
