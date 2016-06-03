! :::::::::::::::::::::::::::: FILPATCH :::::::::::::::::::::::::;
!
!  fill the portion of valbig from rows  nrowst
!                             and  cols  ncolst
!  the patch can also be described by the corners (xlp,ybp) by (xrp,ytp).
!  vals are needed at time t, and level level,
!
!  first fill with  values obtainable from the level level
!  grids. if any left unfilled, then enlarge remaining rectangle of
!  unfilled values by 1 (for later linear interp), and recusively
!  obtain the remaining values from  coarser levels.
!
! :::::::::::::::::::::::::::::::::::::::;:::::::::::::::::::::::;
recursive subroutine filrecur(level,num_eqn,valbig,aux,num_aux,t,mx,my, &
                              nrowst,ncolst,fill_indices)

    use amr_module, only: nghost, xlower, ylower, xupper, yupper, outunit
    use amr_module, only: xperdom, yperdom, spheredom, hxposs, hyposs
    use amr_module, only: intratx, intraty, iregsz, jregsz

    use geoclaw_module, only: sea_level, rho, dry_tolerance

    implicit none

    ! Input
    integer, intent(in) :: level, num_eqn, num_aux, mx, my, nrowst, ncolst
    integer, intent(in) :: fill_indices(4)
    real(kind=8), intent(in) :: t

    ! Output
    real(kind=8), intent(in out) :: valbig(num_eqn,mx,my)
    real(kind=8), intent(in out) :: aux(num_aux,mx,my)

    ! Local storage
    ! Flagging of set cells
    logical :: set
    integer :: i, i_coarse, j_coarse, i_fine, j_fine, n
    integer :: mx_coarse, my_coarse, mx_patch, my_patch
    integer :: unset_indices(4), coarse_indices(4)
    integer :: refinement_ratio_x, refinement_ratio_y
    real(kind=8) :: dx_fine, dy_fine, dx_coarse, dy_coarse
    real(kind=8) :: coarse_rect(4), fill_rect(4)
    real(kind=8) :: h, b, eta_fine, eta1, eta2, up_slope, down_slope
    real(kind=8) :: hv_fine, v_fine, v_new, divide_mass
    real(kind=8) :: h_fine_average, h_fine, h_count, h_coarse
    integer(kind=1) :: flaguse(fill_indices(2)-fill_indices(1)+1,fill_indices(4)-fill_indices(3)+1)
    
    ! Scratch arrays for interpolation
    logical :: fine_flag(num_eqn, fill_indices(2) - fill_indices(1) + 2, fill_indices(4) - fill_indices(3) + 2)

    logical :: reloop
    
    real(kind=8) :: fine_mass(fill_indices(2) - fill_indices(1) + 2, fill_indices(4) - fill_indices(3) + 2)
    real(kind=8) :: eta_coarse(fill_indices(2) - fill_indices(1) + 2, fill_indices(4) - fill_indices(3) + 2)
    real(kind=8) :: vel_max(fill_indices(2) - fill_indices(1) + 2, fill_indices(4) - fill_indices(3) + 2)
    real(kind=8) :: vel_min(fill_indices(2) - fill_indices(1) + 2, fill_indices(4) - fill_indices(3) + 2)
    real(kind=8) :: slope(2, fill_indices(2) - fill_indices(1) + 2, fill_indices(4) - fill_indices(3) + 2)
    integer :: fine_cell_count(fill_indices(2)-fill_indices(1)+2, fill_indices(4)-fill_indices(3)+2)
    integer :: nghost_patch

    ! Stack storage
    !  use stack-based scratch arrays instead of alloc, since dont really
    !  need to save beyond these routines, and to allow dynami_coarse memory resizing
    !  use 1d scratch arrays that are potentially the same size as
    !  current grid, since may not coarsen.
    !  need to make it 1d instead of 2 and do own indexing, since
    !  when pass it in to subroutines they treat it as having di_fineerent
    !  dimensions than the max size need to allocate here
    ! the +2 is to expand on coarse grid to enclose fine
    real(kind=8) :: valcrse((fill_indices(2)-fill_indices(1)+2) * (fill_indices(4)-fill_indices(3)+2) * num_eqn)  
    real(kind=8) :: auxcrse((fill_indices(2)-fill_indices(1)+2) * (fill_indices(4)-fill_indices(3)+2) * num_aux)  

    ! We begin by filling values for grids at level level.
    mx_patch = fill_indices(2) - fill_indices(1) + 1 ! nrowp
    my_patch = fill_indices(4) - fill_indices(3) + 1 ! ncolp

    dx_fine     = hxposs(level)
    dy_fine     = hyposs(level)

    ! Coordinates of edges of patch (xlp,xrp,ybp,ytp)
    fill_rect = [xlower + fill_indices(1) * dx_fine, &
                 xlower + (fill_indices(2) + 1) * dx_fine, &
                 ylower + fill_indices(3) * dy_fine, &
                 ylower + (fill_indices(4) + 1) * dy_fine]

    ! Fill in the patch as much as possible using values at this level
    call intfil(valbig,mx,my,t,flaguse,nrowst,ncolst,fill_indices(1), &
                fill_indices(2),fill_indices(3),fill_indices(4),level,num_eqn,num_aux)

    ! Trimbd returns set = true if all of the entries are filled (=1.).
    ! set = false, otherwise. If set = true, then no other levels are
    ! are required to interpolate, and we return.
    !
    ! Note that the used array is filled entirely in intfil, i.e. the
    ! marking done there also takes  into account the points filled by
    ! the boundary conditions. bc2amr will be called later, after all 4
    ! boundary pieces filled.
    call trimbd(flaguse,mx_patch,my_patch,set,unset_indices)
    ! il,ir,jb,jt = unset_indices(4)

    ! If set is .true. then all cells have been set and we can skip to setting
    ! the remaining boundary cells.  If it is .false. we need to interpolate
    ! some values from coarser levels, possibly calling this routine
    ! recursively.
    if (.not.set) then

        ! Error check 
        if (level == 1) then
            write(outunit,*)" error in filrecur - level 1 not set"
            write(outunit,'("start at row: ",i4," col ",i4)') nrowst,ncolst
            print *," error in filrecur - level 1 not set"
            print *," should not need more recursion "
            print *," to set patch boundaries"
            print '("start at row: ",i4," col ",i4)', nrowst,ncolst
            stop
        endif

        ! We begin by initializing the level level arrays, so that we can use
        ! purely recursive formulation for interpolating.
        dx_coarse  = hxposs(level - 1)
        dy_coarse  = hyposs(level - 1)

        ! Adjust unset_indices to account for the patch
        ! isl, isr, jsb, jst
        unset_indices(1) = unset_indices(1) + fill_indices(1) - 1
        unset_indices(2) = unset_indices(2) + fill_indices(1) - 1
        unset_indices(3) = unset_indices(3) + fill_indices(3) - 1
        unset_indices(4) = unset_indices(4) + fill_indices(3) - 1

        ! Coarsened geometry
        refinement_ratio_x = intratx(level - 1)
        refinement_ratio_y = intraty(level - 1)

        ! New patch rectangle (after we have partially filled it in) but in the
        ! coarse patches [iplo,iphi,jplo,jphi]
        coarse_indices = [(unset_indices(1) - refinement_ratio_x + nghost * refinement_ratio_x) &
                                                / refinement_ratio_x - nghost, &
                          (unset_indices(2) + refinement_ratio_x) / refinement_ratio_x, &
                          (unset_indices(3) - refinement_ratio_y + nghost * refinement_ratio_y) &
                                                / refinement_ratio_y - nghost, &
                          (unset_indices(4) + refinement_ratio_y) / refinement_ratio_y]
        coarse_rect = [xlower + coarse_indices(1) * dx_coarse, &
                       xlower + (coarse_indices(2) + 1) * dx_coarse, &
                       ylower + coarse_indices(3) * dy_coarse, &
                       ylower + (coarse_indices(4) + 1) * dy_coarse]

        ! Coarse grid number of spatial points (nrowc,ncolc)
        mx_coarse   =  coarse_indices(2) - coarse_indices(1) + 1
        my_coarse   =  coarse_indices(4) - coarse_indices(3) + 1

        ! Check to make sure we created big enough scratch arrays
        if (mx_coarse > fill_indices(2) - fill_indices(1) + 3 .or. &
            my_coarse > fill_indices(4) - fill_indices(3) + 3) then

            print *," did not make big enough work space in filrecur "
            print *," need coarse space with mx_coarse,my_coarse ",mx_coarse,my_coarse
            print *," made space for ilo,fill_indices(2),fill_indices(3),fill_indices(4) ",fill_indices
            stop
        endif

        ! Set the aux array values for the coarse grid, this could be done 
        ! instead in intfil using possibly already available bathy data from the
        ! grids
        if (num_aux > 0) then
            nghost_patch = 0
            call setaux(nghost_patch, mx_coarse, my_coarse, &
                        coarse_rect(1), coarse_rect(3), &
                        dx_coarse,dy_coarse,num_aux,auxcrse)
        endif

        ! Fill in the edges of the coarse grid
        if ((xperdom .or. (yperdom .or. spheredom)) .and. sticksout(coarse_indices)) then
            call prefilrecur(level - 1,num_eqn,valcrse,auxcrse,num_aux,t,mx_coarse,my_coarse,1,1,coarse_indices)
        else
            call filrecur(level - 1,num_eqn,valcrse,auxcrse,num_aux,t,mx_coarse,my_coarse,1,1,coarse_indices)
        endif

        ! loop through coarse cells determining intepolation slopes
        ! these will be saved for fine grid loop
        ! prevents finding the same slope possibly lratiox*lratioy times
        ! all fine gid depths will be found before any momentum
        reloop = .false.
        fine_cell_count = 0
        fine_flag = .false.
        fine_mass = 0.d0
        slope = 0.d0
        
        ! Calculate surface elevation eta using dry limiting
        do i_coarse = 1, mx_coarse
            do j_coarse = 1, my_coarse
                h = valcrse(ivalc(1,i_coarse,j_coarse)) / rho(1)
                b = auxcrse(iauxc(i_coarse,j_coarse))

                if (h < dry_tolerance(1)) then
                    eta_coarse(i_coarse,j_coarse) = eta_init(1)
                else
                    eta_coarse(i_coarse,j_coarse) = h + b
                endif
            enddo
        enddo

        ! Calculate limited gradients of coarse grid eta
        do i_coarse = 2, mx_coarse - 1
            do j_coarse = 2, my_coarse - 1 
                
                ! X-Direction
                down_slope = eta_coarse(i_coarse,j_coarse) - eta_coarse(i_coarse-1,j_coarse)
                up_slope = eta_coarse(i_coarse+1,j_coarse) - eta_coarse(i_coarse,j_coarse)
                if (up_slope * down_slope > 0.d0) then
                    slope(1,i_coarse,j_coarse) = min(abs(up_slope), abs(down_slope)) &
                        * sign(1.d0,eta_coarse(i_coarse+1,j_coarse) - eta_coarse(i_coarse-1,j_coarse))
                endif

                ! Y-Direction
                down_slope = eta_coarse(i_coarse,j_coarse) - eta_coarse(i_coarse,j_coarse-1)
                up_slope = eta_coarse(i_coarse,j_coarse+1) - eta_coarse(i_coarse,j_coarse)
                if (up_slope * down_slope > 0.d0) then
                    slope(2,i_coarse,j_coarse) = min(abs(up_slope), abs(down_slope)) &
                        * sign(1.d0,eta_coarse(i_coarse+1,j_coarse) - eta_coarse(i_coarse-1,j_coarse))
                endif
            enddo
        enddo

        ! Loop through patch to be filled, includes multiple coarse cells
        do i_fine = 1, mx_patch
            i_coarse = 2 + (i_fine - (unset_indices(1) - fill_indices(1)) - 1) / refinement_ratio_x
            eta1 = (-0.5d0 + real(mod(i_fine - 1, refinement_ratio_x),kind=8)) &
                                / real(refinement_ratio_x,kind=8)
            do j_fine = 1, my_patch
                j_coarse = 2 + (j_fine - (unset_indices(3) - fill_indices(3)) - 1) / refinement_ratio_y
                eta2 = (-0.5d0 + real(mod(j_fine - 1, refinement_ratio_y),kind=8)) &
                                    / real(refinement_ratio_y,kind=8)

                if (flaguse(i_fine,j_fine) == 0) then
                    ! Interpolate from coarse cells to fine grid for surface
                    fine_cell_count(i_coarse,j_coarse) = &
                                        fine_cell_count(i_coarse,j_coarse) + 1
                    eta_fine = eta_coarse(i_coarse,j_coarse) + eta1 * slope(1,i_coarse,j_coarse) &
                                                             + eta2 * slope(2,i_coarse,j_coarse)
                    h_fine = max(eta_fine - aux(1,i_fine + nrowst - 1, j_fine + ncolst - 1), 0.d0)
                    valbig(1,i_fine + nrowst - 1, j_fine + ncolst - 1) = h_fine * rho(1)
                    fine_mass(i_coarse,j_coarse) = fine_mass(i_coarse,j_coarse) + h_fine
                    
                    ! Flag the corresponding coarse cell as needing relimiting
                    ! if one of the fine cells ends up being dry
                    if (h_fine < dry_tolerance(1)) then
                        fine_flag(1,i_coarse,j_coarse) = .true.
                        reloop = .true.
                    endif
                endif
            enddo
        enddo

        ! Momentum Interpolation
        do n = 2, num_eqn
            do i_coarse = 2, mx_coarse - 1
                do j_coarse = 2, my_coarse - 1

                    ! Determine slopes for interpolation
                    down_slope = (valcrse(ivalc(n,i_coarse,j_coarse)) - valcrse(ivalc(n,i_coarse-1,j_coarse))) / rho(1)
                    up_slope = (valcrse(ivalc(n,i_coarse+1,j_coarse)) - valcrse(ivalc(n,i_coarse,j_coarse))) / rho(1)
                    if (up_slope * down_slope > 0.d0) then
                        slope(1,i_coarse,j_coarse) = &
                                        min(abs(up_slope), abs(down_slope)) * &
                             sign(1.d0, valcrse(ivalc(n,i_coarse+1,j_coarse)) &
                           - valcrse(ivalc(n,i_coarse-1,j_coarse)))
                    endif

                    down_slope = (valcrse(ivalc(n,i_coarse,j_coarse)) - valcrse(ivalc(n,i_coarse,j_coarse-1))) / rho(1)
                    up_slope = (valcrse(ivalc(n,i_coarse,j_coarse+1)) - valcrse(ivalc(n,i_coarse,j_coarse))) / rho(1)
                    if (up_slope * down_slope > 0.d0) then
                        slope(2,i_coarse,j_coarse) = &
                                          min(abs(up_slope), abs(down_slope)) &
                            * sign(1.d0, valcrse(ivalc(n,i_coarse,j_coarse+1)) &
                                - valcrse(ivalc(n,i_coarse,j_coarse-1)))
                    endif

                    ! Set initial values for max/min of current field
                    if (valcrse(ivalc(1,i_coarse,j_coarse)) > dry_tolerance(1)) then
                        vel_max(i_coarse,j_coarse) = valcrse(ivalc(n,i_coarse,j_coarse)) / valcrse(ivalc(1,i_coarse,j_coarse))
                        vel_min(i_coarse,j_coarse) = valcrse(ivalc(n,i_coarse,j_coarse)) / valcrse(ivalc(1,i_coarse,j_coarse))
                    else
                        vel_min(i_coarse,j_coarse) = 0.d0
                        vel_max(i_coarse,j_coarse) = 0.d0
                    endif

                    ! Look for bounds on velocity around each cell
                    ! Necessary since we are interpolating momentum linearly
                    ! but not interpolating depth linearly
                    do i =-1,1,2
                        if (valcrse(ivalc(1,i_coarse + i,j_coarse)) / rho(1) > dry_tolerance(1)) then
                            vel_max(i_coarse,j_coarse) = &
                                max(vel_max(i_coarse,j_coarse), &
                                    valcrse(ivalc(n,i_coarse + i,j_coarse)) &
                                  / valcrse(ivalc(1,i_coarse + i,j_coarse)))
                            vel_min(i_coarse,j_coarse) = &
                                min(vel_min(i_coarse,j_coarse), &
                                    valcrse(ivalc(n,i_coarse + i,j_coarse)) &
                                  / valcrse(ivalc(1,i_coarse + i,j_coarse)))
                        endif
                        if (valcrse(ivalc(1,i_coarse,j_coarse + i)) / rho(1) > dry_tolerance(1)) then
                            vel_max(i_coarse,j_coarse) = &
                                max(vel_max(i_coarse,j_coarse), &
                                      valcrse(ivalc(n,i_coarse,j_coarse + i)) &
                                    / valcrse(ivalc(1,i_coarse,j_coarse + i)))
                            vel_min(i_coarse,j_coarse) = &
                                min(vel_min(i_coarse,j_coarse), &
                                      valcrse(ivalc(n,i_coarse,j_coarse + i)) &
                                    / valcrse(ivalc(1,i_coarse,j_coarse + i)))

                        endif
                    enddo
                enddo
            enddo

            ! Determine momentum in fine cells
            do i_fine = 1, mx_patch
                i_coarse = 2 + (i_fine - (unset_indices(1) - fill_indices(1)) - 1) / refinement_ratio_x
                eta1 = (-0.5d0 + real(mod(i_fine - 1, refinement_ratio_x),kind=8)) / real(refinement_ratio_x,kind=8)
                do j_fine = 1, my_patch
                    j_coarse = 2 + (j_fine - (unset_indices(3) - fill_indices(3)) - 1) / refinement_ratio_y
                    eta2 = (-0.5d0 + real(mod(j_fine - 1, refinement_ratio_y),kind=8)) / real(refinement_ratio_y,kind=8)

                    if (flaguse(i_fine,j_fine) == 0) then
                        ! Cell not already set

                        if (.not.(fine_flag(1,i_coarse,j_coarse))) then
                            ! This cell has no coarse cells that are dry
                            hv_fine = valcrse(ivalc(n,i_coarse,j_coarse)) / rho(1) &
                                            + eta1 * slope(1,i_coarse,j_coarse) &
                                            + eta2 * slope(2,i_coarse,j_coarse)
                            v_fine = hv_fine * rho(1) / valbig(1,i_fine + nrowst - 1, j_fine + ncolst - 1)
                            if (v_fine < vel_min(i_coarse,j_coarse) .or. &
                                v_fine > vel_max(i_coarse,j_coarse)) then

                                fine_flag(n,i_coarse,j_coarse) = .true.
                                reloop = .true.
                            else
                                valbig(n,i_fine+nrowst-1,j_fine+ncolst-1) = hv_fine * rho(1)
                            endif
                        endif
                    endif
                enddo
            enddo

            ! Reset momentum to conserve momentum in the cases where we may have
            ! gained momentum or if velocity bounds were violated
            if (reloop) then
                do i_fine = 1,mx_patch
                    i_coarse = 2 + (i_fine - (unset_indices(1) - fill_indices(1)) - 1) / refinement_ratio_x
                    eta1 = (-0.5d0 + real(mod(i_fine - 1, refinement_ratio_x),kind=8)) / real(refinement_ratio_x,kind=8)
                    do j_fine  = 1,my_patch
                        j_coarse = 2 + (j_fine  - (unset_indices(3) - fill_indices(3)) - 1) / refinement_ratio_y
                        eta2 = (-0.5d0 + real(mod(j_fine - 1, refinement_ratio_y),kind=8)) / real(refinement_ratio_y, kind=8)
                        
                        if (flaguse(i_fine,j_fine) == 0) then
                            if (fine_flag(1,i_coarse,j_coarse) .or. fine_flag(n,i_coarse,j_coarse)) then
                                if (fine_mass(i_coarse,j_coarse) > dry_tolerance(1)) then

                                    h_coarse = valcrse(ivalc(1,i_coarse,j_coarse)) / rho(1)
                                    h_count = real(fine_cell_count(i_coarse,j_coarse),kind=8)
                                    h_fine_average = fine_mass(i_coarse,j_coarse) / h_count
                                    divide_mass = max(h_coarse, h_fine_average)
                                    h_fine = valbig(1, i_fine + nrowst - 1, j_fine + ncolst - 1) / rho(1)
                                    v_new = valcrse(ivalc(n,i_coarse,j_coarse)) / (rho(1) * divide_mass)
                                    
                                    valbig(n,i_fine+nrowst-1,j_fine+ncolst-1) = &
                                        v_new * valbig(1,i_fine+nrowst-1,j_fine+ncolst-1)
                                else
                                    valbig(n,i_fine+nrowst-1,j_fine+ncolst-1) = 0.d0
                                endif
                            endif
                        endif
                    enddo
                enddo
            endif
        enddo
    endif

    ! set bcs, whether or not recursive calls needed. set any part of patch that
    ! stuck out
    call bc2amr(valbig,aux,mx,my,num_eqn,num_aux,dx_fine,dy_fine,level,t,fill_rect(1), &
                fill_rect(2),fill_rect(3),fill_rect(4))

contains

    integer pure function ivalc(n,i,j)
        implicit none
        integer, intent(in) :: n,i,j
        ivalc = n + num_eqn*(i-1) + num_eqn*mx_coarse*(j-1)
    end function ivalc

    ! Index into first component of aux = topo:
    integer pure function iauxc(i,j)
        implicit none
        integer, intent(in) :: i,j
        iauxc = 1 + num_aux*(i-1) + num_aux*mx_coarse*(j-1)
    end function iauxc

    ! logical for checking if this patch sti_coarseks outside of the domain
    logical pure function sticksout(rect)
        implicit none
        integer, intent(in) :: rect(4)
        sticksout = (rect(1) < 0 .or. rect(3) < 0 .or. &
                     rect(2) >= iregsz(level - 1) .or. rect(4) >= jregsz(level - 1))
    end function sticksout

end subroutine filrecur
