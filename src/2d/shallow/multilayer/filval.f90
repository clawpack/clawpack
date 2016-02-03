!
! :::::::::::::::::::::::::::::: FILVAL ::::::::::::::::::::::::::
!
! create and fill coarser (level-1) patch with one extra coarse cell all
! around, plus the ghost cells . will interpolate from this patch to grid mptr
! without needing special boundary code.
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
!
! ------------------------------------------------------------------
subroutine filval(val, mx, my, dx, dy, level, time, valc, auxc, mic, &
                  mjc, xleft, xright, ybot, ytop, nvar, mptr, ilo, ihi, &
                  jlo, jhi, aux, naux, locflip, sp_over_h)

    use amr_module, only: xlower, ylower, intratx, intraty, nghost, xperdom
    use amr_module, only: yperdom, spheredom, xupper, yupper

    use geoclaw_module, only: rho, dry_tolerance, eta_init, varRefTime

    implicit none

    ! Input
    integer, intent(in) :: mx, my, level, mic, mjc, nvar, mptr, ilo, ihi
    integer, intent(in) :: jlo, jhi, naux, locflip
    real(kind=8), intent(in) :: dx, dy, time, xleft, xright, ybot, ytop
    real(kind=8), intent(in) :: valc(nvar,mic,mjc), auxc(naux,mic,mjc)

    ! Output
    real(kind=8), intent(in out) :: sp_over_h
    real(kind=8), intent(in out) :: val(nvar,mx,my), aux(naux,mx,my)

    ! Local storage
    integer :: refinement_ratio_x, refinement_ratio_y, iclo, jclo, ichi, jchi, ng, i, ico, ifine
    integer :: ii, ivar, j, jco, jfine, jj
    real(kind=8) :: coarseval(3), dx_coarse, dy_coarse, xl, xr, yb, yt, area
    real(kind=8) :: dividemass, finemass, hvf, s1m, s1p, slopex, slopey, vel
    real(kind=8) :: velmax, velmin, vf, vnew, xoff, yoff
    logical :: fineflag(3)

    ! External function definitions
    real(kind=8) :: get_max_speed

    refinement_ratio_x = intratx(level - 1)
    refinement_ratio_y = intraty(level - 1)
    dx_coarse  = dx * refinement_ratio_x
    dy_coarse  = dy * refinement_ratio_y
    xl      = xleft  - dx_coarse
    xr      = xright + dx_coarse
    yb      = ybot   - dy_coarse
    yt      = ytop   + dy_coarse

    ! set integer indices for coarser patch enlarged by 1 cell
    ! (can stick out of domain). proper nesting will insure this one
    ! call is sufficient.
    iclo   = ilo / refinement_ratio_x - 1
    jclo   = jlo / refinement_ratio_y - 1
    ichi   = (ihi + 1) / refinement_ratio_x - 1 + 1
    jchi   = (jhi + 1) / refinement_ratio_y - 1 + 1
    ng     = 0

    if (naux == 0) then
        if (xperdom .or. yperdom .or. spheredom) then
            call preintcopy(valc,mic,mjc,nvar,iclo,ichi,jclo,jchi,level - 1,locflip)
        else
            call intcopy(valc,mic,mjc,nvar,iclo,ichi,jclo,jchi,level - 1,1,1)
        endif
    else  
        ! intersect grids and copy all (soln and aux)
        if (xperdom .or. yperdom .or. spheredom) then
            call preicall(valc,auxc,mic,mjc,nvar,naux,iclo,ichi,jclo,jchi, &
                          level - 1,locflip)
        else
            call icall(valc,auxc,mic,mjc,nvar,naux,iclo,ichi,jclo,jchi,level - 1,1,1)
        endif
    endif
    call bc2amr(valc,auxc,mic,mjc,nvar,naux,dx_coarse,dy_coarse,level - 1,time,  &
                xl,xr,yb,yt)

    !-----------------------------
    ! For shallow water over topograpdy, in coarse cells convert from h to eta,
    ! before interpolating:
    !-----------------------------

    ! Prepare slopes - use min-mod limiters
    do j=2, mjc-1
        do i=2, mic-1
            fineflag(1) = .false.
            ! interpolate eta to find depth
            do ii=-1,1
                coarseval(2+ii) = valc(1,i+ii,j) / rho(1) + auxc(1,i+ii,j)
                if (valc(1,i+ii,j) / rho(1) <= dry_tolerance(1)) then
                    coarseval(2+ii)=eta_init(1)
                end if
            end do
            s1p = coarseval(3) - coarseval(2)
            s1m = coarseval(2) - coarseval(1)
            slopex = min(abs(s1p), abs(s1m)) &
                                * sign(1.d0,coarseval(3) - coarseval(1))
            if (s1m*s1p <= 0.d0) slopex=0.d0

            do jj=-1,1
                coarseval(2+jj) = valc(1,i,j+jj) / rho(1) + auxc(1,i,j+jj)
                if (valc(1,i,j+jj) / rho(1) <= dry_tolerance(1)) then
                    coarseval(2+jj)=eta_init(1)
                end if
            end do
            s1p = coarseval(3) - coarseval(2)
            s1m = coarseval(2) - coarseval(1)
            slopey = min(abs(s1p), abs(s1m)) &
                                * sign(1.d0,coarseval(3)-coarseval(1))
            if (s1m*s1p <= 0.d0) slopey=0.d0

            ! Interpolate from coarse cells to fine grid to find depth
            finemass = 0.d0
            do ico = 1,refinement_ratio_x
                do jco = 1,refinement_ratio_y
                    yoff = (real(jco,kind=8) - 0.5d0) / refinement_ratio_y - 0.5d0
                    xoff = (real(ico,kind=8) - 0.5d0) / refinement_ratio_x - 0.5d0
                    jfine = (j-2) * refinement_ratio_y + nghost + jco
                    ifine = (i-2) * refinement_ratio_x + nghost + ico
                    val(1,ifine,jfine) = (coarseval(2) + xoff * slopex &
                                                       + yoff * slopey) * rho(1)
                    val(1,ifine,jfine) = max(0.d0, val(1,ifine,jfine) / rho(1) &
                                            - aux(1,ifine,jfine)) * rho(1)
                    finemass = finemass + val(1,ifine,jfine) / rho(1)
                    if (val(1,ifine,jfine) / rho(1) <= dry_tolerance(1)) then
                        fineflag(1) = .true.
                        val(2,ifine,jfine) = 0.d0
                        val(3,ifine,jfine) = 0.d0
                    end if
                end do
            end do

            !------ Determine Momentum ----------------------------------
            ! finemass is the total mass in all new fine grid cells
            ! all fine mass has been determined for this coarse grid cell
            ! if all fine cells are dry, momentum has already been set
            if (finemass >= dry_tolerance(1)) then
                do ivar = 2,nvar
                    fineflag(ivar)=.false.
                    s1p = (valc(ivar,i+1,j) - valc(ivar,i,j)) / rho(1)
                    s1m = (valc(ivar,i,j) - valc(ivar,i-1,j)) / rho(1)
                    slopex = min(abs(s1p), abs(s1m)) &
                     * sign(1.d0,(valc(ivar,i+1,j) - valc(ivar,i-1,j)) / rho(1))
                    if (s1m*s1p.le.0.d0) slopex=0.d0
                
                    s1p = (valc(ivar,i,j+1) - valc(ivar,i,j)) / rho(1)
                    s1m = (valc(ivar,i,j) - valc(ivar,i,j-1)) / rho(1)
                    slopey = min(abs(s1p), abs(s1m)) &
                     * sign(1.d0,(valc(ivar,i,j+1) - valc(ivar,i,j-1)) / rho(1))
                    if (s1m*s1p.le.0.d0) slopey=0.d0

                    if (valc(1,i,j) / rho(1) > dry_tolerance(1)) then
                        velmax = valc(ivar,i,j)/valc(1,i,j)
                        velmin = valc(ivar,i,j)/valc(1,i,j)
                    else
                        velmax = 0.d0
                        velmin = 0.d0
                    endif
               
                    do ii = -1,1,2
                        if (valc(1,i+ii,j) / rho(1) > dry_tolerance(1)) then
                            vel = valc(ivar,i+ii,j) / valc(1,i+ii,j)
                            velmax = max(vel,velmax)
                            velmin = min(vel,velmin)
                        endif
                        if (valc(1,i,j+ii) / rho(1) > dry_tolerance(1)) then
                            vel = valc(ivar,i,j+ii) / valc(1,i,j+ii)
                            velmax = max(vel,velmax)
                            velmin = min(vel,velmin)
                        endif
                    enddo

                    ! try to set momentum
                    do ico = 1,refinement_ratio_x
                        if (fineflag(1).or.fineflag(ivar)) exit

                        do jco = 1,refinement_ratio_y
                            jfine = (j-2) * refinement_ratio_y + nghost + jco
                            ifine = (i-2) * refinement_ratio_x + nghost + ico
                            yoff = (real(jco,kind=8) - 0.5d0) / refinement_ratio_y - 0.5d0
                            xoff = (real(ico,kind=8) - 0.5d0) / refinement_ratio_x - 0.5d0
                            hvf = valc(ivar,i,j) / rho(1) + xoff * slopex &
                                                          + yoff*slopey
                            vf = hvf / (val(1,ifine,jfine) / rho(1))
                            if (vf > velmax .or. vf < velmin) then
                                fineflag(ivar) = .true.
                                exit
                            else
                                val(ivar,ifine,jfine) = hvf * rho(1)
                            endif
                        enddo
                    enddo

                    ! momentum is set to preserve old momentum or not violate
                    ! generating new extrema in velocities
                    if (fineflag(1) .or. fineflag(ivar)) then
                        ! more mass now, conserve momentum
                        area = real(refinement_ratio_x * refinement_ratio_y,kind=8)
                        dividemass = max(finemass,valc(1,i,j) / rho(1))
                        Vnew = area * valc(ivar,i,j) / (dividemass * rho(1))

                        do ico = 1,refinement_ratio_x
                            do jco = 1,refinement_ratio_y
                                jfine = (j-2) * refinement_ratio_y + nghost + jco
                                ifine = (i-2) * refinement_ratio_x + nghost + ico
                                val(ivar,ifine,jfine) = Vnew * val(1,ifine,jfine)
                            enddo
                        enddo
                    endif

                enddo
            endif

        enddo !end of coarse loop
    enddo !end of coarse loop

    ! overwrite interpolated values with fine grid values, if available.
    call intcopy(val,mx,my,nvar,ilo-nghost,ihi+nghost,jlo-nghost, &
                 jhi+nghost,level,1,1)

    ! scan for max wave speed on newly created grid. this will be used to set appropriate
    ! time step and appropriate refinement in time. For this app not nec to refine by same
    ! amount as refinement in space since refinement at shores where h is shallow has lower
    ! speeds.

    if (varRefTime) then   ! keep consistent with setgrd_geo and qinit_geo
        sp_over_h = get_max_speed(val,mx,my,nvar,aux,naux,nghost,dx,dy)
    endif

end subroutine filval
