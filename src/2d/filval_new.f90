!
! :::::::::::::::::::::::::::::: FILVAL ::::::::::::::::::::::::::
!
! create and fill coarser (lev-1) patch with one extra coarse cell all
! around, plus the ghost cells . will interpolate from this patch to grid mptr
! without needing special boundary code.
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
!
! ------------------------------------------------------------------
subroutine filval(val, mitot, mjtot, hx, hy, lev, time, valc, auxc, mic, &
                  mjc, xleft, xright, ybot, ytop, nvar, mptr, ilo, ihi, &
                  jlo, jhi, aux, naux, locflip, sp_over_h)

    use amr_module, only: xlower, ylower, intratx, intraty, nghost, xperdom
    use amr_module, only: yperdom, spheredom, xupper, yupper

    use geoclaw_module, only: rho, dry_tolerance, eta_init, varRefTime

    implicit none

    ! Input
    integer, intent(in) :: mitot, mjtot, lev, mic, mjc, nvar, mptr, ilo, ihi
    integer, intent(in) :: jlo, jhi, naux, locflip
    real(kind=8), intent(in) :: hx, hy, time, xleft, xright, ybot, ytop
    real(kind=8), intent(in) :: valc(nvar,mic,mjc), auxc(naux,mic,mjc)

    ! Output
    real(kind=8), intent(in out) :: sp_over_h
    real(kind=8), intent(in out) :: val(nvar,mitot,mjtot), aux(naux,mitot,mjtot)

    ! Local storage
    integer :: levc, lratiox, lratioy, iclo, jclo, ichi, jchi, ng, i, ico, ifine
    integer :: ii, ivar, j, jco, jfine, jj
    real(kind=8) :: coarseval(3), hxcrse, hycrse, xl, xr, yb, yt, area
    real(kind=8) :: dividemass, finemass, hvf, s1m, s1p, slopex, slopey, vel
    real(kind=8) :: velmax, velmin, vf, vnew, xoff, yoff
    logical :: fineflag(3)

    ! External function definitions
    real(kind=8) :: get_max_speed


    ! indext into eta array for surface values:
    ! iaddeta(i,j) = loceta + i-1 + mic*(j-1)

    levc    = lev - 1
    lratiox = intratx(levc)
    lratioy = intraty(levc)
    hxcrse  = hx * lratiox
    hycrse  = hy * lratioy
    xl      = xleft  - hxcrse
    xr      = xright + hxcrse
    yb      = ybot   - hycrse
    yt      = ytop   + hycrse

    ! set integer indices for coarser patch enlarged by 1 cell
    ! (can stick out of domain). proper nesting will insure this one
    ! call is sufficient.
    iclo   = ilo / lratiox - 1
    jclo   = jlo / lratioy - 1
    ichi   = (ihi + 1) / lratiox - 1 + 1
    jchi   = (jhi + 1) / lratioy - 1 + 1
    ng     = 0

    if (naux == 0) then
        if (xperdom .or. yperdom .or. spheredom) then
            call preintcopy(valc,mic,mjc,nvar,iclo,ichi,jclo,jchi,levc,locflip)
        else
            call intcopy(valc,mic,mjc,nvar,iclo,ichi,jclo,jchi,levc,1,1)
        endif
    else  
        ! intersect grids and copy all (soln and aux)
        if (xperdom .or. yperdom .or. spheredom) then
            call preicall(valc,auxc,mic,mjc,nvar,naux,iclo,ichi,jclo,jchi, &
                          levc,locflip)
        else
            call icall(valc,auxc,mic,mjc,nvar,naux,iclo,ichi,jclo,jchi,levc,1,1)
        endif
    endif
    call bc2amr(valc,auxc,mic,mjc,nvar,naux,hxcrse,hycrse,levc,time,xl,xr,yb, &
                yt,xlower,ylower,xupper,yupper,xperdom,yperdom,spheredom)

    !-----------------------------
    ! For shallow water over topography, in coarse cells convert from h to eta,
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
            do ico = 1,lratiox
                do jco = 1,lratioy
                    yoff = (real(jco,kind=8) - 0.5d0) / lratioy - 0.5d0
                    xoff = (real(ico,kind=8) - 0.5d0) / lratiox - 0.5d0
                    jfine = (j-2) * lratioy + nghost + jco
                    ifine = (i-2) * lratiox + nghost + ico
                    val(1,ifine,jfine) = (coarseval(2) + xoff * slopex &
                                                       + yoff * slopey) * rho(1)
                    val(1,ifine,jfine) = max(0.d0, val(1,ifine,jfine) / rho(1) &
                                            - aux(1,ifine,jfine)) * rho(1)
                    finemass = finemass + val(1,ifine,jfine)
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
            if (finemass / rho(1) >= dry_tolerance(1)) then
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
                    do ico = 1,lratiox
                        if (fineflag(1).or.fineflag(ivar)) exit

                        do jco = 1,lratioy
                            jfine = (j-2) * lratioy + nghost + jco
                            ifine = (i-2) * lratiox + nghost + ico
                            yoff = (real(jco,kind=8) - 0.5d0) / lratioy - 0.5d0
                            xoff = (real(ico,kind=8) - 0.5d0) / lratiox - 0.5d0
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
                        area = real(lratiox * lratioy,kind=8)
                        dividemass = max(finemass,valc(1,i,j))
                        Vnew = area * valc(ivar,i,j) / dividemass

                        do ico = 1,lratiox
                            do jco = 1,lratioy
                                jfine = (j-2) * lratioy + nghost + jco
                                ifine = (i-2) * lratiox + nghost + ico
                                val(ivar,ifine,jfine) = Vnew * val(1,ifine,jfine)
                            enddo
                        enddo
                    endif

                enddo
            endif

        enddo !end of coarse loop
    enddo !end of coarse loop

    ! overwrite interpolated values with fine grid values, if available.
    call intcopy(val,mitot,mjtot,nvar,ilo-nghost,ihi+nghost,jlo-nghost, &
                 jhi+nghost,lev,1,1)

    ! scan for max wave speed on newly created grid. this will be used to set appropriate
    ! time step and appropriate refinement in time. For this app not nec to refine by same
    ! amount as refinement in space since refinement at shores where h is shallow has lower
    ! speeds.

    if (varRefTime) then   ! keep consistent with setgrd_geo and qinit_geo
        sp_over_h = get_max_speed(val,mitot,mjtot,nvar,aux,naux,nghost,hx,hy)
    endif

end subroutine filval
