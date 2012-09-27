! :::::::::::::::::::::::::::: FIlPATCH :::::::::::::::::::::::::;
!
!  fill the portion of valbig from rows  nrowst
!                             and  cols  ncolst
!  the patch can also be described by the corners (xlp,ybp) by (xrp,ytp).
!  vals are needed at time time, and level level,
!
!  first fill with  values obtainable from the level level
!  grids. if any left unfilled, then enlarge remaining rectangle of
!  unfilled values by 1 (for later linear interp), and recusively
!  obtain the remaining values from  coarser levels.
!
! :::::::::::::::::::::::::::::::::::::::;:::::::::::::::::::::::;
!
! Note that there are internal functions that this subroutine contains
!
recursive subroutine filrecur(level,nvar,valbig,aux,naux,time,mitot,mjtot, &
                                nrowst,ncolst,ilo,ihi,jlo,jhi)

    use amr_module, only: xlower, ylower, hxposs, hyposs, intratx, intraty
    use amr_module, only: nghost, outunit, xperdom, yperdom, spheredom
    use amr_module, only: iregsz, jregsz, xupper, yupper
    
    use geoclaw_module, only: rho, dry_tolerance, eta_init

    implicit none

    ! Input
    integer, intent(in) :: level, nvar, naux, mitot, mjtot, nrowst, ncolst
    integer, intent(in) :: ilo, ihi, jlo, jhi
    real(kind=8), intent(in) :: time

    real(kind=8), intent(in out) :: valbig(nvar,mitot,mjtot)
    real(kind=8), intent(in out) :: aux(naux,mitot,mjtot)

    ! Scratch arrays
    real(kind=8) :: valcrse((ihi-ilo+2)*(jhi-jlo+2)*nvar)  ! NB this is a 1D array
    real(kind=8) :: auxcrse((ihi-ilo+2)*(jhi-jlo+2)*naux)  ! the +2 is to expand on coarse grid to enclose fine

    real(kind=8) :: flaguse(ihi-ilo+1,jhi-jlo+1)

    logical :: fineflag((ihi-ilo+2)*(jhi-jlo+2)*nvar)
    real(kind=8) :: finemass((ihi-ilo+2)*(jhi-jlo+2))
    real(kind=8) :: etacrse((ihi-ilo+2)*(jhi-jlo+2))
    real(kind=8) :: velmax((ihi-ilo+2)*(jhi-jlo+2))
    real(kind=8) :: velmin((ihi-ilo+2)*(jhi-jlo+2))
    real(kind=8) :: slopex((ihi-ilo+2)*(jhi-jlo+2))
    real(kind=8) :: slopey((ihi-ilo+2)*(jhi-jlo+2))
    integer :: icount((ihi-ilo+2)*(jhi-jlo+2))

    ! Local scalars
    integer :: i, j, ic, jc, iff, ii, cc, il, ir, jb, jt, iphi, iplo, isl, isr
    integer :: ivar, jf, jphi, jplo, jsb, jst, levc, ncolc, ncolp, nrowc, nrowp
    integer :: ntot, lratiox, lratioy
    real(kind=8) :: dividemass, eta1, eta2, etafine, flag, hcnt, hcrse, hfine
    real(kind=8) :: hfineave, hxc, hxf, hyc, hyf, s1, s2, Vnew, xlc, xlp, xrc
    real(kind=8) :: xrp, ybc, ybp, ytc, ytp, hvf, vf
    logical :: set, reloop

    ! Geometry
    nrowp = ihi - ilo + 1
    ncolp = jhi - jlo + 1

    hxf = hxposs(level)
    hyf = hyposs(level)
    xlp = xlower + ilo*hxf
    xrp = xlower + (ihi+1)*hxf
    ybp = ylower + jlo*hyf
    ytp = ylower + (jhi+1)*hyf

    ! Begin by filling values for grids at current level.  If all values can be
    ! filled this way we return
    call intfil(valbig,mitot,mjtot,time,flaguse,nrowst,ncolst,ilo,ihi,jlo,jhi, &
                level,nvar,naux)

    ! Trimbd returns set = true if all of the entries are filled (=1.).
    ! set = false, otherwise. If set = true, then no other levels are
    ! are required to interpolate, and we return.
    !
    ! Note that the used array is filled entirely in intfil, i.e. the
    ! marking done there also takes into account the points filled by
    ! the boundary conditions. bc2amr will be called later, after all 4
    ! boundary pieces filled.
    call trimbd(flaguse,nrowp,ncolp,set,il,ir,jb,jt)

    ! Need to fill in using recursive calls to coarser levels to fill remaining
    ! unset points
    if (.not. set) then

        ! Uh oh, reached end of recursion but there are still cells to fill
        if (level == 1) then
             write(outunit,*)" error in filrecur - level 1 not set"
             write(outunit,"('start at row: ',i4,' col ',i4)") nrowst,ncolst
             print *," error in filrecur - level 1 not set"
             print *," should not need more recursion "
             print *," to set patch boundaries"
             write(*,"('start at row: ',i4,' col ',i4)") nrowst,ncolst
             stop
        end if

        ! set = false. we will have to interpolate some values from coarser
        ! levels. We begin by initializing the level level arrays, so that we 
        ! can use purely recursive formulation for interpolating.
        levc = level - 1
        hxc  = hxposs(levc)
        hyc  = hyposs(levc)

        isl  = il + ilo - 1
        isr  = ir + ilo - 1
        jsb  = jb + jlo - 1
        jst  = jt + jlo - 1
    
        ! Coarsen
        lratiox = intratx(levc)
        lratioy = intraty(levc)
        iplo   = (isl-lratiox  +nghost*lratiox)/lratiox - nghost
        jplo   = (jsb-lratioy  +nghost*lratioy)/lratioy - nghost
        iphi   = (isr+lratiox  )/lratiox
        jphi   = (jst+lratioy  )/lratioy

        xlc  =  xlower + iplo*hxc
        ybc  =  ylower + jplo*hyc
        xrc  =  xlower + (iphi+1)*hxc
        ytc  =  ylower + (jphi+1)*hyc

        nrowc   =  iphi - iplo + 1
        ncolc   =  jphi - jplo + 1
        ntot    = nrowc*ncolc*(nvar+naux)

        ! Check to make sure we have enough scratch space
        ! print "(' needed coarse grid size ',2i5,' allocated ',2i5)", nrowc,ncolc, ihi-ilo+2,jhi-jlo+2
        if (nrowc .gt. ihi-ilo+2 .or. ncolc .gt. jhi-jlo+2) then
            print *," did not make big enough work space in filrecur "
            print *," need coarse space with nrowc,ncolc ",nrowc,ncolc
            print *," coarser level is ",levc," nghost ",nghost
            print *," with ratios lratiox lratioy ",lratiox,lratioy
            print *," made space for ilo,ihi,jlo,jhi ",ilo,ihi,jlo,jhi
            print *," isl,isr,jsb,jst ",isl,isr,jsb,jst
            print *,"iplo,iphi,jplo,jphi ",iplo,iphi,jplo,jphi
            print *," orig il,ir,jb,jt ",il,ir,jb,jt
            print *,"filpatch called with mitot,mjtot ",mitot,mjtot
            call outtre(1,.false.,nvar,naux)
            stop
        endif

        ! Fill in the coarse aux array
        if (naux > 0) then
            call setaux(nrowc - 2*nghost, ncolc - 2*nghost, nghost, &
                        nrowc - 2*nghost, ncolc - 2*nghost, xlc + nghost*hxc, &
                        ybc + nghost*hyc, hxc, hyc, naux, auxcrse)
        end if

        ! If we are in a special domain need to do some stuff, otherwise we
        ! recursively call downwards
        if ((xperdom .or. (yperdom .or. spheredom)) .and. sticksout(iplo,iphi,jplo,jphi)) then
            call prefilrecur(levc,nvar,valcrse,auxcrse,naux,time,nrowc,ncolc, &
                            1,1,iplo,iphi,jplo,jphi)
        else
            call filrecur(levc,nvar,valcrse,auxcrse,naux,time,nrowc,ncolc,1,1, &
                            iplo,iphi,jplo,jphi)
        end if

        ! Interpolate back up

        ! Loop through coarse cells determining interpolation slopes, these will
        ! be saved for the fine grid loop.  Prevents finding the same slope
        ! possibly lratiox * lratioy times.  
        
        ! All fine grid depths will be found before any momentum
        reloop = .false.

        ! Find interpolation slope for eta
        do ic=2,nrowc-1
            do jc=2,ncolc-1
                icount(icrse(ic,jc)) = 0
                finemass(icrse(ic,jc)) = 0.d0
                do ivar=1,nvar
                    fineflag(ivalc(ivar,ic,jc)) = .false.
                end do
                
                ! X-Slope
                do i=-1,1
                    etacrse(icrse(ic+i,jc)) = valcrse(ivalc(1,ic+i,jc)) / rho(1) + auxcrse(iauxc(ic+i,jc))
                    if (valcrse(ivalc(1,ic+i,jc)) / rho(1) < dry_tolerance(1)) then
                        etacrse(icrse(ic+i,jc)) = eta_init(1)
                    end if
                end do
                s1 = etacrse(icrse(ic,jc)) - etacrse(icrse(ic-1,jc))
                s2 = etacrse(icrse(ic+1,jc)) - etacrse(icrse(ic,jc))
                if (s1 * s2 <= 0) then
                    slopex(icrse(ic,jc)) = 0.d0
                else
                    slopex(icrse(ic,jc))= min(abs(s1),abs(s2)) &
                            * sign(1.d0, etacrse(icrse(ic+1,jc)) &
                                       - etacrse(icrse(ic-1,jc)))
                end if

                ! Y-Slope
                do j=-1,1
                    etacrse(icrse(ic,jc+j)) = valcrse(ivalc(1,ic,jc+j)) / rho(1) + auxcrse(iauxc(ic,jc+j))
                    if (valcrse(ivalc(1,ic,jc+j)) / rho(1) < dry_tolerance(1)) then
                        etacrse(icrse(ic,jc+j)) = eta_init(1)
                    end if
                enddo
                s1 = etacrse(icrse(ic,jc))  - etacrse(icrse(ic,jc-1))
                s2 = etacrse(icrse(ic,jc+1))- etacrse(icrse(ic,jc))
                if (s1*s2.le.0) then
                    slopey(icrse(ic,jc))= 0.d0
                else
                    slopey(icrse(ic,jc))=min(abs(s1), abs(s2)) &
                            * sign(1.d0, etacrse(icrse(ic,jc+1)) &
                                       - etacrse(icrse(ic,jc-1)))
                end if
            end do  
        end do

        ! Loop through the patch interpolating surface and then calculating new
        ! depth
        do iff = 1,nrowp
            ic = 2 + (iff - (isl - ilo) - 1) / lratiox
            eta1 = (-0.5d0 + real(mod(iff-1, lratiox), kind=8)) / real(lratiox,kind=8)
            do jf  = 1,ncolp
                jc = 2 + (jf - (jsb - jlo) - 1) / lratioy
                eta2 = (-0.5d0 + real(mod(jf-1, lratioy), kind=8)) / real(lratioy,kind=8)
                cc = icrse(ic,jc)
                flag = flaguse(iff,jf)
                if (flag == 0.0) then
                    ! Interp. from coarse cells to fine grid to find eta
                    icount(icrse(ic,jc)) = icount(icrse(ic,jc)) + 1
                    etafine =  etacrse(icrse(ic,jc)) + eta1 * slopex(icrse(ic,jc)) &
                                                     + eta2 * slopey(icrse(ic,jc))
                    hfine = max(etafine - aux(1,iff+nrowst-1,jf+ncolst-1), 0.d0)
                    valbig(1,iff+nrowst-1,jf+ncolst-1) = hfine * rho(1)

                    ! Note that finemass is actually the mass
                    finemass(icrse(ic,jc)) = finemass(icrse(ic,jc)) + valbig(1,iff+nrowst-1,jf+ncolst-1)
                    if (valbig(1,iff+nrowst-1,jf+ncolst-1) / rho(1) < dry_tolerance(1)) then
                        fineflag(ivalc(1,ic,jc)) = .true.
                        reloop = .true.
                    end if
                end if
            end do
        end do
        ! Done interpolating eta and h

        ! Loop through patch interpolating momentum
        do ivar=2,nvar
            do ic=2,nrowc-1
                do jc=2,ncolc-1

                    ! Determine x-slope
                    s1 = valcrse(ivalc(ivar,ic,jc)) / rho(1) - valcrse(ivalc(ivar,ic-1,jc)) / rho(1)
                    s2 = valcrse(ivalc(ivar,ic+1,jc)) / rho(1) - valcrse(ivalc(ivar,ic,jc)) / rho(1)
                    if (s1 * s2 <= 0) then
                        slopex(icrse(ic,jc)) = 0.d0
                    else
                        slopex(icrse(ic,jc)) = min(abs(s1), abs(s2)) &
                            * sign(1.d0, valcrse(ivalc(ivar,ic+1,jc)) / rho(1) &
                                       - valcrse(ivalc(ivar,ic-1,jc)) / rho(1))
                    end if

                    ! Determine y-slope
                    s1 = valcrse(ivalc(ivar,ic,jc)) / rho(1) - valcrse(ivalc(ivar,ic,jc-1)) / rho(1)
                    s2 = valcrse(ivalc(ivar,ic,jc+1)) / rho(1) - valcrse(ivalc(ivar,ic,jc)) / rho(1)
                    if (s1*s2.le.0) then
                        slopey(icrse(ic,jc)) = 0.d0
                    else
                        slopey(icrse(ic,jc)) = min(abs(s1), abs(s2)) &
                            * sign(1.d0, valcrse(ivalc(ivar,ic,jc+1)) / rho(1) &
                                       - valcrse(ivalc(ivar,ic,jc-1)) / rho(1))
                    endif

                    ! Check for new extrema
                    ! look for bounds on velocity to avoid generating new extrema
                    ! necessary since interpolating momentum linearly
                    ! yet depth is not interpolated linearly
                    if (valcrse(ivalc(1,ic,jc)) / rho(1) > dry_tolerance(1)) then
                        velmax(icrse(ic,jc)) = valcrse(ivalc(ivar,ic,jc)) / valcrse(ivalc(1,ic,jc))
                        velmin(icrse(ic,jc)) = valcrse(ivalc(ivar,ic,jc)) / valcrse(ivalc(1,ic,jc))
                    else
                        velmax(icrse(ic,jc)) = 0.d0
                        velmin(icrse(ic,jc)) = 0.d0
                    end if
                    do ii = -1,1,2
                        if (valcrse(ivalc(1,ic+ii,jc))/ rho(1).gt.dry_tolerance(1)) then
                            velmax(icrse(ic,jc)) = max(velmax(icrse(ic,jc)), &
                                            valcrse(ivalc(ivar,ic+ii,jc)) &
                                          / valcrse(ivalc(1,ic+ii,jc)))
                            velmin(icrse(ic,jc)) = min(velmin(icrse(ic,jc)), &
                                            valcrse(ivalc(ivar,ic+ii,jc)) &
                                          / valcrse(ivalc(1,ic+ii,jc)))
                        end if
                        if (valcrse(ivalc(1,ic,jc+ii))/ rho(1).gt.dry_tolerance(1)) then
                            velmax(icrse(ic,jc)) = max(velmax(icrse(ic,jc)), &
                                            valcrse(ivalc(ivar,ic,jc+ii)) &
                                          / valcrse(ivalc(1,ic,jc+ii)))
                            velmin(icrse(ic,jc)) = min(velmin(icrse(ic,jc)), &
                                            valcrse(ivalc(ivar,ic,jc+ii)) &
                                          / valcrse(ivalc(1,ic,jc+ii)))
                        end if
                    end do

                end do
            end do

            ! Determine momentum in fine cells
            do iff = 1,nrowp
                ic = 2 + (iff - (isl - ilo) - 1) / lratiox
                eta1 = (-0.5d0 + real(mod(iff-1, lratiox),kind=8)) / real(lratiox,kind=8)
                do jf  = 1,ncolp
                    jc = 2 + (jf  - (jsb - jlo) - 1) / lratioy
                    eta2 = (-0.5d0 + real(mod(jf -1,lratioy),kind=8)) / real(lratioy,kind=8)

                    flag = flaguse(iff,jf)
                    if (flag .eq. 0.0) then
                        if (.not.(fineflag(ivalc(1,ic,jc)))) then
                            ! This is a normal wet cell, intepolate normally
                            hvf = valcrse(ivalc(ivar,ic,jc)) / rho(1) &
                                            + eta1 * slopex(icrse(ic,jc)) &
                                            + eta2 * slopey(icrse(ic,jc))
                            ! Here we need to divied by rho since valbig stores
                            ! mass
                    
                            vf = hvf * rho(1) / valbig(1,iff+nrowst-1,jf+ncolst-1)
                            if (vf < velmin(icrse(ic,jc)) .or. &
                                vf > velmax(icrse(ic,jc))) then

                                fineflag(ivalc(ivar,ic,jc)) = .true.
                                reloop = .true.
                            else
                                valbig(ivar,iff+nrowst-1,jf+ncolst-1) = hvf * rho(1)
                            end if
                        end if
                    end if
               end do
            end do

            !loop again and reset momentum to conserve momentum
            !in the case of gained momentum, or if velocity bounds violated
            if (reloop) then
                do iff = 1,nrowp
                    ic = 2 + (iff - (isl - ilo) - 1)/lratiox
                    eta1 = (-0.5d0 + real(mod(iff-1,lratiox),kind=8)) / real(lratiox,kind=8)
                    do jf  = 1,ncolp
                        jc = 2 + (jf  - (jsb - jlo) - 1)/lratioy
                        eta2 = (-0.5d0 + real(mod(jf -1,lratioy),kind=8)) / real(lratioy,kind=8)
                        flag = flaguse(iff,jf)
                        if (flag == 0.0) then
                            if (fineflag(ivalc(1,ic,jc)) .or. fineflag(ivalc(ivar,ic,jc))) then
                                if (finemass(icrse(ic,jc)) / rho(1) > dry_tolerance(1)) then
                                    hcrse = valcrse(ivalc(1,ic,jc)) / rho(1)
                                    hcnt = real(icount(icrse(ic,jc)),kind=8)
                                    hfineave = finemass(icrse(ic,jc)) / (hcnt * rho(1))
                                    dividemass = max(hcrse,hfineave)
                                    hfine = valbig(1,iff+nrowst-1,jf+ncolst-1) / rho(1)
                                    Vnew = valcrse(ivalc(ivar,ic,jc)) / (dividemass * rho(1))
                                    valbig(ivar,iff+nrowst-1,jf+ncolst-1) = Vnew * valbig(1,iff+nrowst-1,jf+ncolst-1)
                                else
                                    valbig(ivar,iff+nrowst-1,jf+ncolst-1) = 0.d0
                                endif
                            endif
                        endif
                    enddo
                enddo
            endif
        end do
        ! Done interpolating momenta and other conserved quantities

    end if

    !
    !  set bcs, whether or not recursive calls needed. set any part of patch that stuck out
    !
    call bc2amr(valbig,aux,mitot,mjtot,nvar,naux,hxf,hyf,level,time, &
                xlp,xrp,ybp,ytp,xlower,ylower,xupper,yupper, &
                xperdom,yperdom,spheredom)

contains

    ! Index of variable requested in scratch space
    integer pure function ivalc(ivar,i,j)
        implicit none
        integer, intent(in) :: ivar,i,j
        ivalc = ivar + nvar*(i-1) + nvar*nrowc*(j-1)
    end function ivalc

    ! Index into coarse arrays
    integer pure function icrse(i,j)
        implicit none
        integer, intent(in) :: i,j
        icrse = i + nrowc*(j-1)
    end function icrse

    ! Index into aux coarse arrays
    integer pure function iauxc(i,j)
        implicit none
        integer, intent(in) :: i,j
        iauxc = 1 + naux*(i-1) + naux*nrowc*(j-1)
    end function iauxc

    ! Check if patch sticks out of coarse patch (I think)
    logical pure function sticksout(iplo,iphi,jplo,jphi)
        implicit none
        integer, intent(in) :: iplo,iphi,jplo,jphi
        sticksout = (iplo < 0 .or. &
                     jplo < 0 .or. &
                     iphi >= iregsz(levc) .or. &
                     jphi >= jregsz(levc))
    end function sticksout
    
end subroutine filrecur