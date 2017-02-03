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
! ------------------------------------------------------------------
!
subroutine filval(val, mitot, mjtot, dx, dy, level, time,  mic, &
                  mjc, xleft, xright, ybot, ytop, nvar, mptr, ilo, ihi, &
                  jlo, jhi, aux, naux,  sp_over_h )

    use amr_module, only: xlower, ylower, intratx, intraty, nghost, xperdom
    use amr_module, only: yperdom, spheredom, xupper, yupper, alloc
    use amr_module, only: outunit, NEEDS_TO_BE_SET
    use amr_module, only: newstl, iregsz, jregsz

    use topo_module, only: aux_finalized
    use geoclaw_module, only: dry_tolerance, sea_level
    use refinement_module, only: varRefTime
    use topo_module, only: aux_finalized

    implicit none

    ! Input
    integer, intent(in) :: mitot, mjtot, level, mic, mjc, nvar, mptr, ilo, ihi
    integer, intent(in) :: jlo, jhi, naux
    real(kind=8), intent(in) :: dx, dy, time, xleft, xright, ybot, ytop

    ! Output
    real(kind=8), intent(in out) :: sp_over_h
    real(kind=8), intent(in out) :: val(nvar,mitot,mjtot), aux(naux,mitot,mjtot)

    ! Local storage
    integer :: refinement_ratio_x, refinement_ratio_y, iclo, jclo, ichi, jchi, ng, i, ico, ifine
    integer :: ii, ivar, j, jco, jfine, jj
    real(kind=8) :: valc(nvar,mic,mjc), auxc(naux,mic,mjc)
    real(kind=8) :: coarseval(3), dx_coarse, dy_coarse, xl, xr, yb, yt, area
    real(kind=8) :: dividemass, finemass, hvf, s1m, s1p, slopex, slopey, vel
    real(kind=8) :: velmax, velmin, vf, vnew, xoff, yoff
    logical :: fineflag(3)
    real(kind=8) :: fliparray((mitot+mjtot)*(nvar+naux))
    real(kind=8) :: aux2(naux,mitot,mjtot)
    integer :: nx, ny
    real(kind=8) setflags(mitot,mjtot),maxauxdif
    integer :: jm, im, nm
    logical :: sticksoutxfine, sticksoutyfine,sticksoutxcrse,sticksoutycrse
    logical :: DIAGONAL_CORNER

    ! External function definitions
    real(kind=8) :: get_max_speed

    DIAGONAL_CORNER(im,jm,mic,mjc) = (im .eq. 1   .and. jm .eq. mjc)  .or.     &
                                      (im .eq. mic .and. jm .eq. mjc) .or.     &
                                      (im .eq. 1   .and. jm .eq. 1)    .or.     &
                                      (im .eq. mic .and. jm .eq. 1)

    refinement_ratio_x = intratx(level-1)
    refinement_ratio_y = intraty(level-1)
    dx_coarse  = dx * refinement_ratio_x
    dy_coarse  = dy * refinement_ratio_y
    xl      = xleft  - dx_coarse
    xr      = xright + dx_coarse
    yb      = ybot   - dy_coarse
    yt      = ytop   + dy_coarse

    ! if topo not yet final then aux is set outside filval (in gfixup)
    ! and so aux has real data already, (ie dont overwrite here)

    ! set integer indices for coarser patch enlarged by 1 cell
    ! (can stick out of domain). proper nesting will insure this one
    ! call is sufficient.
    iclo   = ilo / refinement_ratio_x - 1
    jclo   = jlo / refinement_ratio_y - 1
    ichi   = (ihi + 1) / refinement_ratio_x - 1 + 1
    jchi   = (jhi + 1) / refinement_ratio_y - 1 + 1
    ng     = 0


    sticksoutxfine = ( (ilo .lt. 0) .or. (ihi .ge. iregsz(level)))
    sticksoutyfine = ( (jlo .lt. 0) .or. (jhi .ge. jregsz(level)))
    sticksoutxcrse = ((iclo .lt. 0) .or. (ichi .ge. iregsz(level-1)))
    sticksoutycrse = ((jclo .lt. 0) .or. (jchi .ge. jregsz(level-1)))


    if (naux == 0) then
        write(*,*)" in filval/geoclaw with naux=0:  how could this happen?"
        if ((xperdom .and. sticksoutxcrse) .or. (yperdom.and. sticksoutycrse) .or. spheredom) then
            call preintcopy(valc,mic,mjc,nvar,iclo,ichi,jclo,jchi,level-1,fliparray)
        else
            call intcopy(valc,mic,mjc,nvar,iclo,ichi,jclo,jchi,level-1,1,1)
        endif
    else  
        ! intersect grids and copy all (soln and aux)
        if ((xperdom .and. sticksoutxcrse) .or. (yperdom.and. sticksoutycrse) .or. spheredom) then
            call preicall(valc,auxc,mic,mjc,nvar,naux,iclo,ichi,jclo,jchi, &
                          level-1,fliparray)
        else
            call icall(valc,auxc,mic,mjc,nvar,naux,iclo,ichi,jclo,jchi,level-1,1,1)
        endif
        if (aux_finalized .lt. 2) then ! coarse topo was at wrong time. redo
           ! no ghost cells on coarse enlarged patch
           auxc(1,:,:) = NEEDS_TO_BE_SET  ! needs signal for setaux, set everywhere
           call setaux(ng,mic,mjc,xl,yb,dx_coarse,dy_coarse,naux,auxc)
        endif

        if (aux_finalized < 2) then
            call setaux(0,mic,mjc,xl,yb,dx_coarse,dy_coarse,naux,auxc)
        endif
    endif

    call bc2amr(valc,auxc,mic,mjc,nvar,naux,dx_coarse,dy_coarse,level-1,time,  &
                xl,xr,yb,yt)


!  NOTE change in order of code.  Since the interp from coarse to fine needs the aux
!       arrays set already, the fine copy is done first, to set up the aux arrays.
!       we can do this since we have the flag array to test where to overwrite.

!  SO this is no longer overwriting but setting for the first time.
! overwrite interpolated values with fine grid values, if available.
! can only do this if topo stoped moving, otherwise fine grid
! topo is at previous time step.
!! also might need preicallCopy???

       nx = mitot - 2*nghost
       ny = mjtot - 2*nghost

       if (naux .gt. 0) then 
             aux(1,:,:) = NEEDS_TO_BE_SET  ! will indicate fine cells not yet set
             if ((xperdom .and. sticksoutxfine)  .or. (yperdom.and.sticksoutyfine)) then
                call preicall(val,aux,mitot,mjtot,nvar,naux,ilo-nghost,ihi+nghost,  &
                              jlo-nghost,jhi+nghost,level,fliparray)  
             else
                call icall(val,aux,mitot,mjtot,nvar,naux,ilo-nghost,ihi+nghost,  &
                           jlo-nghost,jhi+nghost,level,1,1)   
             endif
             setflags = aux(1,:,:)   ! save since will overwrite in setaux when setting all aux vals
             ! need this so we know where to use coarse grid to set fine solution w/o overwriting
             if (aux_finalized .lt. 2) aux(1,:,:) = NEEDS_TO_BE_SET  ! reset entire aux array since topo moving
               !set remaining aux vals not set by copying from prev existing grids
               call setaux(nghost,nx,ny,xleft,ybot,dx,dy,naux,aux)
       else ! either no aux exists, or cant reuse yet  
          ! if topo not final, then setaux called in gfixup before this routine
          ! so only call intcopy (which copies soln) and not icall.
          if ((xperdom .and. sticksoutxfine)  .or. (yperdom.and.sticksoutyfine)) then
             call preintcopy(val,mitot,mjtot,nvar,ilo-nghost,ihi+nghost,  &
                       jlo-nghost,jhi+nghost,level,1,1,fliparray)
          else
             call intcopy(val,mitot,mjtot,nvar,ilo-nghost,ihi+nghost,  &
                          jlo-nghost,jhi+nghost,level,1,1)   
          endif
       endif
  
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
                coarseval(2+ii) = valc(1,i+ii,j)  + auxc(1,i+ii,j)
                if (valc(1,i+ii,j)  <= dry_tolerance) then
                    coarseval(2+ii)=sea_level
                end if
            end do
            s1p = coarseval(3) - coarseval(2)
            s1m = coarseval(2) - coarseval(1)
            slopex = min(abs(s1p), abs(s1m)) &
                                * sign(1.d0,coarseval(3) - coarseval(1))
            if (s1m*s1p <= 0.d0) slopex=0.d0

            do jj=-1,1
                coarseval(2+jj) = valc(1,i,j+jj) + auxc(1,i,j+jj)
                if (valc(1,i,j+jj) <= dry_tolerance) then
                    coarseval(2+jj)=sea_level
                end if
            end do
            s1p = coarseval(3) - coarseval(2)
            s1m = coarseval(2) - coarseval(1)
            slopey = min(abs(s1p), abs(s1m)) &
                                * sign(1.d0,coarseval(3)-coarseval(1))
            if (s1m*s1p <= 0.d0) slopey=0.d0

            ! Interpolate from coarse cells to fine grid to find depth
            finemass = 0.d0
                do jco = 1,refinement_ratio_y
               do ico = 1,refinement_ratio_x
                    yoff = (real(jco,kind=8) - 0.5d0) / refinement_ratio_y - 0.5d0
                    xoff = (real(ico,kind=8) - 0.5d0) / refinement_ratio_x - 0.5d0
                    jfine = (j-2) * refinement_ratio_y + nghost + jco
                    ifine = (i-2) * refinement_ratio_x + nghost + ico
                    if (setflags(ifine,jfine) .eq. NEEDS_TO_BE_SET) then
                       val(1,ifine,jfine) = (coarseval(2) + xoff * slopex &
                                                          + yoff * slopey)
                       val(1,ifine,jfine) = max(0.d0, val(1,ifine,jfine)  &
                                               - aux(1,ifine,jfine))
                       finemass = finemass + val(1,ifine,jfine)
                       if (val(1,ifine,jfine) <= dry_tolerance) then
                          fineflag(1) = .true.
                          val(2,ifine,jfine) = 0.d0
                          val(3,ifine,jfine) = 0.d0
                       endif
                    endif
                end do
            end do

            !------ Determine Momentum ----------------------------------
            ! finemass is the total mass in all new fine grid cells
            ! all fine mass has been determined for this coarse grid cell
            ! if all fine cells are dry, momentum has already been set
            if (finemass >= dry_tolerance) then
                do ivar = 2,nvar
                    fineflag(ivar)=.false.
                    s1p = (valc(ivar,i+1,j) - valc(ivar,i,j))
                    s1m = (valc(ivar,i,j) - valc(ivar,i-1,j))
                    slopex = min(abs(s1p), abs(s1m)) &
                     * sign(1.d0,(valc(ivar,i+1,j) - valc(ivar,i-1,j)))
                    if (s1m*s1p.le.0.d0) slopex=0.d0
                
                    s1p = (valc(ivar,i,j+1) - valc(ivar,i,j))
                    s1m = (valc(ivar,i,j) - valc(ivar,i,j-1))
                    slopey = min(abs(s1p), abs(s1m)) &
                     * sign(1.d0,(valc(ivar,i,j+1) - valc(ivar,i,j-1)))
                    if (s1m*s1p <= 0.d0) slopey=0.d0

                    if (valc(1,i,j) > dry_tolerance) then
                        velmax = valc(ivar,i,j) / valc(1,i,j)
                        velmin = valc(ivar,i,j) / valc(1,i,j)
                    else
                        velmax = 0.d0
                        velmin = 0.d0
                    endif
               
                    do ii = -1,1,2
                        if (valc(1,i+ii,j) > dry_tolerance) then
                            vel = valc(ivar,i+ii,j) / valc(1,i+ii,j)
                            velmax = max(vel,velmax)
                            velmin = min(vel,velmin)
                        endif
                        if (valc(1,i,j+ii) > dry_tolerance) then
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
                            hvf = valc(ivar,i,j) + xoff * slopex &
                                                          + yoff*slopey
                            vf = hvf / (val(1,ifine,jfine))
                            if (vf > velmax .or. vf < velmin) then
                                fineflag(ivar) = .true.
                                exit
                            else
                                val(ivar,ifine,jfine) = hvf
                            endif
                        enddo
                    enddo

                    ! momentum is set to preserve old momentum or not violate
                    ! generating new extrema in velocities
                    if (fineflag(1) .or. fineflag(ivar)) then
                        ! more mass now, conserve momentum
                        area = real(refinement_ratio_x * refinement_ratio_y,kind=8)
                        dividemass = max(finemass,valc(1,i,j))
                        Vnew = area * valc(ivar,i,j) / (dividemass)

                            do jco = 1,refinement_ratio_y
                            do ico = 1,refinement_ratio_x
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

    ! scan for max wave speed on newly created grid. this will be used to set appropriate
    ! time step and appropriate refinement in time. For this app not nec to refine by same
    ! amount as refinement in space since refinement at shores where h is shallow has lower
    ! speeds.

! CHECK BY CALLING SETAUX AND SETTING ALL, THEN DIFFING
!   aux2(1,:,:) = NEEDS_TO_BE_SET   ! indicates fine cells not yet set
!   call setaux(nghost,nx,ny,xleft,ybot,dx,dy,naux,aux2)
!   maxauxdif = 1.d-13
!   do i = 1, mitot
!   do j = 1, mjtot
!     if (abs(aux(1,i,j)-aux2(1,i,j)) .gt. maxauxdif) then
!        maxauxdif = abs(aux(1,i,j)-aux2(1,i,j))
!        write(*,444)i,j,aux(1,i,j),aux2(1,i,j),maxauxdif
!444     format("i,j = ",2i4," auxs ",2e15.7," maxauxdif ",e12.5)
!     endif
!   end do
!   end do
!   if (maxauxdif .gt. 2.d-13) then
!      write(*,*)" maxauxdif = ",maxauxdif," with mitot,mjtot ",mitot,mjtot, &
!                " on grid ",mptr," level ",level
!   endif
    

    if (varRefTime) then   ! keep consistent with setgrd_geo and qinit_geo
        sp_over_h = get_max_speed(val,mitot,mjtot,nvar,aux,naux,nghost,dx,dy)
    endif

end subroutine filval

! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine dumpaux(aux,naux,mitot,mjtot)
   implicit none
   real(kind=8) :: aux(naux,mitot,mjtot)
   integer :: naux,mitot,mjtot,i,j,iaux

   do j = 1, mjtot 
   do i = 1, mitot 
      write(*,444) i,j,(aux(iaux,i,j),iaux=1,naux)
 444  format(2i4,5e12.5)
   end do
   end do

end subroutine dumpaux
