! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! ::::: Parameters, variables, subroutines related to gauges
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

! Contains:
!   subroutine set_gauges
!     Called initially to read from gauges.data
!   subroutine setbestsrc
!     Called each time regridding is done to determine which patch to 
!     use for interpolating to each gauge location.
!   subroutine print_gauges
!     Called each time step for each grid patch.
!     Refactored dumpgauge routine to interpolate for all gauges on patch.
!
!     Note: by default all components of q are printed at each gauge.
!     To print something different or a different precision, modify 
!     format statement 100 and/or the write statement that uses it.
!   
! Note: Updated for Clawpack 5.3.0:
!   - the dumpgauge and setbestsrc subroutines have been moved to this module 
!     and the dumpgauge subroutine has been refactored and renamed print_gauges.
!   - dumpgauge.f must be removed from Makefiles.
!   - setbestsrc uses quicksort to sort gauge numbers and
!     then figures out which gauges will be updated by grid, and stores this
!     information in new module variables mbestg1, mbestg2.
!   - print_gauges no longer uses binary search to locate first gauge handled
!     by a grid.  Instead loop over gauges specified by mbestg1, mbestg2.

module gauges_module

    implicit none
    save

    logical, private :: module_setup = .false.

    integer, parameter :: OUTGAUGEUNIT=89
    integer :: num_gauges, inum
    real(kind=8), allocatable :: xgauge(:), ygauge(:), t1gauge(:), t2gauge(:)
    integer, allocatable, dimension(:) ::  mbestsrc, mbestorder, &
                          igauge, mbestg1, mbestg2, nextLoc
    
!    integer, parameter :: MAXDATA=1
    integer, parameter :: MAXDATA=1000
    real(kind=8), pointer :: gaugeArray(:,:,:)
    integer, pointer :: levelArray(:,:)

contains

    subroutine set_gauges(restart, nvar, fname)

        use amr_module
        use multilayer_module, only: num_layers

        implicit none

        ! Input
        character(len=*), intent(in), optional :: fname
        logical, intent(in) :: restart
        integer, intent(in) :: nvar

        ! Locals
        integer :: i, ipos, idigit
        integer, parameter :: iunit = 7
        character(len=14) ::  fileName

        if (.not.module_setup) then
            ! Open file
            if (present(fname)) then
                call opendatafile(iunit,fname)
            else
                call opendatafile(iunit,'gauges.data')
            endif

            read(iunit,*) num_gauges

            allocate(xgauge(num_gauges), ygauge(num_gauges))
            allocate(t1gauge(num_gauges), t2gauge(num_gauges))
            allocate(mbestsrc(num_gauges), mbestorder(num_gauges))
            allocate(igauge(num_gauges))
            allocate(mbestg1(maxgr), mbestg2(maxgr))

            allocate(nextLoc(num_gauges))
            print *, 1 + num_layers, nvar, MAXDATA, num_gauges
            allocate(gaugeArray(1 + num_layers + nvar, MAXDATA, num_gauges))
            allocate(levelArray(MAXDATA, num_gauges))
            
            do i=1,num_gauges
                read(iunit,*) igauge(i),xgauge(i),ygauge(i),t1gauge(i),t2gauge(i)
            enddo

            close(iunit)
            
            ! initialize for starters
            mbestsrc = 0
            nextLoc  = 1  ! next location to be filled with gauge info

            do i = 1, num_gauges
               fileName = 'gaugexxxxx.txt'    ! NB different name convention too
               inum = igauge(i)
               do ipos = 10,6,-1              ! do this to replace the xxxxx in the name
                  idigit = mod(inum,10)
                  fileName(ipos:ipos) = char(ichar('0') + idigit)
                  inum = inum / 10
               end do

    !          status unknown since might be a restart run. maybe need to test and rewind?
               if (restart) then
                  open(unit=OUTGAUGEUNIT, file=fileName, status='old',        &
                       position='append', form='formatted')
               else
                  open(unit=OUTGAUGEUNIT, file=fileName, status='unknown',        &
                       position='append', form='formatted')
                  rewind OUTGAUGEUNIT
                  write(OUTGAUGEUNIT,100) igauge(i), xgauge(i), ygauge(i), 8
 100              format("# gauge_id= ",i5," location=( ",1e15.7," ",1e15.7," ) num_eqn= ",i2)
                  write(OUTGAUGEUNIT,101)
 101              format("# Columns: level time h(1) hu(1) hv(1) eta(1) h(2) hu(2), hv(2) eta(2)")
                endif

                close(OUTGAUGEUNIT)

            end do
            module_setup = .true.
        end if

    end subroutine set_gauges


!
! --------------------------------------------------------------------
!
      subroutine setbestsrc()
!
!     Called every time grids change, to set the best source grid patch
!     for each gauge, i.e. the finest level patch that includes the gauge.
!
!     lbase is grid level that didn't change, but since fine
!     grid may have disappeared, we still have to look starting
!     at coarsest level 1.
!
      use amr_module
      implicit none

      integer :: lev, mptr, i, k1, ki

!
! ##  set source grid for each loc from coarsest level to finest.
! ##  that way finest src grid left and old ones overwritten
! ##  this code uses fact that grids do not overlap

! # for debugging, initialize sources to 0 then check that all set
      do i = 1, num_gauges
         mbestsrc(i) = 0
      end do

 
      do 20 lev = 1, lfine  
          mptr = lstart(lev)
 5        do 10 i = 1, num_gauges
            if ((xgauge(i) .ge. rnode(cornxlo,mptr)) .and. &
                (xgauge(i) .le. rnode(cornxhi,mptr)) .and. &  
                (ygauge(i) .ge. rnode(cornylo,mptr)) .and. &
                (ygauge(i) .le. rnode(cornyhi,mptr)) ) then
               mbestsrc(i) = mptr
            endif
 10       continue

          mptr = node(levelptr, mptr)
          if (mptr .ne. 0) go to 5
 20   continue


      do i = 1, num_gauges
        if (mbestsrc(i) .eq. 0) &
            write(6,*)"ERROR in setting grid src for gauge data",i
      end do

!     Sort the source arrays for easy testing during integration
      call qsorti(mbestorder,num_gauges,mbestsrc)

!     After sorting,  
!           mbestsrc(mbestorder(i)) = grid index to be used for gauge i
!     and mbestsrc(mbestorder(i)) is non-decreasing as i=1,2,..., num_gauges

!     write(6,*) '+++ mbestorder: ',mbestorder
!     write(6,*) '+++ mbestsrc: ',mbestsrc

!     Figure out the set of gauges that should be handled on each grid:  
!     after loop below, grid k should handle gauges numbered
!          mbestorder(i) for i = mbestg1(k), mbestg1(k)+1, ..., mbestg2(k)
!     This will be used for looping in print_gauges subroutine.

      ! initialize arrays to default indicating grids that contain no gauges:
      mbestg1 = 0
      mbestg2 = 0

      k1 = 0
      do i=1,num_gauges
          ki = mbestsrc(mbestorder(i))
          if (ki > k1) then
              ! new grid number seen for first time in list
              if (k1 > 0) then
                  ! mark end of gauges seen by previous grid
                  mbestg2(k1) = i-1
!                 write(6,*) '+++ k1, mbestg2(k1): ',k1,mbestg2(k1)
                  endif
              mbestg1(ki) = i
!             write(6,*) '+++ ki, mbestg1(ki): ',ki,mbestg1(ki)
              endif
          k1 = ki
          enddo
      if (num_gauges > 0) then
          ! finalize 
          mbestg2(ki) = num_gauges
!         write(6,*) '+++ ki, mbestg2(ki): ',ki,mbestg2(ki)
          endif


    end subroutine setbestsrc

    ! ==========================================================================

    subroutine update_gauges(q,aux,xlow,ylow,nvar,mitot,mjtot,naux,mptr)

        ! This routine is called each time step for each grid patch, to output
        ! gauge values for all gauges for which this patch is the best one to
        ! use (i.e. at the finest refinement level).
        
        ! It is called after ghost cells have been filled from adjacent grids
        ! at the same level, so bilinear interpolation can be used to to
        ! compute values at any gauge location that is covered by this grid.
        
        ! The grid patch is designated by mptr. We only want to set gauges i
        ! for which mbestsrc(i) == mptr. The array mbestsrc is reset after
        ! each regridding to indicate which grid patch is best to use for each
        ! gauge.
        
        ! This is a refactoring of dumpgauge.f from Clawpack 5.2 Loops over
        ! only the gauges to be handled by this grid, as specified by indices
        ! from mbestg1(mptr) to mbestg2(mptr)
        
        use amr_module
        use multilayer_module, only: num_layers, dry_tolerance, rho
        
        implicit none
        
        integer, intent(in) ::  nvar,mitot,mjtot,naux,mptr
        real(kind=8), intent(in) ::  q(nvar,mitot,mjtot)
        real(kind=8), intent(in) ::  aux(naux,mitot,mjtot)
        real(kind=8), intent(in) ::  xlow,ylow
        
        ! local variables:
        real(kind=8) :: var(maxvar)
        real(kind=8) :: xcent,ycent,xoff,yoff,tgrid,hx,hy
        integer :: level,i,j,iindex,jindex,ivar, ii,i1,i2
        real(kind=8) :: drytol2, topo
        integer :: icell,jcell, index
        integer :: m, k, layer_index
        real(kind=8) :: eta(num_layers), h(num_layers, 4)

        if (num_gauges == 0) then
           return
        endif
        i1 = mbestg1(mptr)
        i2 = mbestg2(mptr)
        if (i1 == 0) then
           ! no gauges to be handled by this grid
           return
        endif

        ! this stuff the same for all gauges on this grid
        tgrid = rnode(timemult,mptr)
        level = node(nestlevel,mptr)
        hx    =  hxposs(level)
        hy    =  hyposs(level)

        do i = i1,i2
            ii = mbestorder(i)
            if (mptr .ne. mbestsrc(ii)) then !!! go to 10  ! this patch not used
                print *, '*** should not happen... i, ii, mbestsrc(ii), mptr:'
                print *, i, ii, mbestsrc(ii), mptr
                stop
            end if
            if (tgrid < t1gauge(ii) .or. tgrid > t2gauge(ii)) then
                ! don't output at this time for gauge i
                cycle
            end if
    !
    !    ## if we did not skip to line 10, we need to output gauge i:
    !    ## prepare to do bilinear interp at gauge location to get vars
    !
    !    *** Note: changed 0.5 to  0.5d0 etc. ****************************
    !
            iindex =  int(.5d0 + (xgauge(ii)-xlow)/hx)
            jindex =  int(.5d0 + (ygauge(ii)-ylow)/hy)
            if ((iindex .lt. nghost .or. iindex .gt. mitot-nghost) .or. &
                (jindex .lt. nghost .or. jindex .gt. mjtot-nghost)) &
              write(*,*)"ERROR in output of Gauge Data "
            xcent  = xlow + (iindex-.5d0)*hx
            ycent  = ylow + (jindex-.5d0)*hy
            xoff   = (xgauge(ii)-xcent)/hx
            yoff   = (ygauge(ii)-ycent)/hy
    !  IF WANT TO USE, MODIFY TO TEST FOR ROUNDOFF LEVEL DIFF
    !       if (xoff .lt. 0.d0 .or. xoff .gt. 1.d0 .or. &
    !           yoff .lt. 0.d0 .or. yoff .gt. 1.d0) then
    !          write(6,*)" BIG PROBLEM in DUMPGAUGE", i
    !       endif

         ! ## Modified below from amrclaw/src/2d/gauges_module.f90 
         ! ## to interpolate only where all four cells are
         ! ## wet, otherwise just take this cell value:

         ! Check for dry cells by comparing h to drytol2, which should be smaller
         ! than drytolerance to avoid oscillations since when h < drytolerance the
         ! velocities are zeroed out which can then lead to increase in h again.

            do m = 1, num_layers
                layer_index = 3 * (m - 1)
                drytol2 = 0.1d0 * dry_tolerance(m)

                h(m, 1) = q(1,iindex,jindex) / rho(m)
                h(m, 2) = q(1,iindex+1,jindex) / rho(m)
                h(m, 3) = q(1,iindex,jindex+1) / rho(m)
                h(m, 4) = q(1,iindex+1,jindex+1) / rho(m)
            
                if ((h(m, 1) < drytol2) .or.  &
                    (h(m, 2) < drytol2) .or.  &
                    (h(m, 3) < drytol2) .or.  &
                    (h(m, 4) < drytol2)) then
                    ! One of the cells is dry, so just use value from grid cell
                    ! that contains gauge rather than interpolating
                    
                    icell = int(1.d0 + (xgauge(ii) - xlow) / hx)
                    jcell = int(1.d0 + (ygauge(ii) - ylow) / hy)
                    do ivar=1,3
                        var(ivar + layer_index) =                    &
                                    q(ivar + layer_index, icell, jcell) / rho(m)
                    enddo
                    ! This is the bottom layer and we should figure out the
                    ! topography
                    if (m == num_layers) then
                        topo = aux(1, icell, jcell)
                    end if
                else
                    ! Linear interpolation between four cells
                    do ivar=1,3
                        var(ivar + layer_index) = &
                           (1.d0 - xoff) * (1.d0 - yoff) * q(ivar + layer_index,     iindex,     jindex) / rho(m)  &
                         + (       xoff) * (1.d0 - yoff) * q(ivar + layer_index, 1 + iindex,     jindex) / rho(m)  &
                         + (1.d0 - xoff) * (       yoff) * q(ivar + layer_index,     iindex, 1 + jindex) / rho(m)  &
                         + (       xoff) * (       yoff) * q(ivar + layer_index, 1 + iindex, 1 + jindex) / rho(m)
                    enddo
                    topo = (1.d0 - xoff) * (1.d0 - yoff) * aux(1,     iindex,     jindex)  &
                         + (       xoff) * (1.d0 - yoff) * aux(1, 1 + iindex,     jindex)  &
                         + (1.d0 - xoff) * (       yoff) * aux(1,     iindex, 1 + jindex)  &
                         + (       xoff) * (       yoff) * aux(1, 1 + iindex, 1 + jindex)
                endif
            end do

            ! Extract surfaces
            eta(num_layers) = var(3 * num_layers - 2) + topo
            do k = num_layers - 1, 1, -1
                eta(k) = var(3 * k - 2) + eta(k + 1)
                if (abs(eta(k)) < 1d-90) eta(k) = 0.d0
            end do

            ! Zero out tiny values to prevent later problems reading data,
            ! as done in valout.f
            do j = 1,3 * num_layers
               if (abs(var(j)) < 1d-90) var(j) = 0.d0
            end do
       
            ! save info for this time 
            index = nextLoc(ii)
     
            levelArray(index,ii) = level
            gaugeArray(1,index,ii) = tgrid
            gaugeArray(2,index,ii) = var(1)
            gaugeArray(3,index,ii) = var(2)
            gaugeArray(4,index,ii) = var(3)
            gaugeArray(5,index,ii) = eta(1)
            gaugeArray(6,index,ii) = var(4)
            gaugeArray(7,index,ii) = var(5)
            gaugeArray(8,index,ii) = var(6)
            gaugeArray(9,index,ii) = eta(2)
            
            nextLoc(ii) = nextLoc(ii) + 1
            if (nextLoc(ii) .gt. MAXDATA) then
              call print_gauges_and_reset_nextLoc(ii, nvar)  
            endif

        enddo ! End of gauge loop
 
    end subroutine update_gauges

    ! ==========================================================================

    subroutine print_gauges_and_reset_nextLoc(gaugeNum, nvar)
    !
    ! Array of gauge data for this gauge reached max capacity
    ! print to file.

        use multilayer_module, only: num_layers

        implicit none
      
        ! Locals
        integer :: gaugeNum, j, inum, k, idigit, ipos, myunit, nvar
        character(len=14) :: fileName
        integer :: mythread, num_values
        integer :: omp_get_thread_num
        character(len=15), parameter :: GAUGE_OUTPUT_FORMAT = "(i5.2, 10e15.7)"
        ! open file for gauge gaugeNum, figure out name
        ! not writing gauge number since it is in file name now
        ! status is old, since file name and info written for
        ! each file in in set_gauges.
        !
        ! NB: output written in different order, losing backward compatibility

        ! Construct file name
        fileName = 'gaugexxxxx.txt'
        inum = igauge(gaugeNum)
        do ipos = 10, 6, -1           ! do this to replace the xxxxx in the name
            idigit = mod(inum, 10)
            fileName(ipos:ipos) = char(ichar('0') + idigit)
            inum = inum / 10
        end do


        mythread = 0
        ! Following line indexes the thread number and determines a output unit
        ! based on this number, if 
!$      mythread = omp_get_thread_num()
        myunit = OUTGAUGEUNIT + mythread

        ! add thread number of outgaugeunit to make a unique unit number.
        ! ok since writing to a unique file. in serial, still using only 
        ! IOUTGAUGEUNIT
        open(unit=myunit, file=fileName, status='old',    &
             position='append', form='formatted')
        
        ! called either because array is full (write MAXDATA amount of gauge 
        ! data) or checkpoint time, so write whatever is in array and reset.
        ! nextLoc has already been increment before this subr. called
        num_values = 1 + nvar + num_layers
        do j = 1, nextLoc(gaugeNum)-1
            write(myunit, GAUGE_OUTPUT_FORMAT) levelArray(j,gaugeNum),    &
                                  (gaugeArray(k, j, gaugeNum), k=1, num_values)
        end do
        nextLoc(gaugeNum) = 1                        
  
        ! close file
        close(myunit)
  
    end subroutine print_gauges_and_reset_nextLoc
  
end module gauges_module
