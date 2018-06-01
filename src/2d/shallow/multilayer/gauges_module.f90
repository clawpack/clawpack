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
!
! Note: Updated for Clawpack 5.4.0
!   - refactor so each gauge writes to its own file, and batches the writes instead of 
!     writing one at a time. This will remove the critical section and should speed up gauges a lot
!   - When array is filled, that gauge will write to file and start over. 
!   - Need to save index so know position in array where left off
!   - At checkpoint times, dump all gauges
!
! Note: Updated for Clawpack 5.4.x
!  - Add gauge formatting capabilities
!
! Changed mjb, May, 2018 hackathon, to remove need for mbestg1 and mbestg2.
! looked at since they depended on maxgr, which is now a variable, not constant.
! also, algorithms didn't scale for O(10^5) grids, and O(100) gauges..

module gauges_module

    implicit none
    save

    logical, private :: module_setup = .false.

    integer, parameter :: OUTGAUGEUNIT=89
    integer :: num_gauges

    integer, parameter :: MAX_BUFFER = 1000

    ! Gauge data types
    type gauge_type
        ! Gauge number
        integer :: gauge_num

        character(len=14) :: file_name

        ! Location in time and space
        real(kind=8) :: x, y, t_start, t_end

        ! Last time recorded
        real(kind=8) :: last_time

        ! Output settings
        integer :: file_format
        real(kind=8) :: min_time_increment
        character(len=10) :: display_format
        logical, allocatable :: q_out_vars(:)
        logical, allocatable :: aux_out_vars(:)
        integer :: num_out_vars

        ! Data buffers - data holds output and time
        real(kind=8), allocatable :: data(:, :)
        integer :: level(MAX_BUFFER)

        ! Where we are in the buffer
        integer :: buffer_index
    end type gauge_type

    ! Gague array
    type(gauge_type), allocatable :: gauges(:)

    ! Gauge source info
    integer, allocatable, dimension(:) ::  mbestsrc, igauge

contains

    subroutine set_gauges(restart, num_eqn, num_aux, fname)

        use utility_module, only: get_value_count

        use multilayer_module, only: num_layers

        implicit none

        ! Input
        logical, intent(in) :: restart
        integer :: num_eqn, num_aux
        character(len=*), intent(in), optional :: fname

        ! Locals
        integer :: i, n, index
        integer :: num, pos, digit
        integer, parameter :: UNIT = 7
        character(len=128) :: header_1
        character(len=32) :: q_column, aux_column, eta_column

        if (.not.module_setup) then

            ! Open file
            if (present(fname)) then
                call opendatafile(UNIT,fname)
            else
                call opendatafile(UNIT,'gauges.data')
            endif

            read(UNIT,*) num_gauges
            allocate(gauges(num_gauges))
            
            ! Initialize gauge source data
            allocate(mbestsrc(num_gauges))
            mbestsrc = 0

            ! Original gauge information
            do i=1,num_gauges
                read(UNIT, *) gauges(i)%gauge_num, gauges(i)%x, gauges(i)%y, &
                              gauges(i)%t_start, gauges(i)%t_end
                gauges(i)%buffer_index = 1
                gauges(i)%last_time = gauges(i)%t_start
            enddo

            ! Read in output formats
            read(UNIT, *)
            read(UNIT, *)
            read(UNIT, *) (gauges(i)%file_format, i=1, num_gauges)
            read(UNIT, *)
            read(UNIT, *)
            read(UNIT, *) (gauges(i)%display_format, i=1, num_gauges)
            read(UNIT, *)
            read(UNIT, *)
            read(UNIT, *) (gauges(i)%min_time_increment, i=1, num_gauges)

            ! Read in q fields
            read(UNIT, *)
            read(UNIT, *)
            do i = 1, num_gauges
                allocate(gauges(i)%q_out_vars(num_eqn))
                read(UNIT, *) gauges(i)%q_out_vars

                ! Count number of vars to be output
                gauges(i)%num_out_vars = 0
                do n = 1, size(gauges(i)%q_out_vars, 1)
                    if (gauges(i)%q_out_vars(n)) then
                        gauges(i)%num_out_vars = gauges(i)%num_out_vars + 1
                    end if
                end do
            end do

            ! Read in aux fields (num_aux > 0 for geoclaw)
            read(UNIT, *)
            read(UNIT, *)
            do i = 1, num_gauges
                allocate(gauges(i)%aux_out_vars(num_aux))
                read(UNIT, *) gauges(i)%aux_out_vars

                ! Count number of vars to be output
                do n = 1, size(gauges(i)%aux_out_vars, 1)
                    if (gauges(i)%aux_out_vars(n)) then
                        gauges(i)%num_out_vars = gauges(i)%num_out_vars + 1
                    end if
                end do
            end do

            ! Count eta as one of the out vars
            gauges(:)%num_out_vars = gauges(:)%num_out_vars + num_layers

            close(UNIT)
            ! Done reading =====================================================

            ! Allocate data buffer - Note extra var out due to eta
            do i = 1, num_gauges
                allocate(gauges(i)%data(gauges(i)%num_out_vars + 2, MAX_BUFFER))
            end do

            ! Create gauge output files
            do i = 1, num_gauges
                gauges(i)%file_name = 'gaugexxxxx.txt'
                num = gauges(i)%gauge_num
                do pos = 10, 6, -1
                    digit = mod(num,10)
                    gauges(i)%file_name(pos:pos) = char(ichar('0') + digit)
                    num = num / 10
                end do

                ! Handle restart
                if (restart) then
                    open(unit=OUTGAUGEUNIT, file=gauges(i)%file_name,       &
                         status='old', position='append', form='formatted')
                else
                    open(unit=OUTGAUGEUNIT, file=gauges(i)%file_name,       &
                         status='unknown', position='append', form='formatted')
                    rewind OUTGAUGEUNIT

                    ! Write header
                    header_1 = "('# gauge_id= ',i5,' " //                   &
                               "location=( ',1e15.7,' ',1e15.7,' ) " //     &
                               "num_var= ',i2)"
                    write(OUTGAUGEUNIT, header_1) gauges(i)%gauge_num,      &
                                                  gauges(i)%x,              &
                                                  gauges(i)%y,              &
                                                  gauges(i)%num_out_vars

                    ! Construct column labels
                    index = 0
                    q_column = "["
                    do n=1, size(gauges(i)%q_out_vars, 1)
                        if (gauges(i)%q_out_vars(n)) then
                            write(q_column(3 * index + 2:4 + 3 * index), "(i3)") n
                            index = index + 1
                        end if
                    end do
                    q_column(3 * index + 2:4 + 3 * index) = "],"

                    ! Add eta
                    index = 0
                    eta_column = "["
                    do n=1, num_layers
                        write(eta_column(3 * index + 2:4 + 3 * index), "(i3)") n
                        index = index + 1
                    end do
                    eta_column(3 * index + 2:4 + 3 * index) = "],"

                    aux_column = "["
                    index = 0
                    do n=1, size(gauges(i)%aux_out_vars, 1)
                        if (gauges(i)%aux_out_vars(n)) then
                            write(aux_column(3 * index + 2:4 + 3 * index), "(i3)") n
                            index = index + 1
                        end if  
                    end do
                    aux_column(3 * index + 2:4 + 3 * index) = "]"

                    write(OUTGAUGEUNIT, *) "# level, time, q", trim(q_column), &
                                                      " eta", trim(eta_column),&
                                                      " aux", trim(aux_column)
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
        mbestsrc = 0

        !! reorder loop for better performance with O(10^5) grids
        !! for each gauge find best source grid for its data
        do 40 i = 1, num_gauges

           do 30 lev = lfine, 1, -1
            mptr = lstart(lev)
 20              if ((gauges(i)%x >= rnode(cornxlo,mptr)) .and. &
                        (gauges(i)%x <= rnode(cornxhi,mptr)) .and. &  
                        (gauges(i)%y >= rnode(cornylo,mptr)) .and. &
                        (gauges(i)%y <= rnode(cornyhi,mptr)) ) then
                        mbestsrc(i) = mptr
                      !! best source found for this gauge, go to next one
                      !! we know its the best because we started at finest level
                      go to 40  ! on to next gauge 
                 else 
                mptr = node(levelptr, mptr)
                    if (mptr .ne. 0) then
                       go to 20  ! try another grid
                    else
                       go to 30  ! try next coarser level grids
            end if
                 end if
 30        continue 

          if (mbestsrc(i) .eq. 0) &
              print *, "ERROR in setting grid src for gauge data", i
 40     continue 

!!!  NO MORE qsort and mbestg arrays. 
!!! Each grid now loops over mbestsrc array to see which gauges it owns.

    end subroutine setbestsrc


    ! ==========================================================================

    subroutine update_gauges(q, aux, xlow, ylow, num_eqn, mitot, mjtot, num_aux, &
                             mptr)
!
!     This routine is called each time step for each grid patch, to output
!     gauge values for all gauges for which this patch is the best one to 
!     use (i.e. at the finest refinement level).  

!     It is called after ghost cells have been filled from adjacent grids
!     at the same level, so bilinear interpolation can be used to 
!     to compute values at any gauge location that is covered by this grid.  

!     The grid patch is designated by mptr.
!     We only want to set gauges i for which mbestsrc(i) == mptr.
!     The array mbestsrc is reset after each regridding to indicate which
!     grid patch is best to use for each gauge.

!     This is a refactoring of dumpgauge.f from Clawpack 5.2 
!     Loops over only the gauges to be handled by this grid, as specified
!     by indices from mbestg1(mptr) to mbestg2(mptr)
!     NO MORE mbestg1 and 2.  Loop through all gauges. No sorting too.

        use amr_module, only: nestlevel, nghost, timemult, rnode, node, maxvar
        use amr_module, only: hxposs, hyposs

        use geoclaw_module, only: rho

        use multilayer_module, only: num_layers, dry_tolerance, eta_init

        implicit none

        ! Input
        integer, intent(in) ::  num_eqn, mitot, mjtot, num_aux, mptr
        real(kind=8), intent(in) :: q(num_eqn, mitot, mjtot)
        real(kind=8), intent(in) :: aux(num_aux, mitot, mjtot)
        real(kind=8), intent(in) :: xlow, ylow

        ! Locals
        real(kind=8) :: var(maxvar * 2)
        real(kind=8) :: xcent, ycent, xoff, yoff, tgrid, hx, hy
        integer :: level, i1, i2, icell, jcell, iindex, jindex
        integer :: i, j, n, m, k, var_index, layer_index, eta_index
        real(kind=8) :: h(4), mod_dry_tolerance, topo
        real(kind=8) :: h_interp(num_layers), eta(num_layers)
        logical :: dry_state

        ! No gauges to record, exit
        if (num_gauges == 0) then
            return
        endif

        ! Grid info
        tgrid = rnode(timemult, mptr)
        level = node(nestlevel, mptr)
        hx = hxposs(level)
        hy = hyposs(level)

        ! Main Gauge Loop ======================================================
        do i = 1, num_gauges
            if (mptr .ne. mbestsrc(i)) cycle 
            if (tgrid < gauges(i)%t_start .or. tgrid > gauges(i)%t_end) then
               cycle
            endif
            ! Minimum increment
            ! TODO Maybe always allow last time output recording?
            if (tgrid - gauges(i)%last_time < gauges(i)%min_time_increment) then
                cycle
            end if

            ! Compute indexing and bilinear interpolant weights
            ! Note: changed 0.5 to  0.5d0 etc.
                iindex =  int(.5d0 + (gauges(i)%x - xlow) / hx)
                jindex =  int(.5d0 + (gauges(i)%y - ylow) / hy)
            if ((iindex < nghost .or. iindex > mitot-nghost) .or. &
                (jindex < nghost .or. jindex > mjtot-nghost)) then
                    print *, "ERROR in output of Gauge Data "
            end if
            xcent  = xlow + (iindex - 0.5d0) * hx
            ycent  = ylow + (jindex - 0.5d0) * hy
                xoff   = (gauges(i)%x - xcent) / hx
                yoff   = (gauges(i)%y - ycent) / hy

            ! Gauge interpolation seems to work, so error test is commented out.
            ! For debugging, use the code below...
            !   Note: we expect 0 <= xoff, yoff <= 1 but if gauge is exactly 
            !   at center of cell these might be off by rounding error

            !if (xoff .lt. -1.d-4 .or. xoff .gt. 1.0001d0 .or. &
            !    yoff .lt. -1.d-4 .or. yoff .gt. 1.0001d0) then
            !   write(6,*) "*** print_gauges: Interpolation problem at gauge ",&
                !               igauge(i)
            !   write(6,*) "    xoff,yoff: ", xoff,yoff
            !endif


            ! Modified below from amrclaw/src/2d/gauges_module.f90 
            ! to interpolate only where all four cells are
            ! wet, otherwise just take this cell value:

            ! Check for dry cells by comparing h to mod_dry_tolerance, which 
            ! should be smaller than drytolerance to avoid oscillations since  
            ! when h < drytolerance the velocities are zeroed out which can then 
            ! lead to increase in h again.

            ! Layer Loop =======================================================
            dry_state = .false.
            
            ! Count for number of variables written to var
            var_index = 0

            do m = 1, num_layers
                layer_index = 3 * (m - 1)
                mod_dry_tolerance = 0.1d0 * dry_tolerance(m)

                h(1) = q(layer_index + 1, iindex, jindex) / rho(m)
                h(2) = q(layer_index + 1, iindex + 1, jindex) / rho(m)
                h(3) = q(layer_index + 1, iindex, jindex + 1) / rho(m)
                h(4) = q(layer_index + 1, iindex + 1, jindex + 1) / rho(m)
            

                if ((h(1) < mod_dry_tolerance) .or.  &
                    (h(2) < mod_dry_tolerance) .or.  &
                    (h(3) < mod_dry_tolerance) .or.  &
                    (h(4) < mod_dry_tolerance)) then

                    ! One of the cells is dry, so just use value from grid cell
                    ! that contains gauge rather than interpolating
                    
                    dry_state = .true.

                    icell = int(1.d0 + (gauges(i)%x - xlow) / hx)
                    jcell = int(1.d0 + (gauges(i)%y - ylow) / hy)

                    ! Note here that we are directly assuming the number of
                    ! layer variables, handle other q variables below outside of
                    ! the layer loop
                    do n = layer_index + 1, layer_index + 3
                        if (gauges(i)%q_out_vars(n)) then
                            var_index = var_index + 1
                            var(var_index) = q(n, icell, jcell) / rho(m)
                        end if
                    end do

                    h_interp(m) = q(layer_index + 1, icell, jcell)
                else
                    ! Linear interpolation between four cells
                    do n = layer_index + 1, layer_index + 3
                        if (gauges(i)%q_out_vars(n)) then
                            var_index = var_index + 1
                            var(var_index) = &
                                    (1.d0 - xoff) * (1.d0 - yoff) * q(n,     iindex,     jindex) / rho(m)  &
                                  + (       xoff) * (1.d0 - yoff) * q(n, 1 + iindex,     jindex) / rho(m)  &
                                  + (1.d0 - xoff) * (       yoff) * q(n,     iindex, 1 + jindex) / rho(m)  &
                                  + (       xoff) * (       yoff) * q(n, 1 + iindex, 1 + jindex) / rho(m)
                        end if
                    enddo
                    
                    ! Explicitly do depth in case the depth is not computed above
                    ! Here we depart from the single-layer case as we would
                    ! otherwise need to compute exactly where the h might have 
                    ! been stored
                    h_interp(m) = &
                                (1.d0 - xoff) * (1.d0 - yoff) * q(layer_index + 1,     iindex,     jindex) / rho(m)  &
                              + (       xoff) * (1.d0 - yoff) * q(layer_index + 1, 1 + iindex,     jindex) / rho(m)  &
                              + (1.d0 - xoff) * (       yoff) * q(layer_index + 1,     iindex, 1 + jindex) / rho(m)  &
                              + (       xoff) * (       yoff) * q(layer_index + 1, 1 + iindex, 1 + jindex) / rho(m)
                endif
            end do  ! End of layer loop ========================================

            ! Handle the rest of the fields
            if (dry_state) then

                ! Compute rest of q fields
                do n = 3 * num_layers + 1, size(gauges(i)%q_out_vars, 1)
                    if (gauges(i)%q_out_vars(n)) then
                        var(var_index) = 0.0
                    end if
                end do

                ! Skip over eta spots
                eta_index = var_index
                var_index = var_index + num_layers

                ! Compute aux fields
                do n=1, size(gauges(i)%aux_out_vars, 1)
                    if (gauges(i)%aux_out_vars(n)) then
                        var_index = var_index + 1
                        var(var_index) = aux(n, icell , jcell)
                    end if
                end do
                
                ! Topography
                topo = aux(1, icell, jcell)

            else

                ! Compute rest of q fields
                do n = 3 * num_layers + 1, size(gauges(i)%q_out_vars, 1)
                    if (gauges(i)%q_out_vars(n)) then
                        var_index = var_index + 1
                        var(var_index) = &
                                (1.d0 - xoff) * (1.d0 - yoff) * q(n,     iindex,     jindex)  &
                              + (       xoff) * (1.d0 - yoff) * q(n, 1 + iindex,     jindex)  &
                              + (1.d0 - xoff) * (       yoff) * q(n,     iindex, 1 + jindex)  &
                              + (       xoff) * (       yoff) * q(n, 1 + iindex, 1 + jindex)
                    end if
                end do

                ! Skip eta spots
                eta_index = var_index
                var_index = var_index + num_layers

                do n = 1, size(gauges(i)%aux_out_vars, 1)
                    if (gauges(i)%aux_out_vars(n)) then
                        var_index = var_index + 1
                        var(var_index) = &
                                (1.d0 - xoff) * (1.d0 - yoff) * aux(n,     iindex,     jindex)  &
                              + (       xoff) * (1.d0 - yoff) * aux(n, 1 + iindex,     jindex)  &
                              + (1.d0 - xoff) * (       yoff) * aux(n,     iindex, 1 + jindex)  &
                              + (       xoff) * (       yoff) * aux(n, 1 + iindex, 1 + jindex)
                    end if
                enddo
                    
                topo = (1.d0 - xoff) * (1.d0 - yoff) * aux(1,     iindex,     jindex)  &
                     + (       xoff) * (1.d0 - yoff) * aux(1, 1 + iindex,     jindex)  &
                     + (1.d0 - xoff) * (       yoff) * aux(1,     iindex, 1 + jindex)  &
                     + (       xoff) * (       yoff) * aux(1, 1 + iindex, 1 + jindex)

            end if

            ! Check to make sure we grabbed all the values
            if (gauges(i)%num_out_vars /= var_index) then
                print *, gauges(i)%num_out_vars, var_index
                print *, gauges(i)%q_out_vars
                print *, gauges(i)%aux_out_vars
                stop "Somehow we did not grab all the values we wanted..."
            end if

            ! Zero out tiny values to prevent later problems reading data,
            ! as done in valout.f
            do j = 1, gauges(i)%num_out_vars
               if (abs(var(j)) < 1d-90) var(j) = 0.d0
            end do
    
            ! Extract surfaces - Note that we do this after the zeroing to
            ! preserve etas that are non-zero
            eta(num_layers) = h_interp(num_layers) + topo
            if (abs(eta(num_layers) - eta_init(num_layers)) < 1d-90) then
                eta(num_layers) = eta_init(num_layers)
            end if
            var(eta_index + num_layers) = eta(num_layers)
            do k = num_layers - 1, 1, -1
                eta(k) = h_interp(k) + eta(k + 1)
                if (abs(eta(k) - eta_init(k)) < 1d-90) eta(k) = eta_init(k)
                var(eta_index + k) = eta(k)
            end do

            ! save info for this time 
            n = gauges(i)%buffer_index
     
            gauges(i)%level(n) = level
            gauges(i)%data(1,n) = tgrid
            do j = 1, gauges(i)%num_out_vars
                gauges(i)%data(1 + j, n) = var(j)
            end do
            
            gauges(i)%buffer_index = n + 1
            if (gauges(i)%buffer_index > MAX_BUFFER) then
                call print_gauges_and_reset_nextLoc(i)
            endif

            gauges(i)%last_time = tgrid

        end do ! End of gauge loop =============================================
 
    end subroutine update_gauges
!
! -------------------------------------------------------------------------
! Write out gauge data for the gauge specified
!
    subroutine print_gauges_and_reset_nextLoc(gauge_num)

        implicit none

        ! Input
        integer, intent(in) :: gauge_num

        ! Locals
        integer :: j, k, myunit
        integer :: omp_get_thread_num, mythread
        character(len=32) :: out_format

        ! Open unit dependent on thread number
        mythread = 0
!$      mythread = omp_get_thread_num()
        myunit = OUTGAUGEUNIT + mythread

        ! ASCII output
        if (gauges(gauge_num)%file_format == 1) then
            ! Construct output format based on number of output variables and
            ! request format
            write(out_format, "(A7, i2, A6, A1)") "(i5.2,",                    &
                                        gauges(gauge_num)%num_out_vars + 1,    &
                                        gauges(gauge_num)%display_format, ")"

            open(unit=myunit, file=gauges(gauge_num)%file_name, status='old', &
                              position='append', form='formatted')
          
            ! Loop through gauge's buffer writing out all available data.  Also
            ! reset buffer_index back to beginning of buffer since we are emptying
            ! the buffer here
            do j = 1, gauges(gauge_num)%buffer_index - 1
                write(myunit, out_format) gauges(gauge_num)%level(j),    &
                    (gauges(gauge_num)%data(k, j), k=1,                  &
                                             gauges(gauge_num)%num_out_vars + 1)
            end do
            gauges(gauge_num)%buffer_index = 1                        

            ! close file
            close(myunit)
        else
            print *, "Unhandled file format ", gauges(gauge_num)%file_format
            stop
        end if

      end subroutine print_gauges_and_reset_nextLoc
  
end module gauges_module
