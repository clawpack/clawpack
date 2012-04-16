module fixedgrids_module

    implicit none
    save

    type fixed_grid_data_type
        real(kind=8), pointer :: data(:,:,:)
    end type fixed_grid_data_type

    ! Number of fixed grids
    integer :: num_fixed_grids
    
    ! Primary data storage, these are each an array of pointers
    integer, allocatable :: num_grid_vars(:,:)
    type(fixed_grid_data_type), allocatable :: early_data_fg(:)
    type(fixed_grid_data_type), allocatable :: often_data_fg(:)
    type(fixed_grid_data_type), allocatable :: late_data_fg(:)
    
    ! Parameters for output of each fixed grid
    integer, allocatable :: arrival_times_output(:),surface_max_output(:)
    integer, allocatable :: num_output_fg(:)
    real(kind=8), allocatable :: t_start_fg(:),t_end_fg(:)

    ! Output tracking
    integer, allocatable :: t_last_output_index_fg(:)
    real(kind=8), allocatable :: t_last_output_fg(:)
    real(kind=8) :: max_fg_time
     
    ! Geometry and time step size
    integer, allocatable :: mx_fg(:),my_fg(:)
    real(kind=8), allocatable :: x_low_fg(:),x_hi_fg(:)
    real(kind=8), allocatable :: y_low_fg(:),y_hi_fg(:)
    real(kind=8), allocatable :: dt_fg(:),dx_fg(:),dy_fg(:)
       
contains

    ! Read in fixed grid settings and setup data structures
    subroutine set_fixed_grids(fname)

        use amr_module, only: parmunit

        implicit none
        
        ! Subroutine arguments
        character(len=25), optional :: fname
        
        ! File opening
        integer, parameter :: unit = 7
        character(len=*), parameter :: fg_line_format = "(2d16.8,1i2,4d16.8,2i4,2d16.8)"
        
        ! Allocation pointer for fixed grid data
        real(kind=8), pointer :: temp_data(:,:,:)

        ! Other locals
        integer :: i

        write(parmunit,*) ' '
        write(parmunit,*) '--------------------------------------------'
        write(parmunit,*) 'SETFIXEDGRIDS:'
        write(parmunit,*) '-----------'

        ! Open file
        if (present(fname)) then
            call opendatafile(unit,fname)
        else
            call opendatafile(unit,'setfixedgrids.data')
        endif

        ! Read in data
        read(7,"(i2)") num_fixed_grids
        write(parmunit,*) '  num_fixed_grids = ',num_fixed_grids
        if (num_fixed_grids == 0) then
            write(parmunit,*) "  No fixed grids specified for output"
            return
        endif
        
        ! Allocate all fixed grid data and info arrays
        allocate(t_start_fg(num_fixed_grids),t_end_fg(num_fixed_grids))
        allocate(num_output_fg(num_fixed_grids))
        allocate(x_low_fg(num_fixed_grids),x_hi_fg(num_fixed_grids))
        allocate(y_low_fg(num_fixed_grids),y_hi_fg(num_fixed_grids))
        allocate(mx_fg(num_fixed_grids),my_fg(num_fixed_grids))
        allocate(dt_fg(num_fixed_grids))
        allocate(dx_fg(num_fixed_grids),dy_fg(num_fixed_grids))
        allocate(arrival_times_output(num_fixed_grids))
        allocate(surface_max_output(num_fixed_grids))
        allocate(num_grid_vars(num_fixed_grids,2))
        allocate(t_last_output_index_fg(num_fixed_grids))
        allocate(t_last_output_fg(num_fixed_grids))
        
        ! These are the data arrays themselves (only pointers)
        allocate(early_data_fg(num_fixed_grids))
        allocate(late_data_fg(num_fixed_grids))
        allocate(often_data_fg(num_fixed_grids))
        
        ! Read in parameters for each fixed grid
        do i=1,num_fixed_grids
            read(unit,*) t_start_fg(i),t_end_fg(i),num_output_fg(i), &
                                   x_low_fg(i),x_hi_fg(i),y_low_fg(i),y_hi_fg(i), &
                                   mx_fg(i),my_fg(i), &
                                   arrival_times_output(i),surface_max_output(i)
            write(parmunit,*) t_start_fg(i),t_end_fg(i),num_output_fg(i), &
                                   x_low_fg(i),x_hi_fg(i),y_low_fg(i),y_hi_fg(i), &
                                   mx_fg(i),my_fg(i), &
                                   arrival_times_output(i),surface_max_output(i)
        enddo
        close(unit)
       
        ! Set some parameters for each grid
        do i=1,num_fixed_grids
            ! Set time step length between outputs
            if (t_end_fg(i) <= t_start_fg(i)) then
                if (num_output_fg(i) > 1) then
                    print *,'SETFIXEDGRIDS: ERROR for fixed grid', i
                    print *,'tstartfg=tendfg yet noutfg>1'
                    print *,'set tendfg > tstartfg or set noutfg = 1'
                    stop
                else
                    dt_fg(i) = 0.d0
                endif
            else
                if (num_output_fg(i) < 2) then
                    print *,'SETFIXEDGRIDS: ERROR for fixed grid', i
                    print *,'tendfg>tstartfg, yet noutfg=1'
                    print *,'set noutfg > 2'
                    stop
                else
                    dt_fg(i) = (t_end_fg(i) - t_start_fg(i)) / real(num_output_fg(i)-1,kind=8)
                endif
            endif
            
            ! Counters for keeping track of output times
            t_last_output_fg(i) = t_start_fg(i) - dt_fg(i)
            t_last_output_index_fg(i) = 0

            ! Set spatial intervals dx and dy on each grid
            if (mx_fg(i) > 1) then
                dx_fg(i) = (x_hi_fg(i) - x_low_fg(i)) / real(mx_fg(i) - 1,kind=8)
            else if (mx_fg(i) == 1) then
                dx_fg(i) = 0.d0
            else
                print *,'SETFIXEDGRIDS: ERROR for fixed grid', i
                print *,'x grid points mxfg<=0, set mxfg>= 1'
            endif
            if (my_fg(i) > 1) then
                dy_fg(i) = (y_hi_fg(i) - y_low_fg(i)) / real(my_fg(i) - 1,kind=8)
            else if (my_fg(i) == 1) then
                dy_fg(i) = 0.d0
            else
                print *,'SETFIXEDGRIDS: ERROR for fixed grid', i
                print *,'y grid points myfg<=0, set myfg>= 1'
            endif

        enddo

      ! Set the number of variables stored for each grid
      ! this should be (the number of variables you want to write out + 1)
      do i=1, num_fixed_grids 
         num_grid_vars(i,1)  = 6  
         num_grid_vars(i,2) = 3 * surface_max_output(i) + arrival_times_output(i)
      enddo

      ! Allocate data arrays, here we use allocate data to the temprorary
      ! pointer and transfer the data to the array of pointers stored in
      ! each of the individual arrays of pointers.  Since we only do this once
      ! this seems like it should not be intolerable in terms of performance
      !
      ! This code used to fill grid data with NaNs to prevent non-filled values
      ! from being misinterpreted, this was not portable and so the arrays are
      ! filled with the intrinsic `huge` instead (largest representable number)
      do i=1,num_fixed_grids
          allocate(temp_data(num_grid_vars(i,1),mx_fg(i),my_fg(i)))
          temp_data = huge(1.d0)
!           call move_alloc(temp_data,early_data_fg(i))
          early_data_fg(i)%data => temp_data
      enddo
      do i=1,num_fixed_grids
          allocate(temp_data(num_grid_vars(i,2),mx_fg(i),my_fg(i)))
          temp_data = huge(1.d0)
!           call move_alloc(temp_data,late_data_fg(i))
          late_data_fg(i)%data => temp_data
      enddo
      do i=1,num_fixed_grids
          allocate(temp_data(num_grid_vars(i,2),mx_fg(i),my_fg(i)))
          temp_data = huge(1.d0)
!           call move_alloc(temp_data,often_data_fg(i))
          often_data_fg(i)%data => temp_data
      enddo
      
      ! Set convenience maximum output for fixed grid time
      max_fg_time = 1d-16

    end subroutine set_fixed_grids

    ! This routine interpolates q and aux on a computational grid
    ! to a fgrid not necessarily aligned with the computational grid
    ! using bilinear interpolation defined on computation grid
    subroutine fgridinterp(fgrid,x_low_fg,y_low_fg, &
                           x_hi_fg,y_hi_fg,dx_fg,dy_fg,mx_fg,my_fg,t,num_vars_fg, &
                           q,meqn, &
                           mx_c,my_c,mbc,dx_c,dy_c,num_vars,x_low_c,y_low_c, &
                           maux,aux, &
                           out_arrival_times,out_surface_max,max_check)
     
        use geoclaw_module, only: drytolerance
     
        implicit none
     
        ! Subroutine arguments
        integer, intent(in) :: mx_fg,my_fg,num_vars_fg
        real(kind=8), intent(in) :: x_low_fg,x_hi_fg,y_low_fg,y_hi_fg,dx_fg,dy_fg
        real(kind=8), intent(inout) :: fgrid(num_vars_fg,1:mx_fg,1:my_fg)
     
        integer, intent(in) :: mbc,mx_c,my_c,num_vars,meqn,maux
        real(kind=8), intent(in) :: x_low_c,y_low_c,dx_c,dy_c
        real(kind=8), intent(in) :: q(meqn,1-mbc:mx_c+mbc,1-mbc:my_c+mbc)
        real(kind=8), intent(in) :: aux(maux,1-mbc:mx_c+mbc,1-mbc:my_c+mbc)
     
        integer, intent(in) :: out_arrival_times,out_surface_max,max_check
        real(kind=8), intent(in) :: t
     
        ! Locals
        integer :: bathy_index,arrival_index,i,j
        integer :: eta_index,eta_index_min,eta_index_max,eta_index_now
        real(kind=8), parameter :: arrival_tol = 1d-2
        real(kind=8) :: x_hi_c,y_hi_c,tol
     
        ! Calculate and extract some important data
        x_hi_c = x_low_c + dx_c * mx_c
        y_hi_c = y_low_c + dy_c * my_c
    
        tol = drytolerance
    
        ! Calculate some of the storage indices
        bathy_index = meqn + 1
        eta_index = meqn + 2
    
        if (max_check > 0) then
            arrival_index = 0
            eta_index_min = 0
            eta_index_max = 0
            eta_index_now = 0
        
            if (out_arrival_times > 0) then
                arrival_index = 1
            endif
        
            if (out_surface_max > 0) then
                eta_index_now = arrival_index + 1
                eta_index_min = arrival_index + 2
                eta_index_max = arrival_index + 3
            endif
        endif
        
        do i=1,mx_fg
            do j=1,my_fg
            
            enddo
        enddo
    end subroutine fgridinterp
    
    ! Advance (output) all fgrids at all times that have not yet been output
    ! but that have been bracketed by computational times.
    subroutine fgrid_advance(t,dt)
        
!         use amr_module
!         use regions_module
!         use qinit_module
!         use gauges_module
        
        implicit none
    
        ! Subroutine arguments
        real(kind=8), intent(in) :: t,dt
        
        ! Local storage
        real(kind=8) :: tc_start,tc_final,out_time
        integer :: n,i,out_start_index,out_end_index,out_flag
        
        ! Store computational step times
        tc_start = t
        tc_final = tc_start + dt
        
        ! Check to see if any fixed grids should be outputed
        do n=1,num_fixed_grids
            if (tc_start > t_start_fg(n) .and. t_last_output_index_fg(n) < num_output_fg(n)) then
                ! fgrid n may need to be written out
                ! find the first output number that has not been written out and
                ! find the first output number on a fixed grid that is >= tc0
                ! which will not be written out
                if (dt_fg(n) > 0.d0) then
                    out_end_index = 1 + max(0,int((tc_start - t_start_fg(n)) / dt_fg(n)))
                else
                    out_end_index = 1
                endif
                out_end_index = min(out_end_index,num_output_fg(n))
                out_start_index = t_last_output_index_fg(n)

                ! write-out fgrid times that are less than tc0, and have not been 
                ! written yet these should be the most accurate values at any given 
                ! point in the fgrid since tc0> output time
                do i=out_start_index,out_end_index
                    out_time = t_start_fg(n) + (i-1) * dt_fg(n)
                    if (out_time < tc_start) then
                        ! Write out solution for fixed grid n
                        out_flag = arrival_times_output(n) * (num_output_fg(n) - t_last_output_index_fg(n))
                        
                        call fgrid_output(early_data_fg(n)%data, &
                                          late_data_fg(n)%data,  &
                                          often_data_fg(n)%data, &
                                          x_low_fg(n),x_hi_fg(n),y_low_fg(n),y_hi_fg(n), &
                                          mx_fg(n),my_fg(n),num_grid_vars(n,:), &
                                          out_time,i,n,arrival_times_output(n),out_flag)
                        
                        t_last_output_fg(n) = out_time
                        t_last_output_index_fg(n) = t_last_output_index_fg(n) + 1
                    endif
                enddo
            endif
        enddo
        
    end subroutine fgrid_advance

    subroutine fgrid_output(fgrid1,fgrid2,fgrid3, &
                        x_low,x_hi,y_low,y_hi,mx,my,num_vars, &
                        t_out,num_out,grid_number, &
                        out_arrival,out_flag)

        use geoclaw_module

        implicit none
    
        ! Subroutine arguments
        integer, intent(in) :: num_vars(2),mx,my
        integer, intent(in) :: num_out,grid_number,out_arrival,out_flag
        real(kind=8), intent(in) :: x_low,x_hi,y_low,y_hi,t_out
        real(kind=8), intent(inout), dimension(num_vars(1),1:mx,1:my) :: fgrid1,fgrid2,fgrid3

        ! File handling
        character(len=30) :: fg_file_name
        integer, parameter :: unit = 90
        integer :: digit,nout,ng,pos,status
    
        ! Other locals
        integer :: i,j,m,deta_min,deta_max,num_columns
        real(kind=8) :: t0,tf,tau
    
        ! Format strings
        character(len=*), parameter :: header_format = "(e18.8,'    time', /," // &
                                                           "i5,'    mx', /,"   // &
                                                           "i5,'    my', /,"   // &
                                                        "e18.8,'    xlow',/"   // &
                                                        "e18.8,'    ylow',/"   // &
                                                        "e18.8,'    xhi',/,"   // &
                                                        "e18.8,'    yhi',/,"   // &
                                                           "i5,'  columns',/)"
        character(len=*), parameter :: data_format = "(8e26.16)"
        character(len=*), parameter :: arrival_header_format = &
                                                          "(i5,'    mx', /," // &
                                                           "i5,'    my', /," // &
                                                        "e18.8,'    xlow',/" // &
                                                        "e18.8,'    ylow',/" // &
                                                        "e18.8,'    xhi',/," // &
                                                        "e18.8,'    yhi',/," // &
                                                           "i5,'  columns',/)"
        character(len=*), parameter :: arrival_data_format = "(1e26.16)"

        ! Open file for writing                               
        fg_file_name = 'fort.fgnn_xxxx'
        ng = grid_number
        do pos=9,8,-1
            digit = mod(ng,10)
            fg_file_name(pos:pos) = char(ichar('0') + digit)
            ng = ng / 10
        enddo
    
        nout = num_out
        do pos=14,11,-1
            digit = mod(nout,10)
            fg_file_name(pos:pos) = char(ichar('0') + digit)
            nout = nout / 10
        enddo
    
        open(unit=unit, file=fg_file_name, iostat=status, status="new", action="write")
        if ( status /= 0 ) then
            print *, "Error opening file fixed grid output file ",fg_file_name
            stop
        endif
    
        ! Determine the number of columns of the output file
        num_columns = num_vars(1) - 1
        if (num_vars(2) > 1) then
            num_columns = num_columns + 2
        endif
        deta_min = out_arrival + 2
        deta_max = out_arrival + 3
    
        ! Write out header
        write(unit,header_format) t_out,mx,my,x_low,y_low,x_hi,y_hi,num_columns

        ! Interpolate the grid in time to hit the requested output time using the
        ! solution in fgrid1 and fgrid2 which represent the solution on the fixed
        ! grid at the two nearest computational times
        do j=1,my
            do i=1,mx
                ! Figure out time interpolant
                t0 = fgrid1(num_vars(1),i,j)
                tf = fgrid2(num_vars(1),i,j)
                tau = (t_out - t0) / (tf - t0)
            
                ! Zero out values that are too small
                do m=1,num_vars(1)-1
                    if (abs(fgrid1(m,i,j)) < 1d-90) then
                        fgrid1(m,i,j) = 0.d0
                    endif
                    if (abs(fgrid2(m,i,j)) < 1d-90) then
                        fgrid2(m,i,j) = 0.d0
                    endif
                enddo
            
                ! Write out interpolants
                if (num_columns == num_vars(1) - 1) then
                    write(unit,data_format) ((1.d0 - tau) * fgrid1(m,i,j) &
                                + tau * fgrid2(m,i,j),m=1,num_vars(1)-1)
                else
                    if (abs(fgrid3(deta_min,i,j)) < 1d-90) then
                        fgrid3(deta_min,i,j) = 0.d0
                    endif
                    if (abs(fgrid3(deta_max,i,j)) < 1d-90) then
                        fgrid3(deta_max,i,j) = 0.d0
                    endif
                
                    write(unit,data_format) ((1.d0 - tau) * fgrid1(m,i,j) &
                            + tau * fgrid2(m,i,j), m=1,num_vars(1)-1), &
                            fgrid3(deta_min,i,j),fgrid3(deta_max,i,j)
                endif
            enddo
        enddo
    
        close(unit)
        print "(' FGRIDOUT: Fixed Grid  ', i2, '  output at time =', e18.8)",grid_number,t_out
    
        ! Output arrival times
        if (out_flag == 1) then
            ! Open output file
            fg_file_name = 'fort.fgnn_arrivaltimes'
            ng = grid_number
            do pos=9,8,-1
                digit = mod(ng,10)
                fg_file_name(pos:pos) = char(ichar('0') + digit)
                ng = ng / 10
            enddo
        
            open(unit=unit, file=fg_file_name, iostat=status, status="new", action="write")
            if ( status /= 0 ) then
                print *,"Error opening file for arrival times file ",fg_file_name
                stop
            endif
        
            write(unit,arrival_header_format) mx,my,x_low,y_low,x_hi,y_hi
            do j=1,my
                do i=1,mx
                    write(unit,arrival_data_format) fgrid3(1,i,j)
                enddo
            enddo
        
            close(unit)
            print "(' FGRIDOUT: Fixed Grid  ', i2, '  arrival times output')",grid_number
        endif
    
    end subroutine fgrid_output

end module
