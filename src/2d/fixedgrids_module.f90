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
          
        ! Construct NaN for filling empty spots
        integer(kind=16) NaN_descriptor
        real(kind=8) :: NaN

        ! Allocation pointer for fixed grid data
        real(kind=8), pointer :: temp_data(:,:,:)

        ! Other locals
        integer :: i

        data NaN_descriptor/B'01111111100000100000000000000000'/
        NaN = transfer(NaN,NaN_descriptor)

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
      ! Also, initialize to Nan to prevent non-filled values from being 
      ! misinterpreted
      do i=1,num_fixed_grids
          allocate(temp_data(num_grid_vars(i,1),mx_fg(i),my_fg(i)))
          temp_data = NaN
!           call move_alloc(temp_data,early_data_fg(i))
          early_data_fg(i)%data => temp_data
      enddo
      do i=1,num_fixed_grids
          allocate(temp_data(num_grid_vars(i,2),mx_fg(i),my_fg(i)))
          temp_data = NaN
!           call move_alloc(temp_data,late_data_fg(i))
          late_data_fg(i)%data => temp_data
      enddo
      do i=1,num_fixed_grids
          allocate(temp_data(num_grid_vars(i,2),mx_fg(i),my_fg(i)))
          temp_data = NaN
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
     
        ! NaN for filling empty spots
        integer(kind=16) NaN_descriptor
        real(kind=8) :: NaN

        data NaN_descriptor/B'01111111100000100000000000000000'/
        NaN = transfer(NaN,NaN_descriptor)
     
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

end module
