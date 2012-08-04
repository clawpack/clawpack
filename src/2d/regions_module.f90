module regions_module

    implicit none
    save
      
    ! Number of regions and level extents
    integer :: num_regions
    integer, allocatable :: min_level_region(:)
    integer, allocatable :: max_level_region(:)
      
    ! Region geometry
    real(kind=8), allocatable :: x_low_region(:)
    real(kind=8), allocatable :: y_low_region(:)
    real(kind=8), allocatable :: x_hi_region(:)
    real(kind=8), allocatable :: y_hi_region(:)
    real(kind=8), allocatable :: t_low_region(:)
    real(kind=8), allocatable :: t_hi_region(:)
      
contains

    subroutine set_regions(fname)

        use amr_module
      
        implicit none
      
        ! Function Arguments
        character(len=*), optional, intent(in) :: fname
      
        ! Locals
        integer, parameter :: unit = 7
        integer :: i

        write(parmunit,*) ' '
        write(parmunit,*) '--------------------------------------------'
        write(parmunit,*) 'SETREGIONS:'
        write(parmunit,*) '-----------'

        if (present(fname)) then
            call opendatafile(unit,fname)
        else
            call opendatafile(unit,'setregions.data')
        endif

        read(unit,"(i2)") num_regions

        if (num_regions == 0) then
            write(parmunit,*) '  No regions specified for refinement'
            return
        endif
    
        ! Allocate data arrays
        allocate(min_level_region(num_regions),max_level_region(num_regions))
        allocate(x_low_region(num_regions),x_hi_region(num_regions))
        allocate(y_low_region(num_regions),y_hi_region(num_regions))
        allocate(t_low_region(num_regions),t_hi_region(num_regions))

        ! Read in data
        do i=1,num_regions
            read(unit,*) min_level_region(i),max_level_region(i), &
                                   t_low_region(i),t_hi_region(i), &
                                   x_low_region(i),x_hi_region(i), &
                                   y_low_region(i),y_hi_region(i)
        enddo
        close(unit)

        ! Output read data into param file
        write(parmunit,*) '  mregions = ',num_regions
        write(parmunit,*) '  minlevel, maxleve, tlow, thi, xlow, xhi, ylow, yhigh values:'

        do i=1,num_regions
            write(parmunit,*) min_level_region(i),max_level_region(i), &
                                        t_low_region(i),t_hi_region(i), &
                                        x_low_region(i),x_hi_region(i), &
                                        y_low_region(i),y_hi_region(i)
        enddo
    end subroutine set_regions

end module
