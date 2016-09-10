! ==============================================================================
!  Variable friction module
!
!    This module contains storage and routines for dealing with variable 
!    friction fields in GeoClaw.
!
!    If variable_friction > 0 then the data located here will override the basic
!    friction coefficient data in GeoClaw and will be ignored.
!
!    Precidence of specification:
!      file based
!      region based
!      depth based
!      constant value from geoclaw_module
!    
! ==============================================================================

module friction_module

    implicit none
    save

    logical, private :: module_setup = .false.

    ! Parameters
    integer, public, parameter :: friction_index = 4

    ! Whether to use variable friction and what type to specify
    logical, public :: variable_friction

    ! Support for region based specification
    type friction_region_type
        ! Bounds of region
        real(kind=8) :: lower(2), upper(2)

        ! Coefficients in the region based on depths
        real(kind=8), pointer :: depths(:)
        real(kind=8), pointer :: manning_coefficients(:)

    end type friction_region_type
    integer, private :: num_friction_regions = 0
    type(friction_region_type), pointer, private :: friction_regions(:)

    ! Support for file based specification
    type friction_file_type
        ! Domain specification
        integer :: num_cells(2)
        real(kind=8) :: lower(2), upper(2), dx(2)

        ! Field information
        real(kind=8) :: no_data_value

        ! Array containing coefficients
        real(kind=8), pointer :: data(:,:)
    end type friction_file_type

    integer, private :: num_friction_files = 0
    type(friction_file_type), pointer, private :: friction_files(:)

contains

    subroutine setup_variable_friction(path)

        use utility_module, only: get_value_count

        implicit none

        ! Input
        character(len=*), optional, intent(in) :: path

        ! Locals
        integer, parameter :: unit = 54
        integer :: m, file_type
        character(len=128) :: line

        if (.not.module_setup) then

            ! Open data file
            if (present(path)) then
                call opendatafile(unit,path)
            else
                call opendatafile(unit,'friction.data')
            endif

            ! Basic switch to turn on variable friction
            read(unit,*) variable_friction
            read(unit,'(a)')

            if (variable_friction) then
                ! Region based friction
                read(unit,'(i2)') num_friction_regions
                allocate(friction_regions(num_friction_regions))
                read(unit,*)
                do m=1,num_friction_regions
                    read(unit,*) friction_regions(m)%lower
                    read(unit,*) friction_regions(m)%upper
                    read(unit,'(a)') line
                    allocate(friction_regions(m)%depths(get_value_count(line)))
                    read(line,*) friction_regions(m)%depths
                    read(unit,'(a)') line
                    allocate(friction_regions(m)%manning_coefficients(get_value_count(line)))
                    read(line,*) friction_regions(m)%manning_coefficients
                    read(unit,*)
                enddo

                ! File based friction
                read(unit,"(i2)") num_friction_files
                allocate(friction_files(num_friction_files))
                do m=1,num_friction_files
                    read(unit,"(a,i1)") line, file_type

                    ! Open and read in friction file
                    print *,"*** WARNING *** File based friction specification unimplemented."
                    friction_files = read_friction_file(line, file_type)
                enddo
            endif

            close(unit)

            module_setup = .true.
        end if

    end subroutine setup_variable_friction

    ! ==========================================================================
    !  set_friction_field - 
    ! ==========================================================================
    subroutine set_friction_field(mx, my, num_ghost, num_aux, xlower, ylower, &
                                  dx, dy, aux)

        use geoclaw_module, only: sea_level

        implicit none

        ! Input
        integer, intent(in) :: mx, my, num_ghost, num_aux
        real(kind=8), intent(in) :: xlower, ylower, dx, dy
        real(kind=8), intent(in out) :: aux(num_aux,                           &
                                            1-num_ghost:mx+num_ghost,&
                                            1-num_ghost:my+num_ghost)

        ! Locals
        integer :: m,i,j,k
        real(kind=8) :: x, y

        if (variable_friction) then
            ! Set region based coefficients
            do m=1, num_friction_regions
                do i=1 - num_ghost, mx + num_ghost
                    do j=1 - num_ghost, my + num_ghost                        
                        x = xlower + (i-0.5d0) * dx
                        y = ylower + (j-0.5d0) * dy
                        if (friction_regions(m)%lower(1) < x .and.   &
                            friction_regions(m)%lower(2) < y .and.   &
                            friction_regions(m)%upper(1) >= x .and.  &
                            friction_regions(m)%upper(2) >= y) then

                            do k=1,size(friction_regions(m)%depths) - 1
                                if (friction_regions(m)%depths(k+1)            &
                                                <= aux(1,i,j) - sea_level.and. &
                                    friction_regions(m)%depths(k)              &
                                                 > aux(1,i,j) - sea_level) then

                                    aux(friction_index,i,j) = &
                                     friction_regions(m)%manning_coefficients(k)
                                endif
                            enddo
                        endif
                    enddo
                enddo
            enddo

            if (num_friction_files > 0) then
                stop "*** ERROR *** File based friction specification not implemented."
            endif
        end if

    end subroutine set_friction_field


    ! ==========================================================================
    !  read_friction_file - Reads an input file containing info on friction
    !    coefficients.
    ! ==========================================================================
    type(friction_file_type) function read_friction_file(path, file_type) result(file)
        
        implicit none

        ! Path to file to be read in
        character(len=*), intent(in) :: path
        integer, intent(in) :: file_type

        ! Locals
        integer, parameter :: unit = 24
!         real(kind=8), parameter :: missing_value = huge(1.d0)
        integer :: ios, missing
        
        ! Open file for reading
        open(unit=unit, file=path, iostat=ios, status="unknown",   &
                action="read", form='formatted')
        if ( ios /= 0 ) then 
            print *,"*** Error *** Could not open friction file ",path
            stop
        endif

        select case(abs(file_type))
            ! ASCII file with x,y,z values on each line.
            ! (progressing from upper left corner across rows, then down)
            ! Assumes a uniform rectangular grid of data values.
            case(1)
                stop "Unimplemented file type for friction files."
            
            ! ================================================================
            ! ASCII file with header followed by z data
            ! (progressing from upper left corner across rows, then down)
            ! one value per line if filetype=2 or
            ! mx values per line if filetype=3
            ! ================================================================
            case(2:3)
                ! Read header
                read(unit,*) file%num_cells(1)
                read(unit,*) file%num_cells(2)
                read(unit,*) file%lower(1)
                read(unit,*) file%lower(2)
                read(unit,*) file%dx(1)
                read(unit,*) file%no_data_value

                ! Calculate upper corner
                file%dx(2) = file%dx(1)
                file%upper(1) = file%lower(1) + file%dx(1) * (file%num_cells(1) - 1)
                file%upper(2) = file%lower(2) + file%dx(2) * (file%num_cells(2) - 1)

                ! Allocate data array
                allocate(file%data(file%num_cells(1),file%num_cells(2)))

                ! Read in data
                missing = 0
                select case(abs(file_type))
                    case(2)
                        stop "Unhandled file type"
!                         do i=1,product(num_cells)
!                             read(iunit,*) field%data(i,)
!                             if (auxinit(i) == no_data_value) then
!                                 missing = missing + 1
!                                 field%data(i) = missing_value
!                             endif
!                         enddo
                    case(3)
                        stop "Unhandled file type"
!                         do j=1,my
!                             read(iunit,*) (auxinit((j-1)*mx + i),i=1,mx)
!                             do i=1,mx
!                                 if (auxinit((j-1)*mx + i) == no_data_value) then
!                                     missing = missing + 1
!                                     auxinit((j-1)*mx + i) = missing_value
!                                 endif
!                             enddo
!                         enddo
                end select
            
            case default
                stop "Unimplemented file type for friction files."
        end select
        
    end function read_friction_file

end module friction_module
