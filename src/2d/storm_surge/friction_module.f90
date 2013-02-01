! ==============================================================================
!  Variable friction module
!
!    This module contains storage and routines for dealing with variable 
!    friction fields in GeoClaw.
!
!    If variable_friction > 0 then the data located here will override the basic
!    friction coefficient data in GeoClaw and will be ignored.
!
!    Types of specifications:
!       0 = No variable friction, constant value used
!       1 = Depth based specification only
!       2 = File and depth based specification, where data from a file is not
!           available the depth based specification will be used.
! ==============================================================================

module friction_module

    implicit none
    save

    ! Parameters
    integer, parameter :: friction_index = 4

    ! Whether to use variable friction and what type to specify
    integer :: variable_friction

    ! Support for depth based specification
    real(kind=8), allocatable :: friction_depths(:)
    real(kind=8), allocatable :: manning_coefficients(:)

    ! Support for file based specification
    type friction_field_type
        ! Domain specification
        integer :: num_cells(2)
        real(kind=8) :: lower(2), upper(2), dx(2)

        ! Field information
        real(kind=8) :: no_data_value

        ! Array containing coefficients
        real(kind=8), pointer :: data(:,:)

    end type friction_field_type

    ! Multi-file support
    integer :: num_friction_files
    type(friction_field_type), pointer :: friction_field(:)

contains

    subroutine setup_variable_friction(path)

        use utility_module, only: get_value_count

        implicit none

        ! Input
        character(len=*), optional, intent(in) :: path

        ! Locals
        integer, parameter :: unit = 54
        character(len=128) :: line

        ! Open data file
        if (present(path)) then
            call opendatafile(unit,path)
        else
            call opendatafile(unit,'friction.data')
        endif

        ! Read in data
        read(unit,*) variable_friction
        if (variable_friction > 0) then
            print *,"You have elected to use a variable friction field, note"
            print *,"that this overrides all friction parameter settings from"
            print *,"the geoclaw data parameters."
        endif
        select case(variable_friction)
            case(0) ! No variable friction
                continue
            case(1:2) ! Specify friction by depth
                read(unit,'(a)') line
                allocate(friction_depths(get_value_count(line)))
                read(line,*) friction_depths
                read(unit,'(a)') line
                allocate(manning_coefficients(get_value_count(line)))
                read(line,*) manning_coefficients

                ! Specify friction by an input file
                if (variable_friction == 2) then
                    stop "Friction fields specified via file is not supported."
                endif

            case default
                stop "Invalid variable friction specification."
        end select
        close(unit)

    end subroutine setup_variable_friction

    ! ==========================================================================
    !  set_friction_field - 
    ! ==========================================================================
    subroutine set_friction_field(num_cells, num_ghost, num_aux, aux)

        use geoclaw_module, only: manning_coefficient, sea_level

        implicit none

        ! Input
        integer, intent(in) :: num_cells(2), num_ghost, num_aux
        real(kind=8), intent(in out) :: aux(num_aux,                           &
                                            1-num_ghost:num_cells(1)+num_ghost,&
                                            1-num_ghost:num_cells(2)+num_ghost)

        ! Locals
        integer :: m,i,j

        select case(variable_friction)
            case(0) 
                ! No variable friction
                aux(friction_index,:,:) = manning_coefficient

            case(1:2) 
                ! Specify friction by depth
                do m=1,size(friction_depths) - 1
                    forall(i=1 - num_ghost:num_cells(1) + num_ghost,           &
                           j=1 - num_ghost:num_cells(2) + num_ghost,           &
                           friction_depths(m+1) <= aux(1,i,j) - sea_level .and.&
                           friction_depths(m) > aux(1,i,j) - sea_level)

                        aux(friction_index,i,j) = manning_coefficients(m)
                    end forall
                end do

                ! Additionally read in data file info
                if (variable_friction == 2) then
                    ! Specify friction by an input file
                    stop "Friction fields specified via file is not supported."
                endif

            case default
                stop "Invalid variable friction specification."
        end select

    end subroutine set_friction_field


    ! ==========================================================================
    !  read_friction_file - Reads an input file containing info on friction
    !    coefficients.
    ! ==========================================================================
    type(friction_field_type) function read_friction_file(path, file_type) result(field)
        
        implicit none

        ! Path to file to be read in
        character(len=*), intent(in) :: path
        integer, intent(in) :: file_type

        ! Locals
        integer, parameter :: unit = 24
        real(kind=8), parameter :: missing_value = huge(1.d0)
        integer :: ios, missing, i
        
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
                read(unit,*) field%num_cells(1)
                read(unit,*) field%num_cells(2)
                read(unit,*) field%lower(1)
                read(unit,*) field%lower(2)
                read(unit,*) field%dx(1)
                read(unit,*) field%no_data_value

                ! Calculate upper corner
                field%dx(2) = field%dx(1)
                field%upper(1) = field%lower(1) + field%dx(1) * (field%num_cells(1) - 1)
                field%upper(2) = field%lower(2) + field%dx(2) * (field%num_cells(2) - 1)

                ! Allocate data array
                allocate(field%data(field%num_cells(1),field%num_cells(2)))

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
