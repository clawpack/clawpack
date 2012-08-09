! ==============================================================================
!  holland_storm_module 
!
! Module contains routines for constructing a wind and pressure field based on
! the Holland hurricane module.  
! 
! Many of these routines have been adapted from PADCIRC version 45.12 03/17/2006
! ==============================================================================
module holland_storm_module

    implicit none
    save

    ! Storm parameters
    type holland_storm_type
        ! Storm time tracking
        integer :: storm_index
        real(kind=8) :: t_next

        ! List of locations of storm
        real(kind=8), allocatable :: track(:,:)

        ! Storm parameters
        real(kind=8), allocatable :: max_wind_radius(:)

    end type holland_storm_type

contains

    ! Setup routine for the holland model
    type(holland_storm_type) pure function storm_setup(storm_data_path) &
        result(storm)

        implicit none

        ! Path to data file
        character(len=*), optional :: storm_data_path

        ! Local storage
        integer, parameter :: data_file = 701
        integer :: lines,io_status

        ! File format string
        character(len=*), parameter :: file_format = "(a)"

        ! Open data file
        if (present(storm_data_path)) then
            open(unit=data_file,file=storm_data_path,status='old',action='read',iostat=io_status)
        else
            open(unit=data_file,file="./hurricane.data",status='old',action='read',iostat=io_status)
        endif
        if (io_status /= 0) stop "Error opening hurricane data file."

        ! Count number of data lines
        lines = 0
        do
            read(data_file,'(a170)',iostat=io_status)
            if (io_status /= 0) exit
            lines = lines + 1
        end do
        rewind(data_file)

        ! Allocate storm parameter file variables
        allocate()

        ! Now re-read the file's contents
        do i=1,lines
            read(data_file,fmt=file_format)
        enddo

        ! Set data index to first valid value for time
        storm%storm_index = 1
        storm%t_next = storm%track(3,storm_index)

    end function storm_setup

    ! ==========================================================================
    !  real(kind=8) pure date_to_seconds(year,months,days,hours,minutes,seconds)
    !    Convert time from year, month, day, hour, min, sec to seconds since the
    !    beginning of the year.
    ! ==========================================================================
    real(kind=8) pure function date_to_seconds(year,months,days,hours,minutes, &
                                               seconds) result(time)
      
        implicit none

        ! Input
        integer, intent(in) :: year, months, days, hours, minutes
        real(kind=8), intent(in) :: seconds

        ! Local storage
        integer :: total_days

        ! Count number of days
        total_days = days

        ! Add days for months that have already passed
        if (months > 1) total_days = total_days + 31
        if (months > 2) then
            if (int(year / 4) * 4 == year) then
                total_days = total_days + 29
            else
                total_days = total_days + 28
            endif
        endif
        if (months > 3)  total_days = total_days + 30
        if (months > 4)  total_days = total_days + 31
        if (months > 5)  total_days = total_days + 30
        if (months > 6)  total_days = total_days + 31
        if (months > 7)  total_days = total_days + 30
        if (months > 8)  total_days = total_days + 31
        if (months > 9)  total_days = total_days + 30
        if (months > 10) total_days = total_days + 31
        if (months > 11) total_days = total_days + 30

        ! Convert everything to seconds since the beginning of the year
        time = real((total_days - 1) * 86400 + hours * 3600 + minutes * 60,kind=8)
        time = time + seconds

    end function date_to_seconds

    ! ==========================================================================
    !  holland_update_storm(t)
    !    Update storm information to time t
    ! ==========================================================================
    subroutine update_storm(t,storm)

        implicit none
        
        ! Input
        real(kind=8), intent(in) :: t
        type(holland_storm_type), intent(in out) :: storm

        ! Update storm variable
        if (t >= storm%t_next) then
            stop "Should update time and index of storm"
        endif

    end subroutine holland_update_storm


    ! ==========================================================================
    !  storm_wind_field()
    ! ==========================================================================
    subroutine storm_wind_field(t,wind_index,wind,storm)

        implicit none

        ! Time of the wind field requested
        real(kind=8), intent(in) :: t

        ! Storm description, need in out here since we may update the storm
        ! if at next time point
        type(holland_storm_type), intent(in out) :: storm

        ! Array storing wind field
        integer, intent(in) :: wind_index
        real(kind=8), intent(in out) :: wind(:,:,:)

        ! First update time point of storm if necessary
        call update_storm(t,storm)

        ! Find

        wind(wind_index:wind_index+1,:,:) = 0.d0
        
    end subroutine storm_wind_field

    ! ==========================================================================
    !  storm_pressure_field()
    ! ==========================================================================
    subroutine storm_pressure_field(t,pressure_index,pressure)

        implicit none

        ! Time requested
        real(kind=8), intent(in) :: t

        ! Array storing pressure field
        integer, intent(in) :: presssure_index
        real(kind=8), intent(in out) :: pressure(:,:,:)

        if (t > )

        pressure(pressure_index,:,:) = 0.d0
        
    end subroutine storm_pressure_field

end module holland_storm_module