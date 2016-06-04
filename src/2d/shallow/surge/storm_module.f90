! ==============================================================================
!  Storm Surge Module - Contains generic routines for dealing with storm surge
!    including AMR parameters and storm fields.  This module includes modules
!    for specific implementations of storms such as the Holland model.
! ==============================================================================
!  Distributed under the terms of the Berkeley Software Distribution (BSD) 
!  license
!                     http://www.opensource.org/licenses/
! ==============================================================================

module storm_module

    use holland_storm_module, only: holland_storm_type
    use constant_storm_module, only: constant_storm_type
    use stommel_storm_module, only: stommel_storm_type

    implicit none

    logical, private :: module_setup = .false.

    ! Log file IO unit
    integer, parameter :: log_unit = 423

    ! Track file IO unit
    integer, parameter :: track_unit = 424

    ! Locations of wind and pressure fields
    integer :: wind_index, pressure_index
    
    ! Source term control and parameters
    logical :: wind_forcing, pressure_forcing

    ! Wind drag law support
    abstract interface
        real(kind=8) pure function drag_function(speed, theta)
            implicit none
            real(kind=8), intent(in) :: speed, theta
        end function drag_function
    end interface
        
    ! Function pointer to wind drag requested
    procedure (drag_function), pointer :: wind_drag

    ! AMR Parameters
    real(kind=8), allocatable :: R_refine(:), wind_refine(:)

    ! Storm object
    integer :: storm_type
    real(kind=8) :: landfall = 0.d0
    real(kind=8) :: ramp_width
    type(holland_storm_type), save :: holland_storm
    type(constant_storm_type), save :: constant_storm
    type(stommel_storm_type), save :: stommel_storm

    ! Store physics here for ease of use
    ! WARNING:  If these do not agree with the storm data objects things will break!
    real(kind=8) :: rho_air, ambient_pressure

    ! Wind drag maximum limit
    real(kind=8), parameter :: WIND_DRAG_LIMIT = 2.d-3

contains
    ! ========================================================================
    !   subroutine set_storm(data_file)
    ! ========================================================================
    ! Reads in the data file at the path data_file.  Sets the following 
    ! parameters:
    !     rho_air = density of air
    !     Pn = Ambient atmospheric pressure
    !
    ! Input:
    !     data_file = Path to data file
    !
    ! ========================================================================
    subroutine set_storm(data_file)

        use utility_module, only: get_value_count

        use holland_storm_module, only: set_holland_storm
        use constant_storm_module, only: set_constant_storm
        use stommel_storm_module, only: set_stommel_storm

        use geoclaw_module, only: pi

        implicit none
        
        ! Input arguments
        character(len=*), optional, intent(in) :: data_file
        
        ! Locals
        integer, parameter :: unit = 13
        integer :: i, drag_law
        character(len=200) :: storm_file_path, line
        
        if (.not.module_setup) then

            ! Open file
            if (present(data_file)) then
                call opendatafile(unit,data_file)
            else
                call opendatafile(unit,'surge.data')
            endif
            
            ! Read in parameters
            ! Physics
            ! TODO: These are currently set directly in the types, should change!
            read(unit,"(1d16.8)") rho_air
            read(unit,"(1d16.8)") ambient_pressure

            ! Forcing terms
            read(unit,*) wind_forcing
            read(unit,*) drag_law
            if (.not.wind_forcing) drag_law = 0
            select case(drag_law)
                case(0)
                    wind_drag => no_wind_drag
                case(1)
                    wind_drag => garret_wind_drag
                case(2)
                    wind_drag => powell_wind_drag
                case default
                    stop "*** ERROR *** Invalid wind drag law."
            end select
            read(unit,*) pressure_forcing
            read(unit,*)

            ! Set some parameters
            read(unit, '(i2)') wind_index
            read(unit, '(i2)') pressure_index
            read(unit, *)
            
            ! AMR parameters
            read(unit,'(a)') line
            if (line(1:1) == "F") then
                allocate(wind_refine(1))
                wind_refine(1) = huge(1.d0)
            else
                allocate(wind_refine(get_value_count(line)))
                read(line,*) (wind_refine(i),i=1,size(wind_refine,1))
            end if
            read(unit,'(a)') line
            if (line(1:1) == "F") then
                allocate(R_refine(1))
                R_refine(1) = -huge(1.d0)
            else
                allocate(R_refine(get_value_count(line)))
                read(line,*) (R_refine(i),i=1,size(R_refine,1))
            end if
            read(unit,*)
            
            ! Storm Setup 
            read(unit,"(i1)") storm_type
            read(unit,"(d16.8)") landfall
            read(unit,*) storm_file_path
            
            close(unit)

            ! Print log messages
            open(unit=log_unit, file="fort.surge", status="unknown", action="write")
            
            write(log_unit,"(a,1d16.8)") "rho_air =          ",rho_air
            write(log_unit,"(a,1d16.8)") "ambient_pressure = ",ambient_pressure
            write(log_unit,*) ""

    !         write(log_unit,"(a,1d16.8)") "wind_tolerance =     ", wind_tolerance
    !         write(log_unit,"(a,1d16.8)") "pressure_tolerance = ", pressure_tolerance
    !         write(log_unit,*) ""

            write(log_unit,*) "Wind Nesting = ", (wind_refine(i),i=1,size(wind_refine,1))
            write(log_unit,*) "R Nesting = ", (R_refine(i),i=1,size(R_refine,1))
            write(log_unit,*) ""

            write(log_unit,*) "Storm Type = ", storm_type
            write(log_unit,*) "  file = ", storm_file_path

            ! Read in hurricane track data from file
            if (storm_type == 0) then
                ! No storm will be used
            else if (storm_type == 1) then
                ! Track file with Holland reconstruction
                call set_holland_storm(storm_file_path,holland_storm,log_unit)
                ! Set rho_air and ambient pressure in storms
            else if (storm_type == 2) then
                ! constant track and holland reconstruction
                call set_constant_storm(storm_file_path,constant_storm,log_unit)
                ! Set rho_air and ambient pressure in storms
            else if (storm_type == 3) then
                ! Stommel wind field
                call set_stommel_storm(storm_file_path,stommel_storm,log_unit)
            else
                print *,"Invalid storm type ",storm_type," provided."
                stop
            endif

            close(log_unit)

            module_setup = .true.
        end if

    end subroutine set_storm

    ! ========================================================================
    !   real(kind=8) function *_wind_drag(wind_speed)
    ! ========================================================================
    !  Calculates the drag coefficient for wind given the given wind speed.
    !  
    !  Input:
    !      wind_speed = Magnitude of the wind in the cell
    !      theta = Angle with primary hurricane direciton
    !
    !  Output:
    !      wind_drag = Coefficient of drag
    ! ==========================================================================
    !  Powell Wind Drag - Sector based wind drag coefficients due to primary
    !    wave direction interaction with wind.  This implementation is based on
    !    the parameterization used in ADCIRC.  For more information see
    !
    !    M.D. Powell (2006). “Final Report to the National Oceanic and 
    !      Atmospheric Administration (NOAA) Joint Hurricane Testbed (JHT) 
    !      Program.” 26 pp.
    !
    real(kind=8) pure function powell_wind_drag(wind_speed, theta)      &
                                         result(wind_drag)
    
        implicit none
        
        ! Input
        real(kind=8), intent(in) :: wind_speed, theta

        ! Locals
        real(kind=8) :: weight, left_drag, right_drag, rear_drag, drag(2)

        ! Calculate sector drags
        if (wind_speed <= 15.708d0) then
            left_drag = 7.5d-4 + 6.6845d-05 * wind_speed
            right_drag = left_drag
            rear_drag = left_drag
        else if (15.708d0 < wind_speed .and. wind_speed <= 18.7d0) then
            left_drag = 1.8d-3
            right_drag = 7.5d-4 + 6.6845d-05 * wind_speed
            rear_drag = right_drag
        else if (18.7d0 < wind_speed .and. wind_speed <= 25.d0) then
            left_drag = 1.8d-3
            right_drag = 2.0d-3
            rear_drag = right_drag
        else if (25.0d0 < wind_speed .and. wind_speed <= 30.d0) then
            left_drag = 1.8d-3 + 5.4d-3 * (wind_speed - 25.d0)
            right_drag = 2.0d-3
            rear_drag = right_drag
        else if (30.0d0 < wind_speed .and. wind_speed <= 35.d0) then
            left_drag = 4.5d-3 - 2.33333d-4 * (wind_speed - 30.d0)
            right_drag = 2.0d-3
            rear_drag = right_drag
        else if (35.0d0 < wind_speed .and. wind_speed <= 45.d0) then
            left_drag = 4.5d-3 - 2.33333d-4 * (wind_speed - 30.d0)
            right_drag = 2.d-3 + 1.d-4 * (wind_speed - 35.d0)
            rear_drag = 2.d-3 - 1.d-4 * (wind_speed - 35.d0)
        else
            left_drag = 1.0d-3
            right_drag = 3.0d-3
            rear_drag = left_drag
        endif

        ! Calculate weights
        ! Left sector =  [240.d0,  20.d0] - Center = 310
        ! Left Right sector = [310, 85] - Width = 145
        ! Right sector = [ 20.d0, 150.d0] - Center = 85
        ! Right Rear sector = [85, 195] - Width = 
        ! Rear sector =  [150.d0, 240.d0] - 195
        ! Rear-Left sector = [85, 195] - Width = 

        ! Left - Right sector
        if (310.d0 < theta .and. theta <= 360.d0) then
            weight = (theta - 310.d0) / 135.d0
            drag = [left_drag, right_drag]

        ! Left - Right sector
        else if (0.d0 < theta .and. theta <= 85.d0) then
            weight = (theta + 50.d0) / 135.d0
            drag = [left_drag, right_drag]

        ! Right - Rear sector
        else if (85.d0 < theta .and. theta <= 195.d0) then
            weight = (theta - 85.d0) / 110.d0
            drag = [right_drag, rear_drag]

        ! Rear - Left sector
        else if (195.d0 < theta .and. theta <= 310.d0) then
            weight = (theta - 195.d0) / 115.d0
            drag = [rear_drag, left_drag]
        endif

        wind_drag = drag(1) * (1.d0 - weight) + drag(2) * weight

        ! Apply wind drag limit - May want to do this...
        ! wind_drag = min(WIND_DRAG_LIMIT, wind_drag)
    
    end function powell_wind_drag


    ! ========================
    !  Garret Based Wind Drag
    ! ========================
    !  This version is a simple limited version of the wind drag
    real(kind=8) pure function garret_wind_drag(wind_speed, theta) result(wind_drag)
    
        implicit none
        
        ! Input
        real(kind=8), intent(in) :: wind_speed, theta
  
        wind_drag = min(WIND_DRAG_LIMIT, (0.75d0 + 0.067d0 * wind_speed) * 1d-3)      
    
    end function garret_wind_drag


    ! ==================================================================
    !  No Wind Drag - Dummy function used to turn off wind drag forcing
    ! ==================================================================
    real(kind=8) pure function no_wind_drag(wind_speed, theta) result(wind_drag)
        implicit none
        real(kind=8), intent(in) :: wind_speed, theta
        wind_drag = 0.d0
    end function no_wind_drag



    ! ==========================================================================
    ! Wrapper functions for all storm types
    function storm_location(t) result(location)

        use amr_module, only: rinfinity

        use holland_storm_module, only: holland_storm_location
        use constant_storm_module, only: constant_storm_location

        implicit none

        ! Input
        real(kind=8), intent(in) :: t

        ! Output
        real(kind=8) :: location(2)

        select case(storm_type)
            case(0)
                location = [rinfinity,rinfinity]
            case(1)
                location = holland_storm_location(t,holland_storm)
            case(2)
                location = constant_storm_location(t,constant_storm)
            case(3)
                location = [rinfinity,rinfinity]
        end select

    end function storm_location

    real(kind=8) function storm_direction(t) result(theta)
        
        use amr_module, only: rinfinity
        use holland_storm_module, only: holland_storm_direction
        use constant_storm_module, only: constant_storm_direction

        implicit none

        ! Input
        real(kind=8), intent(in) :: t

        select case(storm_type)
            case(0)
                theta = rinfinity
            case(1)
                theta = holland_storm_direction(t,holland_storm)
            case(2)
                theta = constant_storm_direction(t,constant_storm)
            case(3)
                theta = rinfinity
        end select

    end function storm_direction


    subroutine set_storm_fields(maux,mbc,mx,my,xlower,ylower,dx,dy,&
                                t,aux)

        use holland_storm_module, only: set_holland_storm_fields
        use constant_storm_module, only: set_constant_storm_fields
        use stommel_storm_module, only: set_stommel_storm_fields

        implicit none

        ! Input arguments
        integer, intent(in) :: maux, mbc, mx, my
        real(kind=8), intent(in) :: xlower, ylower, dx, dy, t
        real(kind=8), intent(in out) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

        select case(storm_type)
            case(0)
                continue
            case(1)
                call set_holland_storm_fields(maux,mbc,mx,my, &
                                    xlower,ylower,dx,dy,t,aux, wind_index, &
                                    pressure_index, holland_storm)
            case(2)
                call set_constant_storm_fields(maux,mbc,mx,my, &
                                    xlower,ylower,dx,dy,t,aux, wind_index, &
                                    pressure_index, constant_storm)
            case(3)
                call set_stommel_storm_fields(maux,mbc,mx,my, &
                                    xlower,ylower,dx,dy,t,aux, wind_index, &
                                    pressure_index, stommel_storm)
        end select

    end subroutine set_storm_fields

    subroutine output_storm_location(t)

        implicit none

        real(kind=8), intent(in) :: t
        
        ! We open this here so that the file flushes and writes to disk
        open(unit=track_unit,file="fort.track",action="write",position='append')

        write(track_unit,"(4e26.16)") t,storm_location(t),storm_direction(t)
        
        close(track_unit)

    end subroutine output_storm_location

end module storm_module
