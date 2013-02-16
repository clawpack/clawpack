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
    type(holland_storm_type) :: holland_storm
    type(constant_storm_type) :: constant_storm
    type(stommel_storm_type) :: stommel_storm

    ! Store physics here for ease of use
    ! WARNING:  If these do not agree with the storm data objects things will break!
    real(kind=8) :: rho_air, ambient_pressure

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
        
        ! Open file
        if (present(data_file)) then
            call opendatafile(unit,data_file)
        else
            call opendatafile(unit,'surge.data')
        endif

        ! Set some parameters
        wind_index = 5
        pressure_index = 7
        
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
            case(3)
                wind_drag => powell_wind_drag
            case default
                stop "*** ERROR *** Invalid wind drag law."
        end select
        read(unit,*) pressure_forcing
        read(unit,*)
        
        ! Set drag law function pointer
!         ! Source term algorithm parameters
!         read(unit,*) wind_tolerance
!         read(unit,*) pressure_tolerance
!         read(unit,*)
        
        ! AMR parameters
        read(unit,'(a)') line
        allocate(wind_refine(get_value_count(line)))
        read(line,*) (wind_refine(i),i=1,size(wind_refine,1))
        read(unit,'(a)') line
        allocate(R_refine(get_value_count(line)))
        read(line,*) (R_refine(i),i=1,size(R_refine,1))
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

    end subroutine set_storm

    ! ========================================================================
    !   real(kind=8) function *_wind_drag(wind_speed)
    ! ========================================================================
    !  Calculates the drag coefficient for wind given the given wind speed.
    !  Based on the modeling from the paper by Weisberg and Zheng (2006).
    !  
    !  Input:
    !      wind_speed = Magnitude of the wind in the cell
    !      theta = Angle with primary hurricane motion
    !
    !  Output:
    !      wind_drag = Coefficient of drag
    ! ==========================================================================
    real(kind=8) pure function powell_wind_drag(wind_speed, theta) result(wind_drag)
    
        implicit none
        
        ! Input
        real(kind=8), intent(in) :: wind_speed, theta

        wind_drag = min(2.d-3, (0.75d0 + 0.067d0 * wind_speed) * 1d-3)
    
    end function powell_wind_drag

    ! This version ignores direction
    real(kind=8) pure function garret_wind_drag(wind_speed, theta) result(wind_drag)
    
        implicit none
        
        ! Input
        real(kind=8), intent(in) :: wind_speed, theta
  
        wind_drag = min(2.d-3, (0.75d0 + 0.067d0 * wind_speed) * 1d-3)      
!         if (wind_speed <= 11.d0) then
!             wind_drag = 1.2d0
!         else if ((wind_speed > 11.d0).and.(wind_speed <= 25.d0)) then
!             wind_drag = 0.49d0 + 0.065d0 * wind_speed
!         else
!             wind_drag = 0.49 + 0.065d0 * 25.d0
!         endif
        
!         wind_drag = wind_drag * 1.d-3
    
    end function garret_wind_drag

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


    subroutine set_storm_fields(maxmx,maxmy,maux,mbc,mx,my,xlower,ylower,dx,dy,&
                                t,aux)

        use holland_storm_module, only: set_holland_storm_fields
        use constant_storm_module, only: set_constant_storm_fields
        use stommel_storm_module, only: set_stommel_storm_fields

        implicit none

        ! Input arguments
        integer, intent(in) :: maxmx, maxmy, maux, mbc, mx, my
        real(kind=8), intent(in) :: xlower, ylower, dx, dy, t
        real(kind=8), intent(in out) :: aux(maux,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc)

        select case(storm_type)
            case(0)
                continue
            case(1)
                call set_holland_storm_fields(maxmx,maxmy,maux,mbc,mx,my, &
                                    xlower,ylower,dx,dy,t,aux, wind_index, &
                                    pressure_index, holland_storm)
            case(2)
                call set_constant_storm_fields(maxmx,maxmy,maux,mbc,mx,my, &
                                    xlower,ylower,dx,dy,t,aux, wind_index, &
                                    pressure_index, constant_storm)
            case(3)
                call set_stommel_storm_fields(maxmx,maxmy,maux,mbc,mx,my, &
                                    xlower,ylower,dx,dy,t,aux, wind_index, &
                                    pressure_index, stommel_storm)
        end select

    end subroutine set_storm_fields

    subroutine output_storm_location(t)

        implicit none

        real(kind=8), intent(in) :: t
        
        ! We open this here so that the file flushes and writes to disk
        open(unit=track_unit,file="fort.track",action="write",position='append')

        write(track_unit,"(3e26.16)") t,storm_location(t)
        
        close(track_unit)

    end subroutine output_storm_location

end module storm_module
