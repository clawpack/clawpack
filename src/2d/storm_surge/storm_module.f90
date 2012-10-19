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

    use geoclaw_module, only: num_layers

    use holland_storm_module, only: holland_storm_type
    use constant_storm_module, only: constant_storm_type
    use stommel_storm_module, only: stommel_storm_type

    implicit none

    ! Log file unit
    integer, parameter :: log_unit = 423

    ! Track file unit
    integer, parameter :: track_unit = 424

    ! Locations of wind and pressure fields
    integer :: wind_index, pressure_index
    
    ! Source term control and parameters
    logical :: wind_forcing, pressure_forcing
    integer :: wind_type
    real(kind=8) :: pressure_tolerance, wind_tolerance

    ! AMR Parameters
    real(kind=8), allocatable :: R_refine(:), wind_refine(:)

    ! Storm object
    integer :: storm_type
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
        integer :: i
        character(len=200) :: storm_file_path, line
        
        ! Open file
        if (present(data_file)) then
            call opendatafile(unit,data_file)
        else
            call opendatafile(unit,'surge.data')
        endif

        ! Set some parameters
        wind_index = 4 + num_layers
        pressure_index = 6 + num_layers
        
        ! Read in parameters
        ! Physics
        ! TODO: These are currently set directly in the types, should change!
        read(unit,"(1d16.8)") rho_air
        read(unit,"(1d16.8)") ambient_pressure

        ! Forcing terms
        read(unit,*) wind_forcing
        read(unit,*) pressure_forcing
        read(unit,*)
        
        ! Source term algorithm parameters
        read(unit,*) wind_tolerance
        read(unit,*) pressure_tolerance
        read(unit,*)
        
        ! AMR parameters
        read(unit,*) line
        allocate(wind_refine(get_value_count(line)))
        read(line,*) (wind_refine(i),i=1,size(wind_refine,1))
        read(unit,*) line
        allocate(R_refine(get_value_count(line)))
        read(line,*) (R_refine(i),i=1,size(R_refine,1))
        read(unit,*)
        
        ! Storm Setup 
        read(unit,"(i1)") storm_type
        read(unit,*) storm_file_path

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
        
        close(unit)

        ! Print log messages
        open(unit=log_unit, file="fort.surge", status="unknown", action="write")
        
        write(log_unit,"(a,1d16.8)") "rho_air =          ",rho_air
        write(log_unit,"(a,1d16.8)") "ambient_pressure = ",ambient_pressure
        write(log_unit,*) ""

        write(log_unit,"(a,1d16.8)") "wind_tolerance =     ", wind_tolerance
        write(log_unit,"(a,1d16.8)") "pressure_tolerance = ", pressure_tolerance
        write(log_unit,*) ""

        write(log_unit,*) "Wind Nesting = ", (wind_refine(i),i=1,size(wind_refine,1))
        write(log_unit,*) "R Nesting = ", (R_refine(i),i=1,size(R_refine,1))
        write(log_unit,*) ""

        write(log_unit,*) "Storm Type = ", storm_type
        write(log_unit,*) "  file = ", storm_file_path

        ! Open track output file if using storm type 1
        if (storm_type == 1) then
            open(unit=track_unit,file="fort.track",status="unknown",action="write")
        endif

    end subroutine set_storm

    ! ========================================================================
    !   double precision function wind_drag(wind_speed)
    ! ========================================================================
    !  Calculates the drag coefficient for wind given the given wind speed.
    !  Based on the modeling from the paper by Weisberg and Zheng (2006).
    !  
    !  Input:
    !      wind_speed = Magnitude of the wind in the cell
    !
    !  Output:
    !      wind_drag = Coefficient of drag
    ! ==========================================================================
    real(kind=8) pure function wind_drag(wind_speed)
    
        implicit none
        
        ! Input
        real(kind=8), intent(in) :: wind_speed
        
        if (wind_speed <= 11.d0) then
            wind_drag = 1.2d0
        else if ((wind_speed > 11.d0).and.(wind_speed <= 25.d0)) then
            wind_drag = 0.49d0 + 0.065d0 * wind_speed
        else
            wind_drag = 0.49 + 0.065d0 * 25.d0
        endif
        
        wind_drag = wind_drag * 10.d-3
    
    end function wind_drag

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
        
        write(track_unit,"(3e26.16)") t,storm_location(t)

    end subroutine output_storm_location

end module storm_module
