! ==============================================================================
!  Storm Surge Module - Contains generic routines for dealing with storm surge
!    including AMR parameters and storm fields.  This module includes modules
!    for specific implementations of storms such as the Holland model.
! ==============================================================================
!                   Copyright (C) Clawpack Developers 2017
!  Distributed under the terms of the Berkeley Software Distribution (BSD) 
!  license
!                     http://www.opensource.org/licenses/
! ==============================================================================

module storm_module

    use model_storm_module, only: model_storm_type
    use data_storm_module, only: data_storm_type

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
    integer :: storm_specification_type
    real(kind=8) :: landfall = 0.d0
    type(model_storm_type), save :: model_storm
    type(data_storm_type), save :: data_storm

    ! Interface to each of the parameterized models
    abstract interface
        subroutine set_model_fields_def(maux, mbc, mx, my, xlower, ylower, &
                                      dx, dy, t, aux, wind_index,              &
                                      pressure_index, storm)

            use model_storm_module, only: model_storm_type

            implicit none
            integer, intent(in) :: maux,mbc,mx,my
            real(kind=8), intent(in) :: xlower,ylower,dx,dy,t
            type(model_storm_type), intent(inout) :: storm
            integer, intent(in) :: wind_index, pressure_index
            real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
        end subroutine set_model_fields_def
    end interface

    abstract interface
        subroutine set_data_fields_def(maux, mbc, mx, my, xlower, ylower,    &
                                      dx, dy, t, aux, wind_index,           &
                                      pressure_index, storm)

            use data_storm_module, only: data_storm_type

            implicit none
            integer, intent(in) :: maux, mbc, mx, my
            real(kind=8), intent(in) :: xlower, ylower, dx, dy, t
            type(data_storm_type), intent(in out) :: storm
            integer, intent(in) :: wind_index, pressure_index
            real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
        end subroutine set_data_fields_def
    end interface

    procedure(set_model_fields_def), pointer :: set_model_fields
    procedure(set_data_fields_def), pointer :: set_data_fields

    ! Display times in days relative to landfall
    logical :: display_landfall_time = .false.

contains
    ! ========================================================================
    !   subroutine set_storm(data_file)
    ! ========================================================================
    ! Reads in the data file at the path data_file.
    !
    ! Input:
    !     data_file = Path to data file
    !
    ! ========================================================================
    subroutine set_storm(data_file)

        use model_storm_module, only: set_model_storm => set_storm
        use model_storm_module, only: set_holland_1980_fields
        use model_storm_module, only: set_holland_2010_fields
        use model_storm_module, only: set_CLE_fields

        ! use data_storm_module, only: set_data_storm => set_storm
        use data_storm_module, only: set_HWRF_fields

        use utility_module, only: get_value_count

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
            read(unit, *) display_landfall_time
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
            read(unit, "(i1)") storm_specification_type
            read(unit, *) storm_file_path
            
            close(unit)

            ! Print log messages
            open(unit=log_unit, file="fort.surge", status="unknown", action="write")

            write(log_unit, *) "Wind Nesting = ", (wind_refine(i),i=1,size(wind_refine,1))
            write(log_unit, *) "R Nesting = ", (R_refine(i),i=1,size(R_refine,1))
            write(log_unit, *) ""

            write(log_unit, *) "Storm specification = ", storm_specification_type
            write(log_unit, *) "  file = ", storm_file_path

            ! Use parameterized storm model
            if (0 < storm_specification_type .and.              &
                    storm_specification_type <=3) then
                select case(storm_specification_type)
                    case(1) ! Holland 1980 model
                        set_model_fields => set_holland_1980_fields
                    case(2) ! Holland 2010 model
                        set_model_fields => set_holland_2010_fields
                    case(3) ! Chavas, Lin, Emmanuel model
                        set_model_fields => set_CLE_fields
                end select
                call set_model_storm(storm_file_path, model_storm,         &
                                     storm_specification_type, log_unit)
            else if (storm_specification_type > 0) then
                print *, "Storm specification model type ",                &
                            storm_specification_type, "not available."
                stop
            end if

            ! Storm will be set based on a gridded set of data
            if (-1 <= storm_specification_type .and.                    &
                      storm_specification_type < 0) then
                select case(storm_specification_type)
                    case(1) ! HWRF Data
                        set_data_fields => set_HWRF_fields
                end select
            else if (storm_specification_type < 0) then
                print *, "Storm specification data type ",               &
                            storm_specification_type, "not available."
                stop
            end if

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
    real(kind=8) pure function powell_wind_drag(wind_speed, theta)     &
                                         result(wind_drag)

        implicit none

        ! Input
        real(kind=8), intent(in) :: wind_speed, theta

        ! Locals
        real(kind=8) :: weight(3), drag(3)

        weight = 0.d0
        drag = garret_wind_drag_limit(wind_speed, 3.5d-3)

        ! Calculate sector weights
        if (0.d0 <= theta .and. theta <= 40.d0) then
            ! Transition from left (3) to right (1)
            weight(3) = 1.d0 - theta / 40.d0
            weight(1) = theta / 40.d0
        else if (40.d0 <= theta .and. theta <= 130.d0) then
            ! Right (1)
            weight(1) = 1.d0
        else if (130.d0 <= theta .and. theta <= 170.d0) then
            ! Transition from right (1) to rear (2)
            weight(1) = 1.d0 - (theta - 130.d0) / 40.d0
            weight(2) = (theta - 130.d0) / 40.d0
        else if (170.d0 <= theta .and. theta <= 220.d0) then
            ! Rear (2)
            weight(2) = 1.d0
        else if (220.d0 <= theta .and. theta <= 260.d0) then
            ! Transition from rear (2) to left (3)
            weight(2) = 1.d0 - (theta - 220.d0) / 40.d0
            weight(3) = (theta - 220.d0) / 40.d0
        else if (260.d0 <= theta .and. theta <= 360.d0) then
            ! Left (3)
            weight(3) = 1.d0
        endif

        ! Calculate wind drag for sector 1 (right)
        if (wind_speed <= 35.d0) then
            drag(1) = min(2.d-3, drag(1))
        else if (35.d0 <= wind_speed .and. wind_speed <= 45.d0) then
            drag(1) = 2.d-3 + 1.d-3 * (wind_speed - 35.d0) / 10.d0
        else if (45.d0 <= wind_speed) then
            drag(1) = 3.d-3
        endif

        ! Calculate wind drag for sector 2 (rear)
        if (wind_speed <= 35.d0) then
            drag(2) = min(2.d-3, drag(2))
        else if (35.d0 <= wind_speed .and. wind_speed <= 45.d0) then
            drag(2) = 2.d-3 - 1.d-3 * (wind_speed - 35.d0) / 10.d0
        else if (45.d0 <= wind_speed) then
            drag(2) = 1.d-3
        endif

        ! Calculate wind drag for sector 3 (left)
        if (drag(3) > 1.8d-3) then
            if (wind_speed <= 25.d0) then
                drag(3) = 1.8d-3
            else if (25.d0 <= wind_speed .and. wind_speed <= 30.d0) then
                drag(3) = 1.8d-3 + 2.7d-3 * (wind_speed - 25.d0) / 5.d0
            else if (30.d0 <= wind_speed .and. wind_speed <= 45.d0) then
                drag(3) = 4.5d-3 - 3.5d-3 * (wind_speed - 30.d0) / 15.d0
            else if (45.d0 <= wind_speed) then
                drag(3) = 1.d-3
            endif
        endif

        ! Use weighted average
        wind_drag = dot_product(weight, drag)

    end function powell_wind_drag


    ! ========================
    !  Garret Based Wind Drag
    ! ========================
    !  This version is a simple limited version of the wind drag
    real(kind=8) pure function garret_wind_drag(wind_speed, theta) result(wind_drag)
    
        implicit none
        
        ! Input
        real(kind=8), intent(in) :: wind_speed, theta
  
        wind_drag = garret_wind_drag_limit(wind_speed, 2.d-3)
    
    end function garret_wind_drag


    ! ===================================
    !  Helper for Garret Based Wind Drag
    ! ===================================
    real(kind=8) pure function garret_wind_drag_limit(wind_speed,      &
                                                      wind_drag_limit) &
                                               result(wind_drag)

        implicit none

        ! Input
        real(kind=8), intent(in) :: wind_speed, wind_drag_limit

        wind_drag = min(wind_drag_limit,                               &
                        1.d-3 * (0.75d0 + 0.067d0 * wind_speed))

    end function garret_wind_drag_limit


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

        use model_storm_module, only: model_location => storm_location
        use data_storm_module, only: data_location => storm_location

        implicit none

        ! Input
        real(kind=8), intent(in) :: t

        ! Output
        real(kind=8) :: location(2)

        if (storm_specification_type == 0) then
            location = [rinfinity,rinfinity]
        else if (storm_specification_type > 0) then
            location = model_location(t, model_storm)
        else if (storm_specification_type < 0) then
            location = data_location(t, data_storm)
        else
            stop "Something may be wrong."
        end if

    end function storm_location

    real(kind=8) function storm_direction(t) result(theta)
        
        use amr_module, only: rinfinity
        use model_storm_module, only: model_direction => storm_direction 
        use data_storm_module, only: data_direction => storm_direction

        implicit none

        ! Input
        real(kind=8), intent(in) :: t

        if (storm_specification_type > 0) then
            theta = model_direction(t, model_storm)
        else if (storm_specification_type < 0) then
            theta = data_direction(t, data_storm)
        else
            theta = rinfinity
        end if

    end function storm_direction


    subroutine set_storm_fields(maux, mbc, mx, my, xlower, ylower, dx, dy,&
                                t, aux)

        implicit none

        ! Input arguments
        integer, intent(in) :: maux, mbc, mx, my
        real(kind=8), intent(in) :: xlower, ylower, dx, dy, t
        real(kind=8), intent(in out) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

        if (storm_specification_type > 0) then
            call set_model_fields(maux,mbc,mx,my, &
                                  xlower,ylower,dx,dy,t,aux, wind_index, &
                                  pressure_index, model_storm)
        end if
        if (storm_specification_type < 0) then
            call set_data_fields(maux,mbc,mx,my, &
                                 xlower,ylower,dx,dy,t,aux, wind_index, &
                                 pressure_index, data_storm)
        end if

    end subroutine set_storm_fields


    subroutine output_storm_location(t)

        implicit none

        real(kind=8), intent(in) :: t
        
        ! We open this here so that the file flushes and writes to disk
        open(unit=track_unit,file="fort.track",action="write",position='append')

        write(track_unit,"(4e26.16)") t, storm_location(t), storm_direction(t)
        
        close(track_unit)

    end subroutine output_storm_location

end module storm_module
