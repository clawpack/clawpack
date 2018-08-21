! ==============================================================================
! model_storm_module 
!
! Module contains routines for constructing a wind and pressure field based on
! the a parameterized model of the wind and pressure fields.  
! 
! ==============================================================================
!                   Copyright (C) Clawpack Developers 2017
!  Distributed under the terms of the Berkeley Software Distribution (BSD) 
!  license
!                     http://www.opensource.org/licenses/
! ==============================================================================
module model_storm_module

    implicit none
    save

    logical, private :: module_setup = .false.
    logical, private :: DEBUG = .false.

    ! Model storm type definition
    type model_storm_type
        ! Fore/hindcast size and current position
        integer :: num_casts

        ! Landfall - This is not used explicitly (t0 = landfall ideally)
        character(len=10) :: landfall

        ! These parameters are located at time points but are interpolated in
        ! time and space when the relevant fields are requested.

        ! Location of storm
        ! Track is a triplet with (time,longitude,latitude)
        real(kind=8), allocatable :: track(:,:)

        ! Storm parameterization
        real(kind=8), allocatable :: max_wind_radius(:)
        real(kind=8), allocatable :: max_wind_speed(:)
        real(kind=8), allocatable :: central_pressure(:)

        ! This is not always provided and is often determined by the last
        ! closed iso-bar of the storm
        real(kind=8), allocatable :: radius(:)

        ! Approximate velocity of storm, approximated via the track points
        ! using a first order difference on the sphere
        real(kind=8), allocatable :: velocity(:, :)

    end type model_storm_type

    ! Interal tracking variables for storm
    integer, private :: last_storm_index

    ! Atmospheric boundary layer, input variable in ADCIRC but always is
    ! set to the following value
    real(kind=8), parameter :: atmos_boundary_layer = 0.9d0

    ! Sampling adjustment from 1 min to 10 min winds
    real(kind=8), parameter :: sampling_time = 0.88d0 

    ! Storm field ramping width - Represents crudely the ramping radial area
    real(kind=8), parameter :: RAMP_WIDTH = 100.0d3

    ! Time tracking tolerance allowance - allows for the beginning of the storm
    ! track to be close to but not equal the start time of the simulation
    real(kind=8), parameter :: TRACKING_TOLERANCE = 1d-10

contains


    ! Setup routine for model storms
    subroutine set_storm(storm_data_path, storm, storm_spec_type, log_unit)

        use geoclaw_module, only: deg2rad, spherical_distance, coordinate_system
        use amr_module, only: t0, rinfinity

        implicit none

        ! Subroutine I/O
        character(len=*), optional :: storm_data_path
        type(model_storm_type), intent(inout) :: storm
        integer, intent(in) :: storm_spec_type, log_unit

        ! Local storage
        integer, parameter :: data_file = 701
        integer :: i, k, io_status
        real(kind=8) :: x(2), y(2), ds, dt

        if (.not. module_setup) then

            ! Open data file
            print *,'Reading storm date file ', storm_data_path
            open(unit=data_file, file=storm_data_path, status='old',        &
                 action='read', iostat=io_status)
            if (io_status /= 0) then
                print *, "Error opening storm data file. status = ", io_status
                stop 
            endif

            read(data_file, "(i4)") storm%num_casts
            read(data_file, "(a)") storm%landfall
            read(data_file, *)

            write(log_unit, "('Data length = ',i3)") storm%num_casts

            ! Allocate storm parameter file variables
            allocate(storm%track(3, storm%num_casts))
            allocate(storm%max_wind_speed(storm%num_casts))
            allocate(storm%max_wind_radius(storm%num_casts))
            allocate(storm%central_pressure(storm%num_casts))
            allocate(storm%radius(storm%num_casts))

            ! Now read in the storm data - note that the units are expected to 
            ! be consistent with:
            ! max_wind_speed = m/s
            ! max_wind_radius = m
            ! central_pressure = Pa
            ! radius = m
            do i=1, storm%num_casts
                read(data_file, *) storm%track(:, i), &
                                   storm%max_wind_speed(i), &
                                   storm%max_wind_radius(i), &
                                   storm%central_pressure(i), &
                                   storm%radius(i)
            enddo

            ! Calculate storm speed 
            allocate(storm%velocity(2, storm%num_casts))
            do i=1,storm%num_casts - 1
                ! Calculate velocity based on great circle distance between

                ! locations of storm
                x = storm%track(2:3,i)
                y = storm%track(2:3,i+1)

                dt = storm%track(1,i + 1) - storm%track(1,i)

                if (coordinate_system == 2) then
                    ds = spherical_distance(x(1), 0.5d0 * (x(2) + y(2)), &
                                            y(1), 0.5d0 * (x(2) + y(2)))
                    storm%velocity(1,i) = sign(ds / dt, y(1) - x(1))
                
                    ds = spherical_distance(0.5d0 * (x(1) + y(1)), x(2), &
                                            0.5d0 * (x(1) + y(1)), y(2))
                    storm%velocity(2, i) = sign(ds / dt, y(2) - x(2))
                else
                    storm%velocity(1, i) = abs(x(2) - x(1)) / dt
                    storm%velocity(2, i) = abs(y(2) - y(1)) / dt
                end if
            end do

            ! Use last approximation for velocity point going forward
            storm%velocity(:, storm%num_casts) = storm%velocity(:,  &
                                                            storm%num_casts - 1)

            if (t0 <= storm%track(1, 1) - TRACKING_TOLERANCE) then
                print *, "Start time", t0, " is outside of the tracking"
                print *, "tolerance range with the track start"
                print *, storm%track(1, 1), "."
                stop
            endif

            ! This is used to speed up searching for correct storm data
            last_storm_index = 2
            last_storm_index = storm_index(t0, storm)
            if (last_storm_index == -1) then
                print *,"Forecast not found for time ",t0,'.'
                stop
            endif

            ! Log everything to the surge log file
            write(log_unit,*) ""
            write(log_unit,*) "Storm Track and Strength"
            write(log_unit,*) ""
            do i=1, storm%num_casts
                write(log_unit,"(8e26.16)") (storm%track(k,i),k=1,3),  &
                                            (storm%velocity(k,i),k=1,2), &
                                             storm%max_wind_radius(i), &
                                             storm%max_wind_speed(i),  &
                                             storm%central_pressure(i)
            enddo

            module_setup = .true.
        end if

    end subroutine set_storm

    ! ==========================================================================
    !  storm_location(t,storm)
    !    Interpolate location of hurricane in the current time interval
    ! ==========================================================================
    function storm_location(t, storm) result(location)

        implicit none

        ! Input
        real(kind=8), intent(in) :: t
        type(model_storm_type), intent(in out) :: storm

        ! Output
        real(kind=8) :: location(2)

        ! Junk storage
        real(kind=8) :: junk(6)

        call get_storm_data(t, storm, location,                 &
                                  junk(1:2), junk(3), junk(4), junk(5), junk(6))

    end function storm_location

    ! ==========================================================================
    !  storm_direction
    !   Angle off of due north that the storm is traveling
    ! ==========================================================================
    real(kind=8) function storm_direction(t, storm) result(theta)

        implicit none

        ! Input
        real(kind=8), intent(in) :: t
        type(model_storm_type), intent(in) :: storm

        ! Locals
        real(kind=8) :: junk(6), velocity(2)

        ! Fetch velocity of storm which has direction encoded in it
        call get_storm_data(t, storm, junk(1:2), velocity, junk(3), junk(4),   &
                                                      junk(5), junk(6))

        ! Unit directional vector
        theta = atan2(velocity(2),velocity(1))

    end function storm_direction

    ! ==========================================================================
    !  storm_index(t,storm)
    !    Finds the index of the next storm data point
    ! ==========================================================================
    integer pure function storm_index(t, storm) result(index)

        implicit none

        ! Input
        real(kind=8), intent(in) :: t
        type(model_storm_type), intent(in) :: storm

        ! Locals
        real(kind=8) :: t0,t1
        logical :: found

        ! Figure out where we are relative to the last time we checked for the
        ! index (stored in last_storm_index)
        
        ! Check if we are already beyond the end of the last forecast time
        if (last_storm_index == storm%num_casts + 1) then
            index = storm%num_casts + 1
        else
            t0 = storm%track(1,last_storm_index - 1)
            t1 = storm%track(1,last_storm_index)

            ! Check to see if we are close enough to the current index to just
            ! use that, tolerance is based on TRACKING_TOLERANCE
            if ((abs(t0 - t) < TRACKING_TOLERANCE) .or.   &
                (abs(t1 - t) < TRACKING_TOLERANCE) .or.   &
                (t0 < t .and. t < t1)) then
                
                index = last_storm_index
            else if ( t1 < t ) then
                found = .false.
                do index=last_storm_index+1,storm%num_casts
                    if (t < storm%track(1,index)) then
                        found = .true.
                        exit
                    endif
                enddo
                ! Assume we have gone past last forecast time
                if (.not. found) then
                    index = storm%num_casts + 1
                endif
            else
                ! Fail gracefully
                if (last_storm_index == 2) then
                    index = -1
                else
                    do index=last_storm_index-1,2,-1
                        if (storm%track(1,index-1) < t) exit
                    enddo
                endif
            endif
        endif

    end function storm_index


    ! ==========================================================================
    !  get_storm_data()
    !    Interpolates in time and returns storm data.
    ! ==========================================================================
    subroutine get_storm_data(t, storm, location, velocity, max_wind_radius,  &
                              max_wind_speed, central_pressure,               &
                              radius)

        use geoclaw_module, only: deg2rad, latlon2xy, xy2latlon

        implicit none

        ! Input
        real(kind=8), intent(in) :: t                 ! Current time
        type(model_storm_type), intent(in) :: storm   ! Storm

        ! Output
        real(kind=8), intent(out) :: location(2), velocity(2)
        real(kind=8), intent(out) :: max_wind_radius, max_wind_speed
        real(kind=8), intent(out) :: central_pressure, radius

        ! Local
        real(kind=8) :: fn(8), fnm(8), weight, tn, tnm, x(2)
        integer :: i

        ! Increment storm data index if needed and not at end of forecast
        i = storm_index(t,storm)
        last_storm_index = i

        ! List of possible error conditions
        if (i <= 1) then
            if (i == 0) then        
                print *,"Invalid storm forecast requested for t = ",t
                print *,"Time requested is before any forecast data."
                print *,"    first time = ",storm%track(1,1)
                print *,"   ERROR = ",i
                stop
            else if (i > storm%num_casts + 2) then
                print *,"Invalid storm indexing, i > num_casts + 2..."
                print *,"This really should not happen, what have you done?"
                stop
            endif
        endif

        ! Interpolate in time for all parameters
        if (i == storm%num_casts + 1) then
            i = i - 1
            ! At last forecast, use last data for storm strength parameters and
            ! velocity, location uses last velocity and constant motion forward

            ! Convert coordinates temporarily to meters so that we can use
            ! the pre-calculated m/s velocities from before
            x = latlon2xy(storm%track(2:3,i),storm%track(2:3,i))
            x = x + (t - storm%track(1,i)) * storm%velocity(:,i)
            
            fn = [xy2latlon(x,storm%track(2:3,i)), &
                  storm%velocity(:,i), storm%max_wind_radius(i), &
                  storm%max_wind_speed(i), storm%central_pressure(i), &
                  storm%radius(i)]
        else
            ! Inbetween two forecast time points (the function storm_index
            ! ensures that we are not before the first data point, i.e. i > 1)
            tn = storm%track(1,i)
            tnm = storm%track(1,i-1)
            weight = (t - tnm) / (tn - tnm)
            fn = [storm%track(2:3,i),storm%velocity(:,i), &
                  storm%max_wind_radius(i),storm%max_wind_speed(i), &
                  storm%central_pressure(i), storm%radius(i)]
            fnm = [storm%track(2:3,i - 1),storm%velocity(:,i - 1), &
                   storm%max_wind_radius(i - 1),storm%max_wind_speed(i - 1), &
                   storm%central_pressure(i - 1), storm%radius(i - 1)]
            fn = weight * (fn - fnm) + fnm
        endif

        ! Set output variables
        location = fn(1:2)
        velocity = fn(3:4)
        max_wind_radius = fn(5)
        max_wind_speed = fn(6)
        central_pressure = fn(7)
        radius = fn(8)

    end subroutine get_storm_data


    ! ==========================================================================
    !  Use the 1980 Holland model to set the storm fields
    ! ==========================================================================
    subroutine set_holland_1980_fields(maux, mbc, mx, my, xlower, ylower,    &
                                       dx, dy, t, aux, wind_index,           &
                                       pressure_index, storm)

        use geoclaw_module, only: g => grav, rho_air, ambient_pressure
        use geoclaw_module, only: coriolis, deg2rad
        use geoclaw_module, only: spherical_distance

        use geoclaw_module, only: rad2deg

        implicit none

        ! Time of the wind field requested
        integer, intent(in) :: maux,mbc,mx,my
        real(kind=8), intent(in) :: xlower,ylower,dx,dy,t

        ! Storm description, need in out here since we may update the storm
        ! if at next time point
        type(model_storm_type), intent(inout) :: storm

        ! Array storing wind and presure field
        integer, intent(in) :: wind_index, pressure_index
        real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

        ! Local storage
        real(kind=8) :: x, y, r, theta, sloc(2), B
        real(kind=8) :: f, mwr, mws, Pc, Pa, dp, wind, tv(2), radius
        real(kind=8) :: mod_mws, trans_speed, ramp
        integer :: i,j

        ! Get interpolated storm data
        call get_storm_data(t, storm, sloc, tv, mwr, mws, Pc, radius)

        ! Other quantities of interest
        Pa = ambient_pressure

        ! Calculate Holland parameters
        ! Subtract translational speed of storm from maximum wind speed
        ! to avoid distortion in the Holland curve fit.  Added back later
        trans_speed = sqrt(tv(1)**2 + tv(2)**2)
        mod_mws = mws - trans_speed

        ! Convert wind speed (10 m) to top of atmospheric boundary layer
        mod_mws = mod_mws / atmos_boundary_layer
        
        ! Calculate central pressure difference
        dp = Pa - Pc
        ! Limit central pressure deficit due to bad ambient pressure,
        ! really should have better ambient pressure...
        if (dp < 100.d0) dp = 100.d0

        ! Calculate Holland parameters and limit the result
        B = rho_air * exp(1.d0) * (mod_mws**2) / dp
        if (B <  1.d0) B = 1.d0
        if (B > 2.5d0) B = 2.5d0

        if (DEBUG) print "('Holland B = ',d16.8)", B
        if (DEBUG) print "('Holland A = ',d16.8)", (mwr / 1000.d0)**B
        
        ! Set fields
        do j=1-mbc,my+mbc
            y = ylower + (j-0.5d0) * dy     ! Degrees latitude
            f = coriolis(y)
            do i=1-mbc,mx+mbc
                x = xlower + (i-0.5d0) * dx   ! Degrees longitude

                ! Calculate storm centric polar coordinate location of grid
                ! cell center, uses Haversine formula
                r = spherical_distance(x, y, sloc(1), sloc(2))
                theta = atan2((y - sloc(2)) * DEG2RAD,(x - sloc(1)) * DEG2RAD)

                ! Set pressure field
                aux(pressure_index,i,j) = Pc + dp * exp(-(mwr / r)**B)

                ! Speed of wind at this point
                wind = sqrt((mwr / r)**B &
                        * exp(1.d0 - (mwr / r)**B) * mws**2.d0 &
                        + (r * f)**2.d0 / 4.d0) - r * f / 2.d0
                
                ! Convert wind velocity from top of atmospheric boundary layer
                ! (which is what the Holland curve fit produces) to wind
                ! velocity at 10 m above the earth's surface

                ! Also convert from 1 minute averaged winds to 10 minute
                ! averaged winds
                wind = wind * atmos_boundary_layer * sampling_time

                ! Velocity components of storm (assumes perfect vortex shape)
                aux(wind_index,i,j)   = -wind * sin(theta)
                aux(wind_index+1,i,j) =  wind * cos(theta)

                ! Add the storm translation speed
                ! Determine translation speed that should be added to final
                ! storm wind speed.  This is tapered to zero as the storm wind
                ! tapers to zero toward the eye of the storm and at long
                ! distances from the storm
                aux(wind_index,i,j) = aux(wind_index,i,j)                 &
                                                    + (abs(wind) / mws) * tv(1)
                aux(wind_index+1,i,j) = aux(wind_index+1,i,j)             &
                                                    + (abs(wind) / mws) * tv(2)

                ! Apply distance ramp down(up) to fields to limit scope
                ramp = 0.5d0 * (1.d0 - tanh((r - radius) / RAMP_WIDTH))
                aux(pressure_index,i,j) = Pa + (aux(pressure_index,i,j) - Pa) &
                                        * ramp
                aux(wind_index:wind_index+1,i,j) =                        &
                                        aux(wind_index:wind_index+1,i,j)  &
                                        * ramp

            enddo
        enddo

    end subroutine set_holland_1980_fields


    ! ==========================================================================
    !  Use the 2010 Holland model to set the storm fields
    ! ==========================================================================
    subroutine set_holland_2010_fields(maux, mbc, mx, my, xlower, ylower,    &
                                       dx, dy, t, aux, wind_index,           &
                                       pressure_index, storm)

        implicit none

        ! Time of the wind field requested
        integer, intent(in) :: maux,mbc,mx,my
        real(kind=8), intent(in) :: xlower,ylower,dx,dy,t

        ! Storm description, need in out here since we may update the storm
        ! if at next time point
        type(model_storm_type), intent(inout) :: storm

        ! Array storing wind and presure field
        integer, intent(in) :: wind_index, pressure_index
        real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

        print *, "This model has not yet been implemented!"
        stop

    end subroutine set_holland_2010_fields


    ! ==========================================================================
    !  Use the Chavas, Lin, and Emmanuel 2016 model
    ! ==========================================================================
    subroutine set_CLE_fields(maux, mbc, mx, my, xlower, ylower,    &
                              dx, dy, t, aux, wind_index,           &
                              pressure_index, storm)

        implicit none

        ! Time of the wind field requested
        integer, intent(in) :: maux,mbc,mx,my
        real(kind=8), intent(in) :: xlower,ylower,dx,dy,t

        ! Storm description, need in out here since we may update the storm
        ! if at next time point
        type(model_storm_type), intent(inout) :: storm

        ! Array storing wind and presure field
        integer, intent(in) :: wind_index, pressure_index
        real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

        print *, "This model has not yet been implemented!"
        stop

    end subroutine set_CLE_fields

end module model_storm_module
