! ==============================================================================
! storm_module 
!
! Module contains routines for constructing a wind and pressure field based on
! the Holland hurricane module.  
! 
! Many of these routines have been adapted from PADCIRC version 45.12 03/17/2006
! ==============================================================================
module holland_storm_module

    implicit none
    save

    ! Holland storm type definition
    type holland_storm_type
        ! Fore/hindcast size and current position
        integer :: index,num_casts

        ! These parameters are located at time points but are interpolated in
        ! time and space when the relevant fields are requested.

        ! Location and speed of storm
        real(kind=8), allocatable :: track(:,:)
        real(kind=8) :: velocity(2)

        ! Storm physics
        real(kind=8) :: ambient_pressure = 101.3d3 ! Pascals
        real(kind=8) :: rho_air = 1.3d0
        real(kind=8), allocatable :: max_wind_radius(:)
        real(kind=8), allocatable :: max_wind_speed(:)
        real(kind=8), allocatable :: central_pressure(:)

    end type holland_storm_type

    ! Interal tracking variables for storm
    real(kind=8), private :: A,B

contains

    ! Setup routine for the holland model
    subroutine set_holland_storm(storm_data_path,storm)

        implicit none

        ! Subroutine I/O
        character(len=*), optional :: storm_data_path
        type(holland_storm_type), intent(in out) :: storm

        ! Local storage
        integer, parameter :: data_file = 701
        integer :: i,lines,io_status,num_casts
        real(kind=8) :: forecast_time,last_time

        ! Reading buffer variables
        integer :: year,month,day,hour,forecast,lat,lon,max_wind_speed
        integer :: central_pressure,RRP,max_wind_radius
        character(len=4) :: cast_type

        ! File format string
        character(len=*), parameter :: file_format = "(8x,i4,i2,i2,i2,6x,a4,2x,i3,1x,i4,3x,i4,3x,i3,2x,i4,47x,i3,2x,i3)"

        ! Open data file
        if (present(storm_data_path)) then
            open(unit=data_file,file=storm_data_path,status='old', &
                 action='read',iostat=io_status)
        else
            open(unit=data_file,file="./hurricane.data",status='old', &
                 action='read',iostat=io_status)
        endif
        if (io_status /= 0) then
            print "(a,i1)", "Error opening storm data file. status = ", io_status
            stop 
        endif

        ! Count number of data lines
        num_casts = 0
        do
            read(data_file,fmt=file_format,iostat=io_status) year,month,day, &
                    hour,cast_type,forecast,lat,lon,max_wind_speed, &
                    central_pressure,RRP,max_wind_radius

            ! Exit loop if we ran into an error or we reached the end of the file
            if (io_status /= 0) exit

            ! Skip counting this line if time is repeated
            forecast_time = date_to_seconds(year,month,day,hour,0,0.d0)
            if (abs(forecast_time - last_time) >= 1.8d3) then
                num_casts = num_casts + 1
            endif
            last_time = forecast_time
        end do
        rewind(data_file)

        ! Allocate storm parameter file variables
        allocate(storm%track(3,num_casts))
        allocate(storm%max_wind_speed(num_casts))
        allocate(storm%max_wind_radius(num_casts))
        allocate(storm%central_pressure(num_casts))

        ! Now re-read the file's contents
        do i=1,num_casts
            read(data_file,fmt=file_format) year, month, day, hour, cast_type, &
                    forecast, lat, lon, max_wind_speed, central_pressure, RRP, &
                    max_wind_radius

            ! Storm position
            ! Conversions:
            !  lat - Convert 10ths of degs to degs
            !  lon - Convert 10ths of degs to degs, also assume W hemisphere
            storm%track(1,i) = date_to_seconds(year,month,day,hour,0,0.d0)
            storm%track(2:3,i) = [real(lat,kind=8) / 10.d0, &
                                 -real(lon,kind=8) / 10.d0]

            ! Storm intensity
            ! Conversions:
            !  max_wind_speed - Convert knots to m/s
            !  radius_max_wind  - convert from nm to m
            !  central_pressure - convert from mbar to Pa
            storm%max_wind_speed(i) = real(max_wind_speed,kind=8) * 0.51444444d0
            storm%max_wind_radius(i) = real(max_wind_radius,kind=8)  * 1.852000003180799d3
            storm%central_pressure(i) = real(central_pressure,kind=8) * 100.d0

        enddo

        ! Set data index to first valid value for time
        storm%index = 1

        ! Initialize stored_velocity
        storm%velocity = 0.d0

    end subroutine set_holland_storm

    ! ==========================================================================
    !  real(kind=8) pure date_to_seconds(year,months,days,hours,minutes,seconds)
    !    Convert time from year, month, day, hour, min, sec to seconds since the
    !    beginning of the year.
    ! ==========================================================================
    pure real(kind=8) function date_to_seconds(year,months,days,hours,minutes, &
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

        use geoclaw_module, only: earth_radius, pi

        implicit none
        
        ! Input
        real(kind=8), intent(in) :: t
        type(holland_storm_type), intent(in out) :: storm

        ! Atmospheric boundary layer, input variable in ADCIRC but always is
        ! set to the following value
        real(kind=8), parameter :: atmos_boundary_layer = 0.9d0

        ! Local storage
        real(kind=8) :: dp, modified_max_wind_speed, trans_speed
        real(kind=8) :: x(2),y(2),ds,dt,dx,dy

        ! Update storm variable

        ! Check if we are past end of fore/hind cast
        ! If this is the case, we do not need to do anything here since the 
        ! Holland parameters will no longer change
        if (storm%index < storm%num_casts) then
            if (t >= storm%track(1,storm%index+1)) then
                ! Update index of next storm
                storm%index = storm%index + 1

                ! Calculate stored storm parameters

                ! Subtract translational speed of storm from maximum wind speed
                ! to avoid distortion in the Holland curve fit.  Added back later
                trans_speed = sqrt(storm%velocity(1)**2 + storm%velocity(2)**2)
                modified_max_wind_speed = storm%max_wind_speed(storm%index) &
                                          - trans_speed

                ! Convert wind speed (10 m) to top of atmospheric boundary layer
                modified_max_wind_speed = modified_max_wind_speed / atmos_boundary_layer
                
                ! Calculate central pressure difference
                dp = storm%ambient_pressure - storm%central_pressure(storm%index)
                ! Limit central pressure deficit due to bad ambient pressure,
                ! really should have better ambient pressure...
!                 if (dp < 100.d0) then
!                     dp = 100.d0
!                 endif

                ! Calculate Holland parameters and limit the result
                B = storm%rho_air * exp(1.d0) * (modified_max_wind_speed**2) / dp
                if (B <  1.d0) B = 1.d0
                if (B > 2.5d0) B = 2.5d0

                ! Calculate Holland A parameter, not needed
                A = (storm%max_wind_speed(storm%index)*1.d3)**B

                ! Update storm velocity
                ! Calculate velocity based on great circle distance between
                ! locations of storm
                x = storm%track(2:3,storm%index)
                y = storm%track(2:3,storm%index+1)

                dt = storm%track(1,storm%index+1) - storm%track(1,storm%index)
                
                dx = (y(1) - x(1)) * pi / 180.d0
                ds = earth_radius * (2.d0 * asin(sqrt(cos(y(2) * pi / 180.d0) &
                                                    * cos(x(2) * pi / 180.d0) &
                                                    * sin(0.5d0 * dx)**2) ) )
                storm%velocity(1) = sign(ds/dt,dx)
                
                dy = (y(2) - y(1)) * pi / 180.d0
                ds = earth_radius * (2.d0 * asin(sqrt(sin( 0.5d0 * dy)**2.0d0))) 
                storm%velocity(2) = sign(ds/dt,dy)

            endif
        endif

    end subroutine update_storm

    ! ==========================================================================
    !  holland_storm_location(t,storm)
    !    Interpolate location of hurricane in the current time interval
    ! ==========================================================================
    function holland_storm_location(t,storm) result(location)

        implicit none

        ! Input
        real(kind=8), intent(in) :: t
        type(holland_storm_type), intent(in out) :: storm

        ! Output
        real(kind=8) :: location(2)

        ! Locals
        real(kind=8) :: weight,t0,t1

        ! Update storm
        call update_storm(t,storm)

        ! Interpolate location
        if (storm%index == storm%num_casts) then
            ! Past end forecast, find position using stored velocity and last
            ! direction
            stop "Location past forecast interval has not been implemented!"
        else
            ! Interpolate in time using lat,long coordinates
            t0 = storm%track(1,storm%index)
            t1 = storm%track(1,storm%index + 1)
            weight = (t - t0) / (t1 - t0)
            location(1) = storm%track(2,storm%index) + weight * &
                (storm%track(2,storm%index + 1) - storm%track(2,storm%index))
            location(2) = storm%track(3,storm%index) + weight * &
                (storm%track(3,storm%index + 1) - storm%track(3,storm%index))
        end if

    end function holland_storm_location

    ! ==========================================================================
    !  set_storm_fields()
    ! ==========================================================================
    subroutine set_storm_fields(maxmx,maxmy,maux,mbc,mx,my,xlower,ylower,dx,dy,&
                                    t,aux, wind_index, pressure_index, storm)

        use geoclaw_module, only: pi, omega, g => grav, rho, num_layers

        implicit none

        ! Time of the wind field requested
        integer, intent(in) :: maxmx,maxmy,maux,mbc,mx,my
        real(kind=8), intent(in) :: xlower,ylower,dx,dy,t

        ! Storm description, need in out here since we may update the storm
        ! if at next time point
        type(holland_storm_type), intent(in out) :: storm

        ! Array storing wind and presure field
        integer, intent(in) :: wind_index, pressure_index
        real(kind=8), intent(inout) :: aux(maux,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc)

        ! Local storage
        real(kind=8) :: x,y,r,sloc(2)
        real(kind=8) :: coriolis
        integer :: i,j

        ! First update storm if necessary
        call update_storm(t,storm)

        ! Calculate storm's current location
        sloc = holland_storm_location(t,storm)

        ! Set wind and pressure field

        aux(wind_index:wind_index+1,:,:) = 0.d0
        aux(pressure_index,:,:) = storm%ambient_pressure
        
    end subroutine set_storm_fields

end module holland_storm_module