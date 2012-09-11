! ==============================================================================
! idealized_storm_module 
!
! Module contains routines for constructing a wind and pressure field based on
! a contant Holland based storm.
! ==============================================================================
module idealized_storm_module

    implicit none
    save

    ! Holland storm type definition
    type ideal_storm_type

        real(kind=8) :: ramp_up_time
        real(kind=8) :: velocity(2)
        real(kind=8) :: R_eye_init(2)
        real(kind=8) :: A
        real(kind=8) :: B
        real(kind=8) :: central_pressure
        real(kind=8) :: ambient_pressure = 101.3d3
        real(kind=8) :: rho_air = 1.3d0

    end type ideal_storm_type

contains

    ! Setup routine for the holland model
    type(ideal_storm_type) function set_ideal_storm(storm_data_path) result(storm)

        implicit none

        ! Path to data file
        character(len=*), optional :: storm_data_path

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
        if (io_status /= 0) stop "Error opening hurricane data file."

        ! Read in hurricane parameters
        read(13,*) storm%ramp_up_time
        read(13,*) storm%velocity
        read(13,*) storm%R_eye_init
        read(13,*) storm%A
        read(13,*) storm%B
        read(13,*) storm%central_pressure

    end function set_ideal_storm

    ! ==========================================================================
    !  holland_update_storm(t)
    !    Update storm information to time t
    ! ==========================================================================
    subroutine update_ideal_storm(t,storm)

        implicit none
        
        ! Input
        real(kind=8), intent(in) :: t
        type(ideal_storm_type), intent(in out) :: storm

        ! Nothing needs to be updated here

    end subroutine update_ideal_storm

    ! ==========================================================================
    !  storm_location(t,storm)
    !    Interpolate location of hurricane in the current time interval
    ! ==========================================================================
    function ideal_storm_location(t,storm) result(location)

        implicit none

        ! Input
        real(kind=8), intent(in) :: t
        type(ideal_storm_type), intent(in out) :: storm

        ! Output
        real(kind=8) :: location(2)

        location = t * storm%velocity + storm%R_eye_init

    end function ideal_storm_location

    ! ==========================================================================
    !  set_storm_fields()
    ! ==========================================================================
    subroutine set_ideal_storm_fields(maxmx,maxmy,maux,mbc,mx,my,xlower,ylower,dx,dy,&
                                    t,aux, wind_index, pressure_index, storm)

        use geoclaw_module, only: coriolis, coriolis_forcing

        implicit none

        ! Time of the wind field requested
        integer, intent(in) :: maxmx,maxmy,maux,mbc,mx,my
        real(kind=8), intent(in) :: xlower,ylower,dx,dy,t

        ! Storm description, need in out here since we may update the storm
        ! if at next time point
        type(ideal_storm_type), intent(in out) :: storm

        ! Array storing wind and presure field
        integer, intent(in) :: wind_index, pressure_index
        real(kind=8), intent(inout) :: aux(maux,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc)

        ! Locals
        integer :: i, j
        real(kind=8) :: x, y, r, R_eye(2), f, C, w

        ! Hurrican eye location
        R_eye = ideal_storm_location(t,storm)
    
        ! Parameter constant
        C = 1d1**2 * storm%A * storm%B * (storm%ambient_pressure - storm%central_pressure) / storm%rho_air
        
        do j=1-mbc,my+mbc
            y = ylower + (j-0.5d0) * dy - R_eye(2)
            ! Coriolis term
            if (coriolis_forcing) f = coriolis(y)
            do i=1-mbc,mx+mbc
                x = xlower + (i-0.5d0) * dx - R_eye(1)
                r = sqrt(x**2+y**2) * 1d-3
            
                ! Set wind
                if (abs(r) < 10d-6) then
                    aux(wind_index:wind_index + 1,i,j) = 0.d0
                else
                    w = sqrt(C * exp(-storm%A/r**storm%B) / r**storm%B + r**2 * f**2 / 4.0) &
                             - r * f / 2.d0
                    r = r * 1d3
                    aux(wind_index,i,j) = -w * y / r
                    aux(wind_index+1,i,j) =  w * x / r
                endif

                ! Set pressure
                if (abs(r) < 10d-3) then
                    aux(pressure_index,i,j) = storm%ambient_pressure
                else
                    aux(pressure_index,i,j) = storm%central_pressure + (storm%ambient_pressure-storm%central_pressure) &
                                    * exp(-1.d3**storm%B * storm%A / abs(r)**storm%B)
                endif
            enddo
        enddo

        ! Ramp up function
        if (t < 0.d0) then
            aux(wind_index  ,:,:) = aux(wind_index  ,:,:) &
                                * exp(-(t/(storm%ramp_up_time*0.45d0))**2)
            aux(wind_index+1,:,:) = aux(wind_index+1,:,:) &
                                * exp(-(t/(storm%ramp_up_time*0.45d0))**2)
            aux(pressure_index,:,:) = aux(pressure_index,:,:) &
                                * exp(-(t/(storm%ramp_up_time*0.45d0))**2)
        endif
        
    end subroutine set_ideal_storm_fields

end module idealized_storm_module