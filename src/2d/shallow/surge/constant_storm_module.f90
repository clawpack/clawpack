! ==============================================================================
! constantized_storm_module 
!
! Module contains routines for constructing a wind and pressure field based on
! a contant Holland based storm.
! ==============================================================================
module constant_storm_module

    implicit none
    save

    ! Holland storm type definition
    type constant_storm_type

        real(kind=8) :: ramp_up_time
        real(kind=8) :: velocity(2)
        real(kind=8) :: R_eye_init(2)
        real(kind=8) :: A
        real(kind=8) :: B
        real(kind=8) :: central_pressure
        real(kind=8) :: ambient_pressure = 101.3d3
        real(kind=8) :: rho_air = 1.3d0

    end type constant_storm_type

    logical, private :: module_setup = .false.
    logical, private, parameter :: DEBUG = .true.

contains

    ! Setup routine for the holland model
    subroutine set_constant_storm(storm_data_path, storm, log_unit)

        use geoclaw_module, only: coordinate_system

        implicit none

        ! Subroutine I/O
        character(len=*), intent(in), optional :: storm_data_path
        type(constant_storm_type), intent(in out) :: storm
        integer, intent(in) :: log_unit

        ! Local storage
        integer, parameter :: unit = 701

        if (.not.module_setup) then
            
            ! Storm type does not work for lat-long coordinate systems
            if (coordinate_system == 2) then
                stop "Constant storm type does not work on lat-long coordinates."
            endif

            ! Open data file
            if (present(storm_data_path)) then
                call opendatafile(unit,storm_data_path)
            else
                call opendatafile(unit,'storm.data')
            endif

            ! Read in hurricane parameters
            read(unit,*) storm%ramp_up_time
            read(unit,*) storm%velocity
            read(unit,*) storm%R_eye_init
            read(unit,*) storm%A
            read(unit,*) storm%B
            read(unit,*) storm%central_pressure

            ! Output to log file
            write(log_unit,*) "Storm Data - Constant Storm"
            write(log_unit,"('Ramp Up Time = ',d16.8)") storm%ramp_up_time
            write(log_unit,"('Velocity = ',2d16.8)") storm%velocity
            write(log_unit,"('Eye initial position = ',2d16.8)") storm%R_eye_init
            write(log_unit,"('Holland parameters (A,B) =',2d16.8)") storm%A, storm%B
            write(log_unit,"('Pressures (Pn,Pc) = ',2d16.8)") storm%ambient_pressure,&
                                                            storm%central_pressure
            write(log_unit,"('Density of Air = ',d16.8)") storm%rho_air

            module_setup = .true.
        end if

    end subroutine set_constant_storm


    ! ==========================================================================
    !  storm_location(t,storm)
    !    Interpolate location of hurricane in the current time interval
    ! ==========================================================================
    function constant_storm_location(t,storm) result(location)

        implicit none

        ! Input
        real(kind=8), intent(in) :: t
        type(constant_storm_type), intent(in out) :: storm

        ! Output
        real(kind=8) :: location(2)

        location = t * storm%velocity + storm%R_eye_init

    end function constant_storm_location

    ! ==========================================================================
    !  storm_direction(t,storm)
    !    Interpolate location of hurricane in the current time interval
    ! ==========================================================================
    real(kind=8) function constant_storm_direction(t,storm) result(theta)

        implicit none

        ! Input
        real(kind=8), intent(in) :: t
        type(constant_storm_type), intent(in out) :: storm

        ! Locals
        real(kind=8) :: junk(2)
        theta = atan2(storm%velocity(2),storm%velocity(1))

    end function constant_storm_direction

    ! ==========================================================================
    !  set_storm_fields()
    ! ==========================================================================
    subroutine set_constant_storm_fields(maux,mbc,mx,my,xlower,ylower,dx,dy,&
                                    t,aux, wind_index, pressure_index, storm)

        use geoclaw_module, only: coriolis, coriolis_forcing

        implicit none

        ! Time of the wind field requested
        integer, intent(in) :: maux,mbc,mx,my
        real(kind=8), intent(in) :: xlower,ylower,dx,dy,t

        ! Storm description, need in out here since we may update the storm
        ! if at next time point
        type(constant_storm_type), intent(in out) :: storm

        ! Array storing wind and presure field
        integer, intent(in) :: wind_index, pressure_index
        real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

        ! Locals
        integer :: i, j
        real(kind=8) :: x, y, r, r_km, R_eye(2), f, C, w

        ! Hurrican eye location
        R_eye = constant_storm_location(t,storm)
    
        ! Parameter constant
        C = storm%A * storm%B * (storm%ambient_pressure - storm%central_pressure) / storm%rho_air
        
        f = 0.d0
        do j=1-mbc,my+mbc
            y = ylower + (j-0.5d0) * dy - R_eye(2)
            
            ! Coriolis term
            if (coriolis_forcing) f = coriolis(y)

            do i=1-mbc,mx+mbc
                x = xlower + (i-0.5d0) * dx - R_eye(1)
                r = sqrt(x**2+y**2)
            
                ! Set wind
                if (abs(r) < 1.0d-4) then
                    aux(wind_index:wind_index + 1,i,j) = 0.d0
                else
                    w = sqrt(C * exp(-storm%A / (1d-3 * r)**storm%B) / (r * 1d-3)**storm%B &
                                + (1d-3 * r * f)**2 / 4.d0) - 1d-3 * r * f / 2.d0

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
            aux(pressure_index,:,:) = storm%ambient_pressure &
                        + (aux(pressure_index,:,:) - storm%ambient_pressure) &
                        * exp(-(t/(storm%ramp_up_time*0.45d0))**2)
        endif
        
    end subroutine set_constant_storm_fields

end module constant_storm_module