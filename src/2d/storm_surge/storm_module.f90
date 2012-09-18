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
    type(constant_storm_type), pointer :: constant_storm => null()

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

        use geoclaw_module, only: pi

        implicit none
        
        ! Input arguments
        character(len=*), optional, intent(in) :: data_file
        
        ! Locals
        integer, parameter :: unit = 13
        integer :: status, i
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
        read(unit,*) storm_type
        read(unit,*) storm_file_path

        ! Read in hurricane track data from file
        if (storm_type == 0) then
            ! No storm will be used
        else if (storm_type == 1) then
            ! Track file with Holland reconstruction
            call set_holland_storm(storm_file_path,holland_storm)
            ! Set rho_air and ambient pressure in storms
        else if (storm_type == 2) then
            ! constant track and holland reconstruction
            call set_constant_storm(storm_file_path,constant_storm)
            ! Set rho_air and ambient pressure in storms
        else if (storm_type == 3) then
            ! Stommel wind field
            stop "Call stommel storm setup routine."
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

        use holland_storm_module, only: holland_storm_location
        use constant_storm_module, only: constant_storm_location

        implicit none

        ! Input
        real(kind=8), intent(in) :: t

        ! Output
        real(kind=8) :: location(2)        

        select case(storm_type)
            case(0)
                location = [0.d0, 0.d0]
            case(1)
                location = holland_storm_location(t,holland_storm)
            case(2)
                location = constant_storm_location(t,constant_storm)
            case(3)
                location = [0.d0, 0.d0]
        end select

    end function storm_location


    subroutine set_storm_fields(maxmx,maxmy,maux,mbc,mx,my,xlower,ylower,dx,dy,&
                                t,aux)

        implicit none

        ! Input arguments
        integer, intent(in) :: maxmx, maxmy, maux, mbc, mx, my
        real(kind=8), intent(in) :: xlower, ylower, dx, dy, t
        real(kind=8), intent(in out) :: aux(maux,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc)

    end subroutine set_storm_fields


    subroutine output_storm_location(t)

        implicit none

        real(kind=8), intent(in) :: t

        write(track_unit,"(3e26.16)") t,storm_location(t)

    end subroutine output_storm_location

end module storm_module



    ! ========================================================================
    !   subroutine hurricane_wind(mbc,mx,my,xlower,ylower,dx,dy,R_eye,wind)
    ! ========================================================================
    ! Calculates an constant 2d field of wind with the strength profile used
    ! from Weisberg and Zheng (2006).
    !
    ! Input:
    !     mbc = Number of ghost cells
    !     mx = Number of grid cells in x direction
    !     my = Number of grid cells in y direction
    !     xlower = Lower coordinate of computational grid in x
    !     ylower = Lower coordinate of computational grid in y
    !     dx = Grid spacing in x
    !     dy = Grid spacing in y
    !     t = Current time
    !
    ! Output:
    !     wind = Velocity of wind (m/s) (aux array 4 (x) and 5 (y))
    !
    ! Hurricane parameters:
    !     A,B = fit parameters for the hurricane
    !     rho_air = density of air
    !     Pn = Ambient atmospheric pressure
    !     Pc = Central pressure of hurricane
    ! ========================================================================
!     subroutine hurricane_wind(maxmx,maxmy,maux,mbc,mx,my,xlower,ylower,dx,dy,t,aux)

!         use geoclaw_module, only: icoriolis,coordinate_system,pi

!         implicit none

!         ! Input arguments
!         integer, intent(in) :: maxmx,maxmy,mbc,mx,my,maux
!         double precision, intent(in) :: xlower,ylower,dx,dy,t

!         ! Output
!         double precision, intent(inout) :: aux(maux,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc)
    
!         ! Local variables
!         integer :: i,j
!         double precision :: x,y,C,r,w,R_eye(2),f,L
        
!         ! If wind forcing is turned off, then set the aux array to zeros
!         if (.not.wind_forcing) then
!             aux(4:5,:,:) = 0.d0
!             return
!         endif
        
!         ! Initialize f in case we don't need it
!         f = 0.d0
        
!         if (wind_type == 0) then
!             stop "Unimplemented wind field type!"
!         else if (wind_type == 1) then
!             ! Hurrican eye location
!             R_eye = t * hurricane_velocity + R_eye_init
        
!             ! Parameter constant
!             C = 1d1**2 * A * B * (Pn-Pc) / rho_air
            
!             ! Set the wind    
!             do j=1-mbc,my+mbc
!                 y = ylower + (j-0.5d0) * dy - R_eye(2)
!                 ! Coriolis term
!                 if (icoriolis == 1) f = coriolis(y)
!                 do i=1-mbc,mx+mbc
!                     x = xlower + (i-0.5d0) * dx - R_eye(1)
!                     r = sqrt(x**2+y**2) * 1d-3
                
!                     if (abs(r) < 10d-6) then
!                         aux(4:5,i,j) = 0.d0
!                     else
!                         w = sqrt(C * exp(-A/r**B) / r**B + r**2 * f**2 / 4.0) &
!                                  - r * f / 2.d0
!                         r = r * 1d3
!                         aux(4,i,j) = -w * y / r
!                         aux(5,i,j) =  w * x / r
!                     endif
!                 enddo
!             enddo
!         ! Stommel wind field
!         else if (wind_type == 2) then
!             ! This corresponds to an effective wind stress of 0.2 in maximum
!             ! amplitude, the division by 1.2 is to account for the variable
!             ! wind speed coefficient
!             L = my*dy
!             do j=1-mbc,my+mbc
!                 y = ylower + (j-0.5d0) * dy
!                 aux(4:5,j,1) = -A * cos(pi*y/L)
!             enddo
!             aux(4:5,:,2) = 0.d0
!         endif
        
!         ! Ramp up
!         if (t < 0.d0) then
!             aux(4:5,:,:) = aux(4:5,:,:) * exp(-(t/(ramp_up_time*0.45d0))**2)
!         endif
        
!     end subroutine hurricane_wind

    ! ========================================================================
    !   subroutine hurricane_pressure(mbc,mx,my,xlower,ylower,dx,dy,R_eye,pressure)
    ! ========================================================================
    ! Calculates an constant 2d field of preesure with the strength profile 
    ! used from Weisberg and Zheng (2006).
    !
    ! Input:
    !     mbc = Number of ghost cells
    !     mx = Number of grid cells in x direction
    !     my = Number of grid cells in y direction
    !     xlower = Lower coordinate of computational grid in x
    !     ylower = Lower coordinate of computational grid in y
    !     dx = Grid spacing in x
    !     dy = Grid spacing in y
    !     R_eye_X = Location of the eye of the hurricane in x
    !     R_eye_Y = Location of the eye of the hurricane in y
    !
    ! Output:
    !     pressure = Atmospheric pressure (mb) (aux array 6)
    !
    ! Hurricane parameters:
    !     A,B = fit parameters for the hurricane
    !     rho_air = density of air
    !     Pn = Ambient atmospheric pressure
    !     Pc = Central pressure of hurricane
    ! ========================================================================
!     subroutine hurricane_pressure(maxmx,maxmy,maux,mbc,mx,my,xlower,ylower,dx,dy,t,aux)

!         implicit none

!         ! Input arguments
!         integer, intent(in) :: maxmx,maxmy,mbc,mx,my,maux
!         double precision, intent(in) :: xlower,ylower,dx,dy,t
    
!         ! Output
!         double precision, intent(inout) :: aux(maux,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc)

!         ! Local variables
!         integer :: i,j
!         double precision :: r,x,y,R_eye(2)
    
!         ! If pressure forcing is turned off, then set the aux array to Pn
!         if (.not.pressure_forcing) then
!             aux(6,:,:) = Pn
!             return
!         endif
    
!         ! Hurrican eye location
!         R_eye = t * hurricane_velocity + R_eye_init
    
!         do i=1-mbc,mx+mbc
!             do j=1-mbc,my+mbc
!                 x = xlower + (i-0.5d0) * dx - R_eye(1)
!                 y = ylower + (j-0.5d0) * dy - R_eye(2)
!                 r = sqrt(x**2+y**2)
                
!                 if (abs(r) < 10d-3) then
!                     aux(6,i,j) = Pn
!                 else
!                     aux(6,i,j) = Pc + (Pn-Pc) * exp(-1.d3**B*A/abs(r)**B)
!                 endif
!             enddo
!         enddo
        
!         ! Ramp up
!         if (t < 0.d0) then
!             aux(6,:,:) = Pn  + (aux(6,:,:) - Pn) * exp(-(t/(ramp_up_time*0.45d0))**2)
!         endif
        
!         ! Convert to Pa instead of millibars
!         aux(6,:,:)  = aux(6,:,:) * 100.d0 
!     end subroutine hurricane_pressure
    
