! ============================================================================
!  Program:     /Users/mandli/src/local/shallow_water/2d/storm_surge/testbeach
!  File:        hurricane
!  Created:     2009-11-11
!  Author:      Kyle Mandli
! ============================================================================
!      Copyright (C) 2009-11-11 Kyle Mandli <mandli@amath.washington.edu>
!
!  Distributed under the terms of the Berkeley Software Distribution (BSD) 
!  license
!                     http://www.opensource.org/licenses/
! ============================================================================

module hurricane_module

    implicit none
    
    ! Wind source terms
    logical :: wind_forcing
    integer :: wind_type,max_wind_nest
    double precision, allocatable :: wind_refine(:)
    double precision :: wind_tolerance
    
    ! Pressure source terms
    logical :: pressure_forcing
    double precision :: pressure_tolerance
    
    ! Refinement parameters
    integer :: max_R_nest
    double precision, allocatable :: R_refine(:)
    
    ! Hurricane Parameters
    double precision :: ramp_up_time,hurricane_velocity(2),R_eye_init(2)
    double precision :: rho_air
    double precision, private :: A,B,Pn,Pc
    
    double precision, allocatable :: time(:),position(:,:)
    double precision, allocatable :: max_wind(:),min_pressure(:)
    
    ! Coriolis location for beta-plane approximation
    !  should be read in but not implemented that way yet
    double precision :: theta_0

contains
    ! ========================================================================
    !   subroutine set_hurricane_parameters(data_file)
    ! ========================================================================
    ! Reads in the data file at the path data_file.  Sets the following 
    ! parameters:
    !     A,B = fit parameters for the hurricane
    !     rho_air = density of air
    !     Pn = Ambient atmospheric pressure
    !     Pc = Central pressure of hurricane
    !
    ! Input:
    !     data_file = Path to data file
    !
    ! ========================================================================
    subroutine set_hurricane(data_file)

        use geoclaw_module, only: pi

        implicit none
        
        ! Input arguments
        character(len=*), optional, intent(in) :: data_file
        
        ! Locals
        integer, parameter :: unit = 13
        integer :: status,i
        character*150 :: hurricane_track_file_path
        
        ! Open file
        if (present(data_file)) then
            call opendatafile(unit,data_file)
        else
            call opendatafile(unit,'hurricane.data')
        endif
        
        ! Read in parameters
        ! Forcing terms
        read(unit,*) wind_forcing
        read(unit,*) pressure_forcing
        read(unit,*)
        
        ! Source term algorithm parameters
        read(unit,*) wind_tolerance
        read(unit,*) pressure_tolerance
        read(unit,*)
        
        ! AMR parameters
        read(unit,*) max_wind_nest
        allocate(wind_refine(max_wind_nest))
        read(unit,*) (wind_refine(i),i=1,max_wind_nest)
        read(unit,*) max_R_nest
        allocate(R_refine(max_R_nest))
        read(unit,*) (R_refine(i),i=1,max_R_nest)
        read(unit,*)
        
        ! Physics
        read(unit,*) rho_air
        read(unit,*) theta_0
        theta_0 = theta_0 * pi / 180d0
        read(13,*)
        
        ! Wind 
        read(unit,*) wind_type
        ! Read in hurricane track data from file
        if (wind_type == 0) then
            read(unit,*) hurricane_track_file_path
            call read_hurricane_track_file(hurricane_track_file_path)
        ! Idealized hurricane
        else if (wind_type == 1) then
            read(unit,*) ramp_up_time
            read(unit,*) hurricane_velocity
            read(unit,*) R_eye_init
            read(unit,*) A
            read(unit,*) B
            read(unit,*) Pn
            read(unit,*) Pc
        ! Stommel wind field
        else if (wind_type == 2) then
            read(13,*) A
        else
            print *,"Invalid wind type ",wind_type," provided."
            stop
        endif
        
        close(unit)

    end subroutine set_hurricane

    subroutine read_hurricane_track_file(track_file)
        
        implicit none
        
        ! Subroutine Arguments
        character(len=*), intent(in) :: track_file
        
        ! Locals
        integer, parameter :: unit = 14
        
    end subroutine read_hurricane_track_file
    
    ! ========================================================================
    !   function eye_location(t)
    ! ========================================================================
    ! Returns the location of the eye of the storm
    ! 
    ! Input:
    !     real(kind=8) t = current time
    ! 
    ! Output:
    !     real(kind=8) eye_location(2) = (x,y) location of the eye
    ! ========================================================================
    pure function eye_location(t)
        implicit none
        
        real(kind=8), intent(in) :: t
        real(kind=8) :: eye_location(2)

        eye_location = t * hurricane_velocity + R_eye_init
        
    end function eye_location
    
    ! ========================================================================
    !   subroutine hurricane_wind(mbc,mx,my,xlower,ylower,dx,dy,R_eye,wind)
    ! ========================================================================
    ! Calculates an idealized 2d field of wind with the strength profile used
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
    subroutine hurricane_wind(maxmx,maxmy,maux,mbc,mx,my,xlower,ylower,dx,dy,t,aux)

        use geoclaw_module, only: icoriolis,coordinate_system,pi

        implicit none

        ! Input arguments
        integer, intent(in) :: maxmx,maxmy,mbc,mx,my,maux
        double precision, intent(in) :: xlower,ylower,dx,dy,t

        ! Output
        double precision, intent(inout) :: aux(maux,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc)
    
        ! Local variables
        integer :: i,j
        double precision :: x,y,C,r,w,R_eye(2),f,L
        
        ! If wind forcing is turned off, then set the aux array to zeros
        if (.not.wind_forcing) then
            aux(4:5,:,:) = 0.d0
            return
        endif
        
        ! Initialize f in case we don't need it
        f = 0.d0
        
        if (wind_type == 0) then
            stop "Unimplemented wind field type!"
        else if (wind_type == 1) then
            ! Hurrican eye location
            R_eye = t * hurricane_velocity + R_eye_init
        
            ! Parameter constant
            C = 1d1**2 * A * B * (Pn-Pc) / rho_air
            
            ! Set the wind    
            do j=1-mbc,my+mbc
                y = ylower + (j-0.5d0) * dy - R_eye(2)
                ! Coriolis term
                if (icoriolis == 1) f = coriolis(y)
                do i=1-mbc,mx+mbc
                    x = xlower + (i-0.5d0) * dx - R_eye(1)
                    r = sqrt(x**2+y**2) * 1d-3
                
                    if (abs(r) < 10d-6) then
                        aux(4:5,i,j) = 0.d0
                    else
                        w = sqrt(C * exp(-A/r**B) / r**B + r**2 * f**2 / 4.0) &
                                 - r * f / 2.d0
                        r = r * 1d3
                        aux(4,i,j) = -w * y / r
                        aux(5,i,j) =  w * x / r
                    endif
                enddo
            enddo
        ! Stommel wind field
        else if (wind_type == 2) then
            ! This corresponds to an effective wind stress of 0.2 in maximum
            ! amplitude, the division by 1.2 is to account for the variable
            ! wind speed coefficient
            L = my*dy
            do j=1-mbc,my+mbc
                y = ylower + (j-0.5d0) * dy
                aux(4:5,j,1) = -A * cos(pi*y/L)
            enddo
            aux(4:5,:,2) = 0.d0
        endif
        
        ! Ramp up
        if (t < 0.d0) then
            aux(4:5,:,:) = aux(4:5,:,:) * exp(-(t/(ramp_up_time*0.45d0))**2)
        endif
        
    end subroutine hurricane_wind

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
    ! ========================================================================
    double precision function wind_drag(wind_speed)
    
        implicit none
        
        ! Input
        double precision, intent(in) :: wind_speed
        
        
        if (wind_speed <= 11.d0) then
            wind_drag = 1.2d0
        else if ((wind_speed > 11.d0).and.(wind_speed <= 25.d0)) then
            wind_drag = 0.49d0 + 0.065d0 * wind_speed
        else
            wind_drag = 0.49 + 0.065d0 * 25.d0
        endif
        
        wind_drag = wind_drag * 10.d-3
    
    end function wind_drag

    ! ========================================================================
    !   subroutine hurricane_pressure(mbc,mx,my,xlower,ylower,dx,dy,R_eye,pressure)
    ! ========================================================================
    ! Calculates an idealized 2d field of preesure with the strength profile 
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
    subroutine hurricane_pressure(maxmx,maxmy,maux,mbc,mx,my,xlower,ylower,dx,dy,t,aux)

        implicit none

        ! Input arguments
        integer, intent(in) :: maxmx,maxmy,mbc,mx,my,maux
        double precision, intent(in) :: xlower,ylower,dx,dy,t
    
        ! Output
        double precision, intent(inout) :: aux(maux,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc)

        ! Local variables
        integer :: i,j
        double precision :: r,x,y,R_eye(2)
    
        ! If pressure forcing is turned off, then set the aux array to Pn
        if (.not.pressure_forcing) then
            aux(6,:,:) = Pn
            return
        endif
    
        ! Hurrican eye location
        R_eye = t * hurricane_velocity + R_eye_init
    
        do i=1-mbc,mx+mbc
            do j=1-mbc,my+mbc
                x = xlower + (i-0.5d0) * dx - R_eye(1)
                y = ylower + (j-0.5d0) * dy - R_eye(2)
                r = sqrt(x**2+y**2)
                
                if (abs(r) < 10d-3) then
                    aux(6,i,j) = Pn
                else
                    aux(6,i,j) = Pc + (Pn-Pc) * exp(-1.d3**B*A/abs(r)**B)
                endif
            enddo
        enddo
        
        ! Ramp up
        if (t < 0.d0) then
            aux(6,:,:) = Pn  + (aux(6,:,:) - Pn) * exp(-(t/(ramp_up_time*0.45d0))**2)
        endif
        
        ! Convert to Pa instead of millibars
        aux(6,:,:)  = aux(6,:,:) * 100.d0 
    end subroutine hurricane_pressure
    
    ! ========================================================================
    !  Calculate the coriolis constant f
    !   If coordinate_system == 1 then
    !       A beta-plane approximation is used and y should be in meters
    !   if coordinate_system == 2 then
    !       Grid is in lat-long and y should be in degrees which is converted
    !       to radians
    ! ========================================================================
    double precision function coriolis(y)
    
        use geoclaw_module, only: coordinate_system,earth_radius,pi,icoriolis
    
        implicit none
        
        ! Input
        double precision, intent(in) :: y
        
        ! Locals
        double precision :: theta
        
        ! Angular speed of earth = 2.d0*pi/(86400.d0) 
        double precision, parameter :: OMEGA = 7.2722052166430395d-05
        
        ! Assume beta plane approximation and y is in meters
        if (icoriolis == 1) then
            if (coordinate_system == 1) then
                theta = y / 111d3 * pi / 180d0 + theta_0
                coriolis = 2.d0 * OMEGA * (sin(theta_0) + (theta - theta_0)     &
                                                        * cos(theta_0))
            else if (coordinate_system == 2) then        
                theta = pi*y/180.d0
                coriolis = 2.d0 * OMEGA * sin(theta)
            else
                print *,"Unknown coordinate system, unable to calculate coriolis."
                coriolis = 0.d0
            endif
        ! Stommel problem coriolis
        else if (icoriolis == 2) then
            ! Seems to correlate with a theta_0 = 43.5 degrees latitude
            coriolis = 1d-4 + 1d-11 * (y - 500.0d3)
        else
            coriolis = 0.d0
        endif
    end function coriolis
    
end module hurricane_module