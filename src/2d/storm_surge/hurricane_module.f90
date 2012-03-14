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
    
    ! Algorithm parameters
    logical :: momentum_refinement
    integer :: max_speed_nest,max_R_nest
    double precision, allocatable :: speed_refine(:), R_refine(:)
    
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
    subroutine set_hurricane_params(data_file)

        use geoclaw_module, only: pi

        implicit none
        
        ! Input arguments
        character(len=*), optional, intent(in) :: data_file
        
        ! Locals
        integer :: ios,i
        character*150 :: hurricane_track_file_path
        
        ! Open file
        if (present(data_file)) then
            open(unit=13,file=data_file,iostat=ios,status="old",action="read",access="sequential")
        else
            open(unit=13,file='hurricane.data',iostat=ios,status="old",action="read",access="sequential")
        endif
        if ( ios /= 0 ) stop "Error opening file data_file"
        
        ! Read in parameters
        ! Forcing terms
        read(13,*) wind_forcing
        read(13,*) pressure_forcing
        read(13,*)
        
        ! Source term algorithm parameters
        read(13,*) wind_tolerance
        read(13,*) pressure_tolerance
        read(13,*)
        
        ! AMR parameters
        read(13,*) max_wind_nest
        allocate(wind_refine(max_wind_nest))
        read(13,*) (wind_refine(i),i=1,max_wind_nest)
        read(13,*) max_speed_nest
        allocate(speed_refine(max_speed_nest))
        read(13,*) (speed_refine(i),i=1,max_speed_nest)
        read(13,*) momentum_refinement
        read(13,*) max_R_nest
        allocate(R_refine(max_R_nest))
        read(13,*) (R_refine(i),i=1,max_R_nest)
        read(13,*)
        
        ! Physics
        read(13,*) rho_air
        read(13,*) theta_0
        theta_0 = theta_0 * pi / 180d0
        read(13,*)
        
        ! Wind 
        read(13,*) wind_type
        ! Read in hurricane track data from file
        if (wind_type == 0) then
            read(13,*) hurricane_track_file_path
            call read_hurricane_track_file(hurricane_track_file_path)
        ! Idealized hurricane
        else if (wind_type == 1) then
            read(13,*) ramp_up_time
            read(13,*) hurricane_velocity
            read(13,*) R_eye_init
            read(13,*) A
            read(13,*) B
            read(13,*) Pn
            read(13,*) Pc
        ! Stommel wind field
        else if (wind_type == 2) then
            read(13,*) A
        else
            print *,"Invalid wind type ",wind_type," provided."
            stop
        endif
        
        close(13)

    end subroutine set_hurricane_params
    
    !
    !
    subroutine read_hurricane_track_file(file_path)

        implicit none
    
        ! Input
        character*150, intent(in) :: file_path
    
        ! Local 
        integer :: i,ios,num_points
        integer, parameter :: log_unit = 42
        integer, parameter :: data_unit = 43
        character*50 :: line_format
    
        ! Line format
        line_format = "(5d16.8)"
    
        ! Open log file
        open(unit=log_unit, file='./fort.hurricane', iostat=ios, status="unknown", action="write")
        if ( ios /= 0 ) stop "Error opening file ./fort.hurricane"
    
        ! Open hurricane data file
        open(unit=data_unit, file=file_path, iostat=ios, status="old", action="read")
        if ( ios /= 0 ) then
            print *, "Error opening file ", file_path
            stop
        endif
    
        ! Read header
        read(data_unit,'(i3)') num_points
    
        ! Allocate arrays
        allocate(time(num_points),position(num_points,2))
        allocate(max_wind(num_points),min_pressure(num_points))
    
        do i=1,num_points
            read(data_unit,*) time(i),position(i,1),position(i,2),max_wind(i),min_pressure(i)
        enddo
    
    end subroutine read_hurricane_track_file

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

        use geoclaw_module, only: icoriolis,icoordsys,pi

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
        
        if (wind_type == 1) then
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
    !   If icoordsys == 1 then
    !       A beta-plane approximation is used and y should be in meters
    !   if icoordsys == 2 then
    !       Grid is in lat-long and y should be in degrees which is converted
    !       to radians
    ! ========================================================================
    double precision function coriolis(y)
    
        use geoclaw_module, only: icoordsys,Rearth,pi,icoriolis
    
        implicit none
        
        ! Input
        double precision, intent(in) :: y
        
        ! Locals
        double precision :: theta
        
        ! Angular speed of earth = 2.d0*pi/(86400.d0) 
        double precision, parameter :: OMEGA = 7.2722052166430395d-05
        
        ! Assume beta plane approximation and y is in meters
        if (icoriolis == 1) then
            if (icoordsys == 1) then
                theta = y / 111d3 * pi / 180d0 + theta_0
                coriolis = 2.d0 * OMEGA * (sin(theta_0) + (theta - theta_0)     &
                                                        * cos(theta_0))
            else if (icoordsys == 2) then        
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