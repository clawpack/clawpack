module parameters_module

    implicit none
    
    ! Physics constants
    double precision, parameter :: PI = 3.141592654d0, g = 9.8d0
    double precision :: rho_air,rho(2)
    double precision :: r, one_minus_r,manning
    
    ! Algorithm settings
    double precision :: dry_tolerance
    integer :: eigen_method,inundation_method
    logical :: entropy_fix

    ! Initial condition settings
    integer :: init_type,wave_family
    double precision :: init_location,eta_left(2),eta_right(2)
    double precision :: u_left(2),u_right(2),epsilon,sigma,bathy_type
    double precision :: bathy_location,bathy_left,bathy_right
    double precision :: x0,basin_depth,x1,shelf_depth,shelf_slope
    
    ! Wind forcing settings
    integer, private :: wind_type
    double precision, private :: A,B,Pn,Pc,ramp_up_time,R_eye_init
    double precision, private :: hurricane_velocity,omega,N,t_length

contains

    ! ========================================================================
    !  Setup all parameters in the module
    ! ========================================================================
    subroutine setup_parameters(data_file)

        implicit none

        ! Input
        character(len=*), optional, intent(in) :: data_file

        ! Locals
        integer :: ios,i

        ! Open file
        if (present(data_file)) then
            call opendatafile(13,data_file)
        else
            call opendatafile(13,'problem.data')
        endif

        ! Physics constants
        read(13,*) rho_air
        read(13,*) rho
        r = rho(1) / rho(2)
        one_minus_r = 1.d0 - r
        read(13,*) manning
        read(13,*)
        
        ! Algorithm settings
        read(13,*) dry_tolerance
        read(13,*) eigen_method
        read(13,*) inundation_method
        read(13,*) entropy_fix
        read(13,*)
        
        ! Initial conditions
        read(13,*) init_type
        read(13,*) init_location
        read(13,*) wave_family
        read(13,*) eta_left
        read(13,*) eta_right
        read(13,*) u_left
        read(13,*) u_right
        read(13,*) epsilon
        read(13,*) sigma
        read(13,*)
        
        ! Bathymetry
        read(13,*) bathy_type
        if (bathy_type == 1) then
            read(13,*) bathy_location
            read(13,*) bathy_left
            read(13,*) bathy_right
        else if (bathy_type == 2) then
            read(13,*) x0
            read(13,*) basin_depth
            read(13,*) x1
            read(13,*) shelf_depth
            shelf_slope = (basin_depth - shelf_depth) / (x0 - x1)
        endif
        read(13,*)
        
        ! Wind settings
        read(13,*) wind_type
        select case(wind_type)
            case(0) ! No wind
            case(1) ! Constant wind
                read(13,*) A
            case(2) ! Hurricane wind
                read(13,*) A
                read(13,*) B
                read(13,*) Pn
                read(13,*) Pc
                read(13,*) ramp_up_time
                read(13,*) hurricane_velocity
                read(13,*) R_eye_init
            case(3) ! Oscillatory wind
                read(13,*) A
                read(13,*) N
                read(13,*) omega
                read(13,*) t_length
            case(4) ! Gaussian bump
                read(13,*) A
                read(13,*) hurricane_velocity
                read(13,*) R_eye_init
                read(13,*) sigma
            case default
                print *,"Invalid wind field type ", wind_type    
                stop
        end select
        
    end subroutine setup_parameters    
    
    ! Set wind speed based on parameters in this module
    subroutine set_wind(maxmx,mbc,mx,xlower,dx,t,wind)

        implicit none

        ! Input arguments
        integer, intent(in) :: maxmx,mbc,mx
        double precision, intent(in) :: xlower,dx,t

        ! Output
        double precision, intent(inout) :: wind(1-mbc:maxmx+mbc)
    
        ! Local variables
        integer :: i
        double precision :: x,L,time,C,R_eye,r
    
        ! Parse time input, if -99.d99 assume we are at the initial time step and
        ! use ramp_up_time as the initial time if we are trying to do a hurricane
        ! wind field and ramp up, otherwise leave t alone
        if (t == -99.d99) then
            if (wind_type == 2) then
                time = -ramp_up_time
            else
                time = 0.d0
            endif
        else
            time = t
        endif
    
        select case(wind_type)
            ! No wind
            case(0)
                wind = 0.d0
            ! Constant wind
            case(1)
                wind = A
            ! Hurricane wind
            case(2)
                ! Hurrican eye location
                R_eye = t * hurricane_velocity + R_eye_init
        
                ! Parameter constant
                C = 1.d1 * 1.d3**(B/2.d0) * sqrt(A*B*(Pn-Pc)/(rho_air))
        
                ! Set wind field
                do i=1-mbc,mx+mbc
                    ! Distance from center of storm
                    x = xlower + (i-0.5d0) * dx - R_eye
                    r = abs(x)     
                    ! If sufficiently far from eye, calculate wind velocity
                    if (r >= 1.d-3) then
                        wind(i) = sign(1.d0,x) * C * sqrt(exp(-1.d3**B*A/r**B)/r**B)
                    endif
                enddo
                
                ! Ramp up
                if (t < 0.d0) then
                    wind = wind * exp(-(t/(ramp_up_time*0.45d0))**2)
                endif
        
            ! Oscillating wind field
            case(3)
                L = xlower + mx * dx
                do i=1-mbc,mx+mbc
                    x = xlower + (i-0.5d0) * dx
                    wind(i) = A * sin(PI*n*(x)/L) * sin(2*PI*omega/t_length*t)
                enddo
            ! Moving gaussian wind field
            case(4)
                R_eye = t * hurricane_velocity + R_eye_init
                do i=1-mbc,mx+mbc
                    x = xlower + (i-0.5d0) * dx - R_eye
                    wind(i) = A * exp(-((x-R_eye_init)/sigma)**2)
                enddo
        end select
    
    end subroutine set_wind

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


end module parameters_module


subroutine setprob()

    use parameters_module

    implicit none
    
    call setup_parameters()

end subroutine setprob