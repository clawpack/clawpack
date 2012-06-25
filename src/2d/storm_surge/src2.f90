subroutine src2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt)

    use multilayer_module, only: rho,eta,num_layers
    use hurricane_module
    use geoclaw_module

    implicit none
    
    ! Input parameters
    integer, intent(in) :: maxmx,maxmy,meqn,mbc,mx,my,maux
    double precision, intent(in) :: xlower,ylower,dx,dy,t,dt
    
    ! Ouput
    double precision, intent(inout) :: q(meqn,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc)
    double precision, intent(inout) :: aux(maux,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc)
    
    ! Locals
    integer :: i,j
    double precision :: yc
    ! Friction
    double precision :: h(2),g,coeff,tol,speed,D
    ! Wind and pressure
    double precision :: tau,wind_speed,P_atmos_x,P_atmos_y
    ! Coriolis source term
    ! angular velocity of earth = 2.d0*pi/(86400.d0) 
!     double precision, parameter :: OMEGA = 7.2722052166430395d-05
!     double precision, parameter :: OMEGA = 7.2722052166430395d-03 ! Incorrect value but used for testing
    ! For beta-plane approximation, center of domain's latitude 
    double precision :: fdt,theta,a11,a12,a21,a22,hu,hv
    
    integer :: mn,n

    ! Check for NANs in solution:
!     call check4nans(maxmx,maxmy,meqn,mbc,mx,my,q,t,2)

    g=grav
    coeff = coeffmanning
    tol = 1.d-30  ! To prevent divide by zero in gamma

    ! friction--------------------------------------------------------
    if (coeffmanning > 0.d0 .and. friction_forcing > 0) then
        ! Constant coefficient of friction
        if (ifriction == 1) then
            D = coeff
            do i=1,mx
                do j=1,my
                    ! Check to see which layer we are doing this on
                    if (num_layers > 1) then
                        h(1) = q(1,i,j) / rho(1)
                        h(2) = q(4,i,j) / rho(2)
                    else
                        h(1) = q(1,i,j)
                        h(2) = 0.d0
                    endif
            
                    ! Bottom layer wet, apply to bottom layer
                    if (h(2) > tol) then
                        ! Exactly integrate and modify bottom layer momentum
                        q(5,i,j) = q(5,i,j) * exp(-D*dt)
                        q(6,i,j) = q(6,i,j) * exp(-D*dt)
                
                    ! Only top layer wet, apply to top layer only
                    else if (h(1) > tol) then
                        ! Set bottom layer momentum to zero
                        if (num_layers > 1) q(5:6,i,j) = 0.d0
                        
                        ! Exactly integrate and modify top layer momentum
                        q(2,i,j) = q(2,i,j) * exp(-D*dt)
                        q(3,i,j) = q(3,i,j) * exp(-D*dt)
                
                    ! Neither layer wet, set momentum to zero
                    else
                        q(2:3,i,j) = 0.d0
                        if (num_layers > 1) q(5:6,i,j) = 0.d0
                    endif
                enddo
            enddo
        ! Manning-N friction
        else if (ifriction == 2) then
            do i=1,mx
                do j=1,my
                    ! Check to see which layer we are doing this on
                    if (num_layers > 1) then
                        h(1) = q(1,i,j) / rho(1)
                        h(2) = q(4,i,j) / rho(2)
                    else
                        h(1) = q(1,i,j)
                        h(2) = 0.d0
                    endif
            
                    ! Bottom layer wet, apply to bottom layer
                    if (h(2) > tol) then
                        ! Extract speed of bottom layer
                        speed = sqrt(q(5,i,j)**2 + q(6,i,j)**2) / q(4,i,j)
                                    
                        ! Calculate drag coefficient 
                        D = coeff**2 * g * sum(h)**(-7/3) * speed
                
                        ! Exactly integrate and modify bottom layer momentum
                        q(5,i,j) = q(5,i,j) * exp(-D*dt)
                        q(6,i,j) = q(6,i,j) * exp(-D*dt)
                
                    ! Only top layer wet, apply to top layer only
                    else if (h(1) > tol) then
                        ! Set bottom layer momentum to zero
                        if (num_layers > 1) q(5:6,i,j) = 0.d0
                
                        ! Extract speed of top layer
                        speed = sqrt(q(2,i,j)**2 + q(3,i,j)**2) / q(1,i,j)
                
                        ! Calculate drag coefficient
                        D = coeff**2 * g * sum(h)**(-7/3) * speed
                        
                        ! Exactly integrate and modify top layer momentum
                        q(2,i,j) = q(2,i,j) * exp(-D*dt)
                        q(3,i,j) = q(3,i,j) * exp(-D*dt)
                
                    ! Neither layer wet, set momentum to zero
                    else
                        q(2:3,i,j) = 0.d0
                        if (num_layers > 1) q(5:6,i,j) = 0.d0
                    endif
                enddo
            enddo
        endif
    endif
    ! ----------------------------------------------------------------

    ! coriolis--------------------------------------------------------
    if (icoriolis > 0) then
        do j=1,my
            yc = ylower + (j-.5d0)*dy
            fdt = coriolis(yc) * dt
            do i=1,mx
                !dq/dt = 2w*sin(latitude)*[0 1 ; -1 0] q = Aq
                !e^Adt = [a11 a12; a21 a22] + I
                a11 = 1.d0 - 0.5d0 * fdt**2 + fdt**4 / 24.d0
                a12 = fdt - fdt**3 / 6.0d0
                a21 = -fdt + fdt**3 / 6.0d0
                a22 = a11
                hu = q(2,i,j)
                hv = q(3,i,j)
                !q = e^Adt * q0
                q(2,i,j) = hu * a11 + hv * a12
                q(3,i,j) = hu * a21 + hv * a22
                if (num_layers > 1) then
                    hu = q(5,i,j)
                    hv = q(6,i,j)
                    q(5,i,j) = hu * a11 + hv * a12
                    q(6,i,j) = hu * a21 + hv * a22
                endif
            enddo
        enddo
    endif
    ! ----------------------------------------------------------------

    ! wind -----------------------------------------------------------
    if (wind_forcing) then
        ! Here we have to take into account geoclaw which we use for the 
        ! single layer case.  It needs to divide by the water's density
        if (num_layers > 1) then
            do i=1,mx
                do j=1,my
                    if (q(1,i,j) / rho(1) > drytolerance) then
                        wind_speed = sqrt(aux(4,i,j)**2 + aux(5,i,j)**2)
                        if (wind_speed > wind_tolerance) then
                            tau = wind_drag(wind_speed) * rho_air * wind_speed
                            q(2,i,j) = q(2,i,j) + dt * tau * aux(4,i,j)
                            q(3,i,j) = q(3,i,j) + dt * tau * aux(5,i,j)
                        endif
                    endif
                enddo
            enddo
        else
            do i=1,mx
                do j=1,my
    !                 if (abs(q(i,j,1)) > drytolerance) then
                        wind_speed = sqrt(aux(4,i,j)**2 + aux(5,i,j)**2)
                        tau = wind_drag(wind_speed) * rho_air * wind_speed
                        q(2,i,j) = q(2,i,j) + dt * tau * aux(4,i,j) / rho(1)
                        q(3,i,j) = q(3,i,j) + dt * tau * aux(5,i,j) / rho(1)
    !                 endif
                enddo
            enddo
        endif
    endif
    ! ----------------------------------------------------------------

    ! atmosphere -----------------------------------------------------
    if (pressure_forcing) then
        do i=1,mx
            do j=1,my                  
                h = 0.d0  
                if (num_layers > 1) then
                    h(1) = q(1,i,j) / rho(1)
                    h(2) = q(4,i,j) / rho(2)
                else
                    h(1) = q(1,i,j)
                endif
                
                ! Calculate gradient of Pressure
                P_atmos_x = (aux(6,i+1,j) - aux(6,i-1,j)) / (2.d0*dx)
                P_atmos_y = (aux(6,i,j+1) - aux(6,i,j-1)) / (2.d0*dy)
                if (abs(P_atmos_x) < pressure_tolerance) then
                    P_atmos_x = 0.d0
                endif
                if (abs(P_atmos_y) < pressure_tolerance) then
                    P_atmos_y = 0.d0
                endif
                
                if (num_layers > 1) then
                    if (h(1) > drytolerance) then
                        q(2,i,j) = q(2,i,j) - dt * h(1) * P_atmos_x
                        q(3,i,j) = q(3,i,j) - dt * h(1) * P_atmos_y
                    else if (h(2) > drytolerance) then
                        q(5,i,j) = q(5,i,j) - dt * h(2) * P_atmos_x
                        q(6,i,j) = q(6,i,j) - dt * h(2) * P_atmos_y
                    endif
                else
                    if (h(1) > drytolerance) then
                        q(2,i,j) = q(2,i,j) - dt * h(1) * P_atmos_x / rho(1)
                        q(3,i,j) = q(3,i,j) - dt * h(1) * P_atmos_y / rho(1)
                    endif
                endif
            enddo
        enddo
    endif
    ! ----------------------------------------------------------------
end subroutine src2