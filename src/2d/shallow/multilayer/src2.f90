subroutine src2(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt)
      
    use storm_module, only: wind_forcing, pressure_forcing
    use storm_module, only: wind_drag
    use storm_module, only: wind_index, pressure_index
    use storm_module, only: storm_direction, storm_location

    use friction_module, only: friction_index

    use geoclaw_module, only: g => grav, RAD2DEG, pi, rho_air, ambient_pressure
    use geoclaw_module, only: coriolis_forcing, coriolis, rho
    use geoclaw_module, only: friction_forcing, friction_depth
    use geoclaw_module, only: manning_coefficient
    use geoclaw_module, only: manning_break, num_manning
    use geoclaw_module, only: spherical_distance, coordinate_system

    use amr_module, only: mcapa, cflv1
      
    use multilayer_module, only: num_layers, dry_tolerance
    
    use friction_module, only: variable_friction, friction_index

    implicit none
    
    ! Input parameters
    integer, intent(in) :: meqn,mbc,mx,my,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy,t,dt
    
    ! Output
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    ! Locals
    integer :: i, j, m, k, bottom_index, bottom_layer, layer, nman
    logical :: found
    real(kind=8) :: xm, xc, xp, ym, yc, yp, dx_meters, dy_meters
    real(kind=8) :: h(num_layers), hu, hv, u, v, hu0, hv0, coeff
    real(kind=8) :: tau, wind_speed, theta, phi, psi, P_gradient(2), S(2), sloc(2)
    real(kind=8) :: fdt, Ddt, a(2,2)

    ! Physics parameters
    real(kind=8), parameter :: manning_power = 7.d0 / 3.d0
    real(kind=8), parameter :: H_break = 2.d0
    real(kind=8), parameter :: theta_f = 10.d0
    real(kind=8), parameter :: gamma_f = 4.d0 / 3.d0

    ! Algorithm parameters
    ! Parameter controls when to zero out the momentum at a depth in the
    ! friction source term
    real(kind=8), parameter :: depth_tolerance = 1.0d-30

    character(len=*), parameter :: CFL_FORMAT =   &
                        "('*** WARNING *** Courant number  =', d12.4," // &
                         "'  is larger than input cfl_max = ', d12.4," // &
                         "' in src2.')"

    ! ----------------------------------------------------------------
    ! Friction source term - Use backward Euler to solve
    ! Hybrid friction formula with a spatially varying Manning's-N factor
    if (friction_forcing) then
        do j=1,my
            do i=1,mx
                ! Extract depths
                forall (m=1:num_layers)
                    h(m) = q(3 * (m-1) + 1,i,j) / rho(m)
                end forall

                ! Extract appropriate momentum, also zero momentum in dry layers
                m = num_layers
                found = .false.
                do while(.not.found .and. m > 0)
                    if (h(m) > depth_tolerance) then
                        ! Extract momentum components and exit loop
                        bottom_layer = m
                        bottom_index = 3 * (m - 1)
                        hu = q(bottom_index + 2, i, j) / rho(m)
                        hv = q(bottom_index + 3, i, j) / rho(m)
                        found = .true.
                    else
                        ! Set almost dry layers momentum to zero
                        q(3 * (m - 1) + 2, i, j) = 0.d0
                        q(3 * (m - 1) + 3, i, j) = 0.d0
                    endif
                    m = m - 1
                end do

                if (.not.found) then
                    cycle
                end if

                if (variable_friction) then
                    coeff = aux(friction_index, i, j)
                else
                    do nman = num_manning, 1, -1
                        if (aux(1,i,j) < manning_break(nman)) then
                            coeff = manning_coefficient(nman)
                        endif
                    enddo
                end if

                ! Apply friction if in shallow enough water
                if (sum(h) <= friction_depth) then
                    tau = g * coeff**2 / h(bottom_layer)**(manning_power) &
                            * sqrt(hu**2 + hv**2) &
                            * (1.d0 + (H_break / h(bottom_layer))**theta_f)**(gamma_f / theta_f)
                            
                    q(bottom_index + 2,i,j) = q(bottom_index + 2,i,j) / (1.d0 + dt * tau)
                    q(bottom_index + 3,i,j) = q(bottom_index + 3,i,j) / (1.d0 + dt * tau)
                endif
            enddo
        enddo
    endif
    ! End of friction source term

    ! Coriolis source term
    ! Use backward Euler to solve, q_t = A q -> q^n+1 = (I + dt * A)^-1 q^n
    ! or use first 4 terms of matrix exponential
    if (coriolis_forcing) then
        do j=1,my
            yc = ylower + (j - 0.5d0) * dy
            fdt = coriolis(yc) * dt ! Calculate f dependent on coordinate system

            ! Calculate matrix components - backward Euler
            !a(1,:) = [1.d0,  fdt]
            !a(2,:) = [-fdt, 1.d0]
            !a = a / (1.d0 + fdt**2)

            ! Matrix exponential
            a(1,1) = 1.d0 - 0.5d0 * fdt**2 + fdt**4 / 24.d0
            a(1,2) =  fdt - fdt**3 / 6.d0
            a(2,1) = -fdt + fdt**3 / 6.d0
            a(2,2) = a(1,1)

            do i=1,mx
                do m=1,num_layers
                    layer = 3 * (m-1)
                    hu = q(layer + 2,i,j)
                    hv = q(layer + 3,i,j)
                    ! Note that these are not really hu, hv but rho * hu and 
                    ! rho * hv actually
                    q(2,i,j) = hu * a(1,1) + hv * a(1,2)
                    q(3,i,j) = hu * a(2,1) + hv * a(2,2)
                end do
            enddo
        enddo
    endif
    ! End of coriolis source term

    ! wind -----------------------------------------------------------
    if (wind_forcing) then
        ! Force only the top layer of water, assumes top most layer is last
        ! to go dry
        do j=1,my
            yc = ylower + (j - 0.5d0) * dy
            do i=1,mx
                xc = xlower + (i - 0.5d0) * dx
                if (q(1,i,j) / rho(1) > dry_tolerance(1)) then
                    psi = atan2(yc - sloc(2), xc - sloc(1))
                    if (theta > psi) then
                        phi = (2.d0 * pi - theta + psi) * RAD2DEG
                    else
                        phi = (psi - theta) * RAD2DEG 
                    endif
                    wind_speed = sqrt(aux(wind_index,i,j)**2        &
                                    + aux(wind_index+1,i,j)**2)
                    tau = wind_drag(wind_speed, phi) * rho_air * wind_speed
                    q(2,i,j) = q(2,i,j) + dt * tau * aux(wind_index,i,j)
                    q(3,i,j) = q(3,i,j) + dt * tau * aux(wind_index+1,i,j)
                endif
            enddo
        enddo
    endif
    ! ----------------------------------------------------------------

    ! Atmosphere Pressure --------------------------------------------
    ! Not yet handled in Riemann solver
    if (pressure_forcing) then
        do j=1,my  
            ym = ylower + (j - 1.d0) * dy
            yc = ylower + (j - 0.5d0) * dy
            yp = ylower + j * dy
            do i=1,mx  
                xm = xlower + (i - 1.d0) * dx
                xc = xlower + (i - 0.5d0) * dx
                xp = xlower + i * dx
                
                if (coordinate_system == 2) then
                    ! Convert distance in lat-long to meters
                    dx_meters = spherical_distance(xp,yc,xm,yc)
                    dy_meters = spherical_distance(xc,yp,xc,ym)
                else
                    dx_meters = dx
                    dy_meters = dy
                endif

                ! Calculate gradient of Pressure
                P_gradient(1) = (aux(pressure_index,i+1,j) &
                               - aux(pressure_index,i-1,j)) / (2.d0 * dx_meters)
                P_gradient(2) = (aux(pressure_index,i,j+1) &
                               - aux(pressure_index,i,j-1)) / (2.d0 * dy_meters)

                ! Modify momentum in each layer
                do k=1,num_layers
                    layer = 3 * (k - 1)
                    h(k) = q(layer + 1, i, j) / rho(k)
                    if (h(k) > dry_tolerance(k)) then
                        q(layer + 2, i, j) = q(layer + 2, i, j)   &
                                    - dt * h(k) * P_gradient(1)
                        q(layer + 3, i, j) = q(layer + 3, i, j)   &
                                    - dt * h(k) * P_gradient(2)
                    end if
                end do
            enddo
        enddo
    endif

end subroutine src2