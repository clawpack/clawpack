subroutine src2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt)
      
    use storm_module, only: wind_forcing, pressure_forcing
    use storm_module, only: rho_air, wind_drag
    use storm_module, only: pressure_tolerance, wind_tolerance
    use storm_module, only: wind_index, pressure_index

    use geoclaw_module, only: g => grav, dry_tolerance
    use geoclaw_module, only: coriolis_forcing, coriolis
    use geoclaw_module, only: friction_forcing, manning_coefficient 
    use geoclaw_module, only: friction_depth, num_layers, rho
    use geoclaw_module, only: spherical_distance, coordinate_system

    implicit none
    
    ! Input parameters
    integer, intent(in) :: maxmx,maxmy,meqn,mbc,mx,my,maux
    double precision, intent(in) :: xlower,ylower,dx,dy,t,dt
    
    ! Output
    double precision, intent(inout) :: q(meqn,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc)
    double precision, intent(inout) :: aux(maux,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc)

    ! Locals
    integer :: i, j, m, bottom_index, bottom_layer, layer_index
    logical :: found
    real(kind=8) :: h(num_layers), hu, hv, gamma, fdt(my), a(2,2)
    real(kind=8) :: P_atmos_x, P_atmos_y, tau, wind_speed
    real(kind=8) :: xm, xc, xp, ym, yc, yp, dx_meters, dy_meters

    real(kind=8) :: D, speed

    ! Chezy formulation for friction    
!     real(kind=8), parameter :: D_min = 3.0d-3
    real(kind=8), parameter :: H_break = 2.d0
    real(kind=8), parameter :: theta_f = 10.d0
    real(kind=8), parameter :: gamma_f = 1.3333d0

    ! Algorithm parameters
    ! Parameter controls when to zero out the momentum at a depth in the
    ! friction source term
    real(kind=8), parameter :: depth_tolerance = 0.5

    ! Coriolis source term
    if (coriolis_forcing) then
        
        ! Calculate coriolis term
        forall(j=1:my)
            fdt(j) = coriolis(ylower + (j-0.5d0) * dy) * dt
        end forall

        ! Loop through grid points
        do j=1,my
            do i=1,mx
                do m=1,num_layers
                    layer_index = 3 * (m-1)
                    h(m) = q(layer_index + 1,i,j) / rho(1)

                    if (h(m) > dry_tolerance(m)) then

                        ! Calculate matrix components
                        a(1,1) = 1.d0 - 0.5d0 * fdt(j)**2 + fdt(j)**4 / 24.d0
                        a(1,2) =  fdt(j) - fdt(j)**3 / 6.d0
                        a(2,1) = -fdt(j) + fdt(j)**3 / 6.d0
                        a(2,2) = a(1,1)

                        hu = q(layer_index + 2,i,j) / rho(m)
                        hv = q(layer_index + 3,i,j) / rho(m)

                        q(layer_index + 2,i,j) = &
                                            (hu * a(1,1) + hv * a(1,2)) * rho(m)
                        q(layer_index + 3,i,j) = &
                                            (hu * a(2,1) + hv * a(2,2)) * rho(m)
                    end if
                enddo
            enddo
        enddo
    endif
    ! End of coriolis source term

    ! wind -----------------------------------------------------------
    if (wind_forcing) then
        ! Force only the top layer of water, assumes top most layer is last
        ! to go dry
        aux(pressure_index+1,i,j) = 0.d0
        aux(pressure_index+2,i,j) = 0.d0
        do j=1,my
            do i=1,mx
                if (q(1,i,j) / rho(1) > dry_tolerance(1)) then
                    wind_speed = sqrt(aux(wind_index,i,j)**2 &
                                    + aux(wind_index+1,i,j)**2)
                    if (wind_speed > wind_tolerance) then
                        tau = wind_drag(wind_speed) * rho_air * wind_speed
                        aux(pressure_index+1,i,j) = dt * tau * aux(wind_index,i,j)
                        aux(pressure_index+2,i,j) = dt * tau * aux(wind_index+1,i,j)
                        q(2,i,j) = q(2,i,j) + dt * tau * aux(wind_index,i,j)
                        q(3,i,j) = q(3,i,j) + dt * tau * aux(wind_index+1,i,j)
                    endif
                endif
            enddo
        enddo
    endif
    ! ----------------------------------------------------------------

    ! Atmosphere Pressure --------------------------------------------
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
                    dx_meters = spherical_distance([xp,yc],[xm,yc])
                    dy_meters = spherical_distance([xc,yp],[xc,ym])
                else
                    dx_meters = dx
                    dy_meters = dy
                endif

                ! Extract depths
                forall (m=1:num_layers)
                    h(m) = q(3 * (m-1) + 1,i,j) / rho(m)
                end forall

                ! Calculate gradient of Pressure
                P_atmos_x = (aux(pressure_index,i+1,j) &
                                    - aux(pressure_index,i-1,j)) / (2.d0 * dx_meters)
                P_atmos_y = (aux(pressure_index,i,j+1) &
                                    - aux(pressure_index,i,j-1)) / (2.d0 * dy_meters)
                if (abs(P_atmos_x) < pressure_tolerance) then
                    P_atmos_x = 0.d0
                endif
                if (abs(P_atmos_y) < pressure_tolerance) then
                    P_atmos_y = 0.d0
                endif

                ! Modify momentum in each layer
                do m=1,num_layers
                    if (h(m) > dry_tolerance(m)) then
                        layer_index = 3 * (m - 1)

                        q(layer_index + 2, i, j) = q(layer_index + 2, i, j) &
                                                        - dt * h(m) * P_atmos_x
                        q(layer_index + 3, i, j) = q(layer_index + 3, i, j) &
                                                        - dt * h(m) * P_atmos_y
                    end if
                enddo
            enddo
        enddo
    endif
    ! ----------------------------------------------------------------

    ! Friction source term
    if (friction_forcing > 0) then

        ! Parameter checks
        if (friction_forcing == 2) stop "Variable friction unimplemented!"
        if (depth_tolerance > friction_depth) then
            stop "Parameter depth_tolerance > friction_depth!"
        endif

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
                endif

                ! Apply friction source term only if in shallower water, this
                ! used backward Euler for the source term
                if (sum(h) <= friction_depth) then
                    speed = sqrt(hu**2 + hv**2) / h(bottom_layer)

                    ! Non-limited version
                    ! D = manning_coefficient**2 * g * sum(h)**(-7/3) * speed
                    ! Limited Chezy-type coefficient
                    D = g * manning_coefficient**2 * h(bottom_layer)**(-4/3) * speed * &
                        (1.d0 + (H_break / h(bottom_layer))**theta_f )**(gamma_f / theta_f)
                    gamma = 1.d0 + dt * D

                    q(bottom_index + 2:bottom_index + 3, i ,j) = &
                        q(bottom_index + 2:bottom_index + 3, i ,j) / gamma

                endif
            enddo
        enddo
    endif
    ! End of friction source term

end subroutine src2
