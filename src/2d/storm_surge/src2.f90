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
    real(kind=8) :: h(num_layers), hu, hv, gamma, dgamma, fdt, a(2,2)
    real(kind=8) :: P_atmos_x, P_atmos_y, tau, wind_speed
    real(kind=8) :: xm, xc, xp, ym, yc, yp, dx_meters, dy_meters

    ! Algorithm parameters
    ! Parameter controls when to zero out the momentum at a depth in the
    ! friction source term
    real(kind=8), parameter :: depth_tolerance = 1.0d-30

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

                ! Apply friction source term only if in shallower water
                if (sum(h) <= friction_depth) then
                    ! Calculate source term
                    gamma = sqrt(hu**2 + hv**2) * (g * manning_coefficient**2) &
                                                / (h(bottom_layer)**(7/3))
                    dgamma = 1.d0 + dt * gamma
                    q(bottom_index + 2, i, j) = q(bottom_index + 2, i, j) / dgamma
                    q(bottom_index + 3, i, j) = q(bottom_index + 3, i, j) / dgamma
                endif
            enddo
        enddo
    endif
    ! End of friction source term

    ! Coriolis source term
    ! TODO: May want to remove the internal calls to coriolis as this could 
    !       lead to slow downs.
    if (coriolis_forcing) then
        do j=1,my
            yc = ylower + (j - 0.5d0) * dy
            fdt = coriolis(yc) * dt ! Calculate f dependent on coordinate system

            ! Calculate matrix components
            a(1,1) = 1.d0 - 0.5d0 * fdt**2 + fdt**4 / 24.d0
            a(1,2) =  fdt - fdt**3 / 6.d0
            a(2,1) = -fdt + fdt**3 / 6.d0
            a(2,2) = a(1,1)

            do i=1,mx
                do m=1,num_layers
                    layer_index = 3 * (m-1)
                    q(layer_index + 2,i,j) = q(layer_index + 2, i, j) * a(1,1) &
                                        + q(layer_index + 3, i, j) * a(1,2)
                    q(layer_index + 3,i,j) = q(layer_index + 2, i, j) * a(2,1) &
                                        + q(layer_index + 3, i, j) * a(2,2)
                enddo
            enddo
        enddo
    endif
    ! End of coriolis source term

    ! wind -----------------------------------------------------------
    if (wind_forcing) then
        ! Force only the top layer of water, assumes top most layer is last
        ! to go dry
        do j=1,my
            do i=1,mx
                if (q(1,i,j) / rho(1) > dry_tolerance(1)) then
                    wind_speed = sqrt(aux(wind_index,i,j)**2 &
                                    + aux(wind_index+1,i,j)**2)
                    if (wind_speed > wind_tolerance) then
                        tau = wind_drag(wind_speed) * rho_air * wind_speed
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
                else
!                     print "('grad(P) = ',2d16.8)",P_atmos_x,P_atmos_y
!                     print "('dt * h = ',1d16.8)", h(1)*dt
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

end subroutine src2
