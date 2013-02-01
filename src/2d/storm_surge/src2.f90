subroutine src2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt)
      
    use storm_module, only: wind_forcing, pressure_forcing
    use storm_module, only: rho_air, wind_drag
    use storm_module, only: pressure_tolerance, wind_tolerance
    use storm_module, only: wind_index, pressure_index

    use friction_module, only: friction_index

    use geoclaw_module, only: g => grav, dry_tolerance
    use geoclaw_module, only: coriolis_forcing, coriolis
    use geoclaw_module, only: friction_forcing
    use geoclaw_module, only: spherical_distance, coordinate_system

    use geoclaw_module, only: pi, omega

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
    real(kind=8) :: h, hu, hv, fdt, gamma, dgamma, a(2,2)
    real(kind=8) :: P_atmos_x, P_atmos_y, tau, wind_speed
    real(kind=8) :: xm, xc, xp, ym, yc, yp, dx_meters, dy_meters

    ! Physics parameters
    real(kind=8), parameter :: rho = 1025.d0
    real(kind=8), parameter :: H_break = 2.d0
    real(kind=8), parameter :: theta_f = 10.d0
    real(kind=8), parameter :: gamma_f = 4.d0 / 3.d0

    ! Algorithm parameters
    ! Parameter controls when to zero out the momentum at a depth in the
    ! friction source term
    real(kind=8), parameter :: depth_tolerance = 1.0d-30

    ! wind -----------------------------------------------------------
    if (wind_forcing) then
        ! Force only the top layer of water, assumes top most layer is last
        ! to go dry
        do j=1,my
            do i=1,mx
                if (q(1,i,j) > dry_tolerance) then
                    wind_speed = sqrt(aux(wind_index,i,j)**2 &
                                    + aux(wind_index+1,i,j)**2)
                    if (wind_speed > wind_tolerance) then
                        tau = wind_drag(wind_speed) * rho_air * wind_speed / rho
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
                    dx_meters = spherical_distance(xp,yc,xm,yc)
                    dy_meters = spherical_distance(xc,yp,xc,ym)
                else
                    dx_meters = dx
                    dy_meters = dy
                endif

                ! Extract depths
                h = q(1,i,j)

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
                if (h > dry_tolerance) then
                    q(2, i, j) = q(2, i, j) - dt * h * P_atmos_x / rho
                    q(3, i, j) = q(3, i, j) - dt * h * P_atmos_y / rho
                end if
            enddo
        enddo
    endif
    
    ! Coriolis source term
    ! Use backward Euler to solve, q_t = A q -> q^n+1 = (I + dt * A)^-1 q^n
    if (coriolis_forcing) then
        do j=1,my
            yc = ylower + (j - 0.5d0) * dy
            fdt = coriolis(yc) * dt ! Calculate f dependent on coordinate system

            ! Calculate matrix components
            a(1,1) = -fdt**2 / (fdt**2 + 1.d0) + 1.d0
            a(1,2) = -fdt / (fdt**2 + 1.d0)
            a(2,1) = fdt / (fdt**2 + 1.d0)
            a(2,2) = 1 / (fdt**2 + 1.d0)

            do i=1,mx
                hu = q(2,i,j)
                hv = q(3,i,j)
                q(2,i,j) = hu * a(1,1) + hv * a(1,2)
                q(3,i,j) = hu * a(2,1) + hv * a(2,2)
            enddo
        enddo
    endif
    ! End of coriolis source term

    ! ----------------------------------------------------------------
    ! Friction source term - Use backward Euler to solve
    ! Hybrid friction formula with a spatially varying Manning's-N factor
    if (friction_forcing) then
        do j=1,my
            do i=1,mx
                if (q(1,i,j) < depth_tolerance) then
                    q(2:3,i,j) = 0.d0
                else
                    tau = g * aux(friction_index,i,j)**2 / q(1,i,j)**(7.d0 / 3.d0) &
                          * (1 + (H_break / q(1,i,j))**theta_f)**(gamma_f / theta_f) &
                          * sqrt(q(2,i,j)**2 + q(3,i,j)**2)
                    q(2,i,j) = q(2,i,j) / (1.d0 + dt * tau)
                    q(3,i,j) = q(3,i,j) / (1.d0 + dt * tau)
                endif
            enddo
        enddo
    endif

    ! End of friction source term

end subroutine src2
