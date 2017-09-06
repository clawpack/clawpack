subroutine src2(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt)
      
    use geoclaw_module, only: g => grav, coriolis_forcing, coriolis
    use geoclaw_module, only: friction_forcing, friction_depth
    use geoclaw_module, only: manning_coefficient
    use geoclaw_module, only: manning_break, num_manning
    use geoclaw_module, only: spherical_distance, coordinate_system
    use geoclaw_module, only: RAD2DEG, pi, dry_tolerance
    use geoclaw_module, only: rho_air
      
    use storm_module, only: wind_forcing, pressure_forcing, wind_drag
    use storm_module, only: wind_index, pressure_index
    use storm_module, only: storm_direction, storm_location

    use friction_module, only: variable_friction, friction_index

    implicit none
    
    ! Input parameters
    integer, intent(in) :: meqn,mbc,mx,my,maux
    double precision, intent(in) :: xlower,ylower,dx,dy,t,dt
    
    ! Output
    double precision, intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    double precision, intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    ! Locals
    integer :: i, j, nman
    real(kind=8) :: h, hu, hv, gamma, dgamma, y, fdt, a(2,2), coeff
    real(kind=8) :: xm, xc, xp, ym, yc, yp, dx_meters, dy_meters
    real(kind=8) :: u, v, hu0, hv0
    real(kind=8) :: tau, wind_speed, theta, phi, psi, P_gradient(2), S(2)
    real(kind=8) :: Ddt, sloc(2)

    ! Algorithm parameters
    ! Parameter controls when to zero out the momentum at a depth in the
    ! friction source term
    real(kind=8), parameter :: depth_tolerance = 1.0d-30

    ! Physics
    ! Nominal density of water
    real(kind=8), parameter :: rho = 1025.d0

    ! Friction source term
    if (friction_forcing) then
        do j=1,my
            do i=1,mx
                ! Extract appropriate momentum
                if (q(1,i,j) < depth_tolerance) then
                    q(2:3,i,j) = 0.d0
                else
                    ! Apply friction source term only if in shallower water
                    if (q(1,i,j) <= friction_depth) then
                        if (.not.variable_friction) then
                            do nman = num_manning, 1, -1
                                if (aux(1,i,j) .lt. manning_break(nman)) then
                                    coeff = manning_coefficient(nman)
                                endif
                            enddo
                        else
                            coeff = aux(friction_index,i,j)
                        endif
                        
                        ! Calculate source term
                        gamma = sqrt(q(2,i,j)**2 + q(3,i,j)**2) * g     &   
                              * coeff**2 / (q(1,i,j)**(7.d0/3.d0))
                        dgamma = 1.d0 + dt * gamma
                        q(2, i, j) = q(2, i, j) / dgamma
                        q(3, i, j) = q(3, i, j) / dgamma
                    endif
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
            y = ylower + (j - 0.5d0) * dy
            fdt = coriolis(y) * dt ! Calculate f dependent on coordinate system

            ! Calculate matrix components
            a(1,1) = 1.d0 - 0.5d0 * fdt**2 + fdt**4 / 24.d0
            a(1,2) =  fdt - fdt**3 / 6.d0
            a(2,1) = -fdt + fdt**3 / 6.d0
            a(2,2) = a(1,1)

            do i=1,mx
                q(2,i,j) = q(2, i, j) * a(1,1) + q(3, i, j) * a(1,2)
                q(3,i,j) = q(2, i, j) * a(2,1) + q(3, i, j) * a(2,2)
            enddo
        enddo
    endif
    ! End of coriolis source term

    ! wind -----------------------------------------------------------
    if (wind_forcing) then
        ! Need storm location and direction for sector based wind drag
        sloc = storm_location(t)
        theta = storm_direction(t)
        do j=1,my
            yc = ylower + (j - 0.5d0) * dy
            do i=1,mx
                xc = xlower + (i - 0.5d0) * dx
                if (q(1,i,j) > dry_tolerance) then
                    psi = atan2(yc - sloc(2), xc - sloc(1))
                    if (theta > psi) then
                        phi = (2.d0 * pi - theta + psi) * RAD2DEG
                    else
                        phi = (psi - theta) * RAD2DEG 
                    endif
                    wind_speed = sqrt(aux(wind_index,i,j)**2        &
                                    + aux(wind_index+1,i,j)**2)
                    tau = wind_drag(wind_speed, phi) * rho_air * wind_speed / rho
                    q(2,i,j) = q(2,i,j) + dt * tau * aux(wind_index,i,j)
                    q(3,i,j) = q(3,i,j) + dt * tau * aux(wind_index+1,i,j)
                endif
            enddo
        enddo
    endif
    ! ----------------------------------------------------------------

    ! Atmosphere Pressure --------------------------------------------
    ! Handled in Riemann solver
    ! if (pressure_forcing) then
    !     do j=1,my  
    !         ym = ylower + (j - 1.d0) * dy
    !         yc = ylower + (j - 0.5d0) * dy
    !         yp = ylower + j * dy
    !         do i=1,mx  
    !             xm = xlower + (i - 1.d0) * dx
    !             xc = xlower + (i - 0.5d0) * dx
    !             xp = xlower + i * dx
                
    !             if (coordinate_system == 2) then
    !                 ! Convert distance in lat-long to meters
    !                 dx_meters = spherical_distance(xp,yc,xm,yc)
    !                 dy_meters = spherical_distance(xc,yp,xc,ym)
    !             else
    !                 dx_meters = dx
    !                 dy_meters = dy
    !             endif

    !             ! Extract depths
    !             h = q(1,i,j)

    !             ! Calculate gradient of Pressure
    !             P_gradient(1) = (aux(pressure_index,i+1,j) &
    !                            - aux(pressure_index,i-1,j)) / (2.d0 * dx_meters)
    !             P_gradient(2) = (aux(pressure_index,i,j+1) &
    !                            - aux(pressure_index,i,j-1)) / (2.d0 * dy_meters)

    !                 ! Modify momentum in each layer
    !             if (h > dry_tolerance) then
    !                 q(2, i, j) = q(2, i, j) - dt * h * P_gradient(1) / rho
    !                 q(3, i, j) = q(3, i, j) - dt * h * P_gradient(2) / rho
    !             end if
    !         enddo
    !     enddo
    ! endif

end subroutine src2
