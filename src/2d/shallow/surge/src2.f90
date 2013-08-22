subroutine src2(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt)
      
    use storm_module, only: wind_forcing, pressure_forcing
    use storm_module, only: rho_air, wind_drag, ambient_pressure
    use storm_module, only: wind_index, pressure_index
    use storm_module, only: storm_direction, storm_location

    use friction_module, only: friction_index

    use geoclaw_module, only: g => grav, dry_tolerance, RAD2DEG, pi
    use geoclaw_module, only: coriolis_forcing, coriolis
    use geoclaw_module, only: friction_forcing, manning_coefficient
    use geoclaw_module, only: spherical_distance, coordinate_system

    use amr_module, only: mcapa, cflv1

    implicit none
    
    ! Input parameters
    integer, intent(in) :: meqn,mbc,mx,my,maux
    double precision, intent(in) :: xlower,ylower,dx,dy,t,dt
    
    ! Output
    double precision, intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    double precision, intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    ! Locals
    integer :: i, j
    real(kind=8) :: xm, xc, xp, ym, yc, yp, dx_meters, dy_meters
    real(kind=8) :: h, hu, hv, u, v, hu0, hv0
    real(kind=8) :: tau, wind_speed, theta, phi, psi, P_gradient(2), S(2), sloc(2)
    real(kind=8) :: fdt, Ddt, a(2,2)

    ! Physics parameters
    real(kind=8), parameter :: fric_coefficient = 7.d0 / 3.d0
    real(kind=8), parameter :: rho = 1025.d0
    real(kind=8), parameter :: H_break = 2.d0
    real(kind=8), parameter :: theta_f = 10.d0
    real(kind=8), parameter :: gamma_f = 4.d0 / 3.d0

    ! Algorithm parameters
    ! Here to prevent a divide by zero in friction term
    real(kind=8), parameter :: friction_tolerance = 1.0d-30

    character(len=*), parameter :: CFL_FORMAT =   &
                        "('*** WARNING *** Courant number  =', d12.4," // &
                         "'  is larger than input cfl_max = ', d12.4," // &
                         "' in src2.')"

    real(kind=8) :: dtdx, dtdy, cfl

    ! CFL check parameters
    dtdx = dt / dx
    dtdy = dt / dy

    ! Get storm direction for this time
    theta = storm_direction(t)
    sloc = storm_location(t)

    ! Primary loops over grid
!     do j=1,my
!         ym = ylower + (j - 1.d0) * dy
!         yc = ylower + (j - 0.5d0) * dy
!         yp = ylower + j * dy
        
!         fdt = coriolis(yc) * dt * 0.5d0
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

!             ! State
!             h = q(1,i,j)
!             if (h < dry_tolerance) then
!                 q(2,i,j) = 0.d0
!                 q(3,i,j) = 0.d0
!             else
!                 hu0 = q(2,i,j)
!                 hv0 = q(3,i,j)
!                 S = 0.d0

!                 ! Wind
!                 psi = atan2(yc - sloc(2), xc - sloc(1))
!                 if (theta > psi) then
!                     phi = (2.d0 * pi - theta + psi) * RAD2DEG
!                 else
!                     phi = (psi - theta) * RAD2DEG 
!                 endif
!                 wind_speed = sqrt(aux(wind_index,i,j)**2 + aux(wind_index+1,i,j)**2)
!                 tau = wind_drag(wind_speed, phi) * rho_air * wind_speed / rho
!                 S = S + tau * aux(wind_index:wind_index+1,i,j)

!                 ! Pressure
!                 P_gradient(1) = (aux(pressure_index,i+1,j) &
!                                 - aux(pressure_index,i-1,j)) / (2.d0 * dx_meters)
!                 P_gradient(2) = (aux(pressure_index,i,j+1) &
!                                 - aux(pressure_index,i,j-1)) / (2.d0 * dy_meters)
!                 S = S - h * P_gradient / rho
!                 aux(pressure_index+1,i,j) = sqrt(S(1)**2 + S(2)**2)
! !                 aux(pressure_index+1,i,j) = wind_drag(wind_speed, phi)
! !                 if (20.d0 < phi .and. phi <= 150.d0) then
! !                     aux(pressure_index+2,i,j) = 1.d0
! !                 else if (150.d0 < phi .and. phi <= 240.d0) then
! !                     aux(pressure_index+2,i,j) = 2.d0
! !                 else
! !                     aux(pressure_index+2,i,j) = 3.d0
! !                 endif
! !                 if (310.d0 < phi .and. phi <= 360.d0) then
! !                     aux(pressure_index+2,i,j) = 1.d0
! !                 else if (0.d0 < phi .and. phi <= 85.d0) then
! !                     aux(pressure_index+2,i,j) = 2.d0
! !                 ! Right - Rear sector
! !                 else if (85.d0 < phi .and. phi <= 195.d0) then
! !                     aux(pressure_index+2,i,j) = 3.d0

! !                 ! Rear - Left sector
! !                 else if (195.d0 < phi .and. phi <= 310.d0) then
! !                     aux(pressure_index+2,i,j) = 4.d0
! !                 endif

!                 ! Friction and Coriolis
!                 Ddt = g * aux(friction_index,i,j)**2 / h**(7.d0 / 3.d0) &
!                         * (1.d0 + (H_break / h)**theta_f)**(gamma_f / theta_f) &
!                         * sqrt(hu0**2 + hv0**2) * dt * 0.5d0

!                 ! Matrix I + dt / 2 A
!                 a(1,:) = [1.d0 - Ddt, fdt]
!                 a(2,:) = [-fdt, 1.d0 - Ddt]
!                 hu = a(1,1) * hu0 + a(1,2) * hv0 + dt * S(1)
!                 hv = a(2,1) * hu0 + a(2,2) * hv0 + dt * S(2)

!                 ! Matrix (I - dt A / 2)^-1
!                 a(1,:) = [1.d0 + Ddt, -fdt]
!                 a(2,:) = [fdt, 1.d0 + Ddt]
!                 a = a / ((1.d0 + Ddt)**2 + fdt**2)

!                 ! Update q
!                 q(2,i,j) = hu * a(1,1) + hv * a(1,2)
!                 q(3,i,j) = hu * a(2,1) + hv * a(2,2)

!                 ! Check CFL estimate
!                 if (mcapa > 0) then
!                     dtdx = dtdx / aux(mcapa,i,j)
!                     dtdy = dtdy / aux(mcapa,i,j)
!                 endif
!                 u = q(2,i,j) / q(1,i,j)
!                 v = q(3,i,j) / q(1,i,j)
!                 cfl = max(dtdx * (u - sqrt(g * h)), dtdx * (u + sqrt(g * h)), &
!                           dtdy * (v - sqrt(g * h)), dtdy * (v + sqrt(g * h)))
!                 if (cfl > cflv1) then
!                     print CFL_FORMAT, cfl,cflv1
!                     print "('h,u,v(',i3,',',i3,') = ',4d12.4)", i,j,h,u,v
!                 endif

!                 aux(pressure_index+2,i,j) = sqrt((q(2,i,j) - hu0 - dt*S(1))**2  &
!                                                + (q(3,i,j) - hv0 - dt*S(2))**2)
!             endif
!         enddo
!     enddo

    ! wind -----------------------------------------------------------
    if (wind_forcing) then
        ! Force only the top layer of water, assumes top most layer is last
        ! to go dry
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
!                     if (wind_speed > wind_tolerance) then
!                         tau = wind_drag(wind_speed) * rho_air * wind_speed / rho
                        q(2,i,j) = q(2,i,j) + dt * tau * aux(wind_index,i,j)
                        q(3,i,j) = q(3,i,j) + dt * tau * aux(wind_index+1,i,j)
!                     endif
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
!                 if (abs(aux(pressure_index,i,j) - ambient_pressure)      &
!                                                     >= pressure_tolerance) then
                P_gradient(1) = (aux(pressure_index,i+1,j) &
                               - aux(pressure_index,i-1,j)) / (2.d0 * dx_meters)
                P_gradient(2) = (aux(pressure_index,i,j+1) &
                               - aux(pressure_index,i,j-1)) / (2.d0 * dy_meters)

                    ! Modify momentum in each layer
                if (h > dry_tolerance) then
                    q(2, i, j) = q(2, i, j) - dt * h * P_gradient(1) / rho
                    q(3, i, j) = q(3, i, j) - dt * h * P_gradient(2) / rho
                end if
!                 endif
            enddo
        enddo
    endif

    ! ----------------------------------------------------------------
    ! Friction source term - Use backward Euler to solve
    ! Hybrid friction formula with a spatially varying Manning's-N factor
    if (friction_forcing) then
        do j=1,my
            do i=1,mx
                if (q(1,i,j) < friction_tolerance) then
                    q(2:3,i,j) = 0.d0
                else
                    tau = g * aux(friction_index,i,j)**2 / q(1,i,j)**(fric_coefficient) &
                            * (1.d0 + (H_break / q(1,i,j))**theta_f)**(gamma_f / theta_f) &
                            * sqrt(q(2,i,j)**2 + q(3,i,j)**2)
                    q(2,i,j) = q(2,i,j) / (1.d0 + dt * tau)
                    q(3,i,j) = q(3,i,j) / (1.d0 + dt * tau)

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
                hu = q(2,i,j)
                hv = q(3,i,j)
                q(2,i,j) = hu * a(1,1) + hv * a(1,2)
                q(3,i,j) = hu * a(2,1) + hv * a(2,2)
            enddo
        enddo
    endif
    ! End of coriolis source term

end subroutine src2
