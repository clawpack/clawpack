! This routine should be a simplified version of src2
! which applies source terms for a 1-d slice of data along the
! edge of a grid.  This is called only from qad where the conservative
! fix-up is applied and is used to apply source terms over partial
! time steps to the coarse grid cell values used in solving Riemann
! problems at the interface between coarse and fine grids.
subroutine src1d(meqn,mbc,mx1d,q1d,maux,aux1d,t,dt)
      
    use storm_module, only: wind_forcing, pressure_forcing
    use storm_module, only: rho_air, wind_drag
    use storm_module, only: pressure_tolerance
      
    use geoclaw_module, only: g => grav, coriolis_forcing, coriolis
    use geoclaw_module, only: friction_forcing, manning_coefficient 
    use geoclaw_module, only: friction_depth, num_layers, rho
    use geoclaw_module, only: omega, coordinate_system

    implicit none

    ! Input
    integer, intent(in) :: meqn, mbc, mx1d, maux
    real(kind=8), intent(in) :: t, dt
    real(kind=8), intent(inout) :: q1d(meqn, mx1d), aux1d(maux, mx1d)

    ! Local storage
    integer :: i, m, bottom_index, bottom_layer, layer_index
    logical :: found
    real(kind=8) :: h(num_layers), hu, hv, gamma, dgamma, y, fdt, a(2,2)

    ! Algorithm parameters
    ! Parameter controls when to zero out the momentum at a depth in the
    ! friction source term
    real(kind=8), parameter :: depth_tolerance = 1.0d-30
    
    ! Friction forcing
    if (friction_forcing > 0) then
        if (friction_forcing == 2) stop "Variable friction unimplemented!"

        do i=1,mx1d

            ! Extract depths
            forall (m=1:num_layers)
                h(m) = q1d(3 * (m-1) + 1,i) / rho(m)
            end forall

            ! Extract appropriate momentum, also zero momentum in dry layers
            m = num_layers
            found = .false.
            do while(.not.found .and. m > 0)
                if (h(m) > depth_tolerance) then
                    ! Extract momentum components and exit loop
                    bottom_layer = m
                    bottom_index = 3 * (m - 1)
                    hu = q1d(bottom_index + 2, i) / rho(m)
                    hv = q1d(bottom_index + 3, i) / rho(m)
                    found = .true.
                else
                    ! Set almost dry layers momentum to zero
                    q1d(3 * (m - 1) + 2, i) = 0.d0
                    q1d(3 * (m - 1) + 3, i) = 0.d0
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
                q1d(bottom_index + 2, i) = q1d(bottom_index + 2, i) / dgamma
                q1d(bottom_index + 3, i) = q1d(bottom_index + 3, i) / dgamma
            endif
        enddo
    endif
    
    ! Only lat-long coordinate system supported here right now
    if (coriolis_forcing .and. coordinate_system == 2) then

        do i=1,mx1d
            ! aux(3,:,:) stores the y coordinates multiplied by deg2rad
            fdt = 2.d0 * omega * sin(aux1d(3,i)) * dt

            ! Calculate matrix components
            a(1,1) = 1.d0 - 0.5d0 * fdt**2 + fdt**4 / 24.d0
            a(1,2) =  fdt - fdt**3 / 6.d0
            a(2,1) = -fdt + fdt**3 / 6.d0
            a(2,2) = a(1,1)
    
            do m=1,num_layers
                q1d(layer_index + 2,i) = q1d(layer_index + 2,i) * a(1,1) &
                                       + q1d(layer_index + 3,i) * a(1,2)
                q1d(layer_index + 3,i) = q1d(layer_index + 2,i) * a(2,1) &
                                       + q1d(layer_index + 3,i) * a(2,2)
            enddo
        enddo
        
    endif
    
    ! = Wind Forcing =========================================================
    if (wind_forcing) then
        ! Force only the top layer of water
        do i=1,mx1d
            if (q1d(1,i) / rho(1) > drytolerance) then
                wind_speed = sqrt(aux1d(4,i)**2 + aux1d(5,i)**2)
                if (wind_speed > wind_tolerance) then
                    tau = wind_drag(wind_speed) * rho_air * wind_speed
                    q1d(2,i) = q1d(2,i) + dt * tau * aux1d(4,i)
                    q1d(3,i) = q1d(3,i) + dt * tau * aux1d(5,i)
                endif
            endif
        enddo
    endif
    ! ========================================================================
    
    ! == Pressure Forcing ====================================================
    if (pressure_forcing .or. .false.) then
        stop "Not sure how to proceed, need direction or the right dx or dy."
    endif

end subroutine src1d
