! This routine should be a simplified version of src2
! which applies source terms for a 1-d slice of data along the
! edge of a grid.  This is called only from qad where the conservative
! fix-up is applied and is used to apply source terms over partial
! time steps to the coarse grid cell values used in solving Riemann
! problems at the interface between coarse and fine grids.
subroutine src1d(meqn,mbc,mx1d,q1d,maux,aux1d,t,dt)
      
    use geoclaw_module, only: g => grav, coriolis_forcing, coriolis
    use geoclaw_module, only: friction_forcing, friction_depth
    use geoclaw_module, only: omega, coordinate_system, manning_coefficient
    use geoclaw_module, only: manning_break, num_manning
                              

    implicit none

    ! Input
    integer, intent(in) :: meqn, mbc, mx1d, maux
    real(kind=8), intent(in) :: t, dt
    real(kind=8), intent(inout) :: q1d(meqn, mx1d), aux1d(maux, mx1d)

    ! Local storage
    integer :: i, nman
    logical :: found
    real(kind=8) :: h, hu, hv, gamma, dgamma, y, fdt, a(2,2), coeff

    ! Algorithm parameters
    ! Parameter controls when to zero out the momentum at a depth in the
    ! friction source term
    real(kind=8), parameter :: depth_tolerance = 1.0d-30
    
    ! Friction forcing
    if (friction_forcing) then

        do i=1,mx1d

            ! Extract depths
            h = q1d(1,i)
            hu = q1d(2,i)
            hv = q1d(3,i)

            ! If depth is near-zero, set momentum to zero
            if (h < depth_tolerance) then
                q1d(2:3,i) = 0.d0
                cycle 
            endif
            
            ! Apply friction source term only if in shallower water
            if (h <= friction_depth) then
                do nman = num_manning, 1, -1
                    if (aux1d(1,i) .lt. manning_break(nman)) then
                        coeff = manning_coefficient(nman)
                    endif
                enddo

                ! Calculate source term
                gamma = sqrt(hu**2 + hv**2) * (g * coeff**2) / h**(7.d0/3.d0)
                dgamma = 1.d0 + dt * gamma
                q1d(2, i) = q1d(2, i) / dgamma
                q1d(3, i) = q1d(3, i) / dgamma
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
    
            q1d(2,i) = q1d(2,i) * a(1,1) + q1d(3,i) * a(1,2)
            q1d(3,i) = q1d(2,i) * a(2,1) + q1d(3,i) * a(2,2)
        enddo
    endif

end subroutine src1d
