subroutine src2(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt)
      
    use geoclaw_module, only: g => grav, coriolis_forcing, coriolis
    use geoclaw_module, only: friction_index, friction_forcing, friction_depth
    use geoclaw_module, only: num_layers, rho

    implicit none
    
    ! Input parameters
    integer, intent(in) :: meqn,mbc,mx,my,maux
    double precision, intent(in) :: xlower,ylower,dx,dy,t,dt
    
    ! Output
    double precision, intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    double precision, intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    ! Locals
    integer :: i, j, m, bottom_index, bottom_layer, layer_index
    logical :: found
    real(kind=8) :: h(num_layers), hu, hv, gamma, dgamma, y, fdt, a(2,2)

    ! Algorithm parameters
    ! Parameter controls when to zero out the momentum at a depth in the
    ! friction source term
    real(kind=8), parameter :: depth_tolerance = 1.0d-30

    ! Friction source term
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
                endif

                ! Apply friction source term only if in shallower water
                if (sum(h) <= friction_depth) then
                    ! Calculate source term
                    gamma = sqrt(hu**2 + hv**2) * g  &
                                    * aux(friction_index,i,j)**2 &
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
            y = ylower + (j - 0.5d0) * dy
            fdt = coriolis(y) * dt ! Calculate f dependent on coordinate system

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

end subroutine src2