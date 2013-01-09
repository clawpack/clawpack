subroutine src2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt)
<<<<<<< HEAD
      
    use geoclaw_module, only: g => grav, coriolis_forcing, coriolis
    use geoclaw_module, only: friction_index, friction_forcing, friction_depth

    implicit none
    
    ! Input parameters
    integer, intent(in) :: maxmx,maxmy,meqn,mbc,mx,my,maux
    double precision, intent(in) :: xlower,ylower,dx,dy,t,dt
    
    ! Output
    double precision, intent(inout) :: q(meqn,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc)
    double precision, intent(inout) :: aux(maux,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc)

    ! Locals
    integer :: i, j
    real(kind=8) :: h, hu, hv, gamma, dgamma, y, fdt, a(2,2)

    ! Algorithm parameters
    ! Parameter controls when to zero out the momentum at a depth in the
    ! friction source term
    real(kind=8), parameter :: depth_tolerance = 1.0d-30

    ! Friction source term
    if (friction_forcing) then
        do j=1,my
            do i=1,mx
                ! Extract appropriate momentum, also zero momentum in dry layers
                if (q(1,i,j) < depth_tolerance) then
                    q(2:3,i,j) = 0.d0
                else
                    ! Apply friction source term only if in shallower water
                    if (q(1,i,j) <= friction_depth) then
                        ! Calculate source term
                        gamma = sqrt(q(2,i,j)**2 + q(3,i,j)**2) * g     &   
                                * aux(friction_index,i,j)**2        &
                                / (q(1,i,j)**(7.d0/3.d0))
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
=======

  use geoclaw_module, only: g => grav, coriolis_forcing, coriolis,dry_manning_coefficient
  use geoclaw_module, only: friction_index, friction_forcing, friction_depth
  !  use geoclaw_module, only: num_layers, rho


  implicit none

  ! Input parameters
  integer, intent(in) :: maxmx,maxmy,meqn,mbc,mx,my,maux
  double precision, intent(in) :: xlower,ylower,dx,dy,t,dt

  ! Output
  double precision, intent(inout) :: q(meqn,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc)
  double precision, intent(inout) :: aux(maux,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc)


  ! Locals
  integer :: i, j
  real(kind=8) :: h, hu, hv, gamma, dgamma, y, fdt, a(2,2), coeff


  ! Algorithm parameters
  ! Parameter controls when to zero out the momentum at a depth in the
  ! friction source term
  real(kind=8), parameter :: depth_tolerance = 1.0d-30

  ! Friction source term
  if (friction_forcing) then
     do j=1,my
        do i=1,mx

           ! Extract depths
           h = q(1,i,j)

! orig code - does not match exactly
!!$           ! Extract appropriate momentum, also zero momentum in dry layers
!!$           if (h > depth_tolerance) then
!!$              ! Extract momentum components and exit loop
!!$              hu = q(2, i, j) 
!!$              hv = q(3, i, j) 
!!$              if (h <= friction_depth) then
!!$                 ! Calculate source term
!!$                 gamma = sqrt(hu**2 + hv**2) * g &
!!$                      * aux(friction_index,i,j)**2 &
!!$                      / (h**(7.d0/3.d0))
!!$                 dgamma = 1.d0 + dt * gamma
!!$                 q(2, i, j) = q(2, i, j) / dgamma
!!$                 q(3, i, j) = q(3, i, j) / dgamma
!!$              endif
!!$
!!$           else
!!$              ! Set almost dry layers momentum to zero
!!$              q(2, i, j) = 0.d0
!!$              q(3, i, j) = 0.d0
!!$           endif

                     ! Apply friction source term only if in shallower water

! 4-x code
               coeff = dry_manning_coefficient
               if (h.le.friction_depth) then
!                 # apply friction source term only in shallower water
                  hu=q(2,i,j)
                  hv=q(3,i,j)

                  if (h.lt.depth_tolerance) then
                     q(2,i,j)=0.d0
                     q(3,i,j)=0.d0
                  else
                     gamma= dsqrt(hu**2 + hv**2)*(g*coeff**2)/(h**(7.d0/3.d0))
                     dgamma=1.d0 + dt*gamma
                     q(2,i,j)= q(2,i,j)/dgamma
                     q(3,i,j)= q(3,i,j)/dgamma
                  endif
               endif

        enddo
     enddo
  endif
  ! End of friction source term

  ! Coriolis source term
  ! TODO: May want to remove the internal calls to coriolis as this could
  ! lead to slow downs.
  if (coriolis_forcing) then
     do j=1,my
        y = ylower + (j - 0.5d0) * dy
        fdt = coriolis(y) * dt ! Calculate f dependent on coordinate system

        ! Calculate matrix components
        a(1,1) = 1.d0 - 0.5d0 * fdt**2 + fdt**4 / 24.d0
        a(1,2) = fdt - fdt**3 / 6.d0
        a(2,1) = -fdt + fdt**3 / 6.d0
        a(2,2) = a(1,1)

        do i=1,mx
             
           q(2,i,j) = q(2, i, j) * a(1,1) + q(3, i, j) * a(1,2)
           q(3,i,j) = q(2, i, j) * a(2,1) + q(3, i, j) * a(2,2)

>>>>>>> mjb/omp-tests
        enddo
     enddo
  endif
  ! End of coriolis source term

end subroutine src2
