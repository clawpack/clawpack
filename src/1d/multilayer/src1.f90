subroutine src1(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt)

    use parameters_module, only: manning,g,rho

    implicit none

    ! Input parameters
    integer, intent(in) :: maxmx,meqn,mbc,mx,maux
    double precision, intent(in) :: xlower,dx,t,dt
    
    ! Output
    double precision, intent(inout) :: q(1-mbc:maxmx+mbc, meqn)
    double precision, intent(inout) :: aux(1-mbc:maxmx+mbc, maux)
    
    ! Locals
    integer :: i
    double precision, parameter :: TOLERANCE = 1d-30
    double precision :: gamma,dgamma,h,u
    
    ! Friction
    if (manning > TOLERANCE) then
        do i=1,mx
            h = q(i,3) / rho(2)
            u = q(i,4) / rho(2)
    
            if (h > TOLERANCE) then
                q(i,3) = 0.d0
                q(i,4) = 0.d0
            else
                gamma = u * g * manning**2 / h**(4/3)
                dgamma = 1.d0 + dt * gamma
                q(i,4) = q(i,4) / dgamma * rho(2)
            endif
        enddo
    endif
    
end subroutine src1
