subroutine b4step1(maxmx,mbc,mx,meqn,q,xlower,dx,t,dt,maux,aux)

    use parameters_module

    implicit none
    
    ! Arguments
    integer, intent(in) :: maxmx,mbc,mx,meqn,maux
    double precision, intent(in) :: xlower,dx,t,dt
    double precision, intent(inout) :: q(1-mbc:maxmx+mbc,meqn)
    double precision, intent(inout) :: aux(1-mbc:maxmx+mbc,maux)
    
    ! Locals
    integer :: i,m,layer_index
    double precision :: h(2),u(2)
    double precision, parameter :: RICH_TOLERANCE = 0.95d0
    
    ! negative values to 0
    do m=1,4,2
        do i=1,mx
            if (q(i,m) < 0.d0) then
                q(i,m) = 0.d0
            endif
        enddo
    enddo
    
    ! Calculate wind field
    call set_wind(maxmx,mbc,mx,xlower,dx,t,aux(:,2))
    
    ! Calculate kappa
    do i=1,mx
        do m=1,2
            layer_index = 2*(m-1)
            h(m) = q(i,layer_index+1) / rho(m)
            if (h(m) < dry_tolerance) then
                u(m) = 0.d0
            else
                u(m) = q(i,layer_index+2) / q(i,layer_index+1)
            endif
        enddo
        aux(i,5) = (u(1)-u(2))**2 / (g*one_minus_r*sum(h))
        if ((aux(i,5) > RICH_TOLERANCE).and.(.not.h(2) < dry_tolerance)) then
            print "(a,i4,a,d16.8)","Hyperbolicity may have failed, kappa(",i,") = ",aux(i,5)
        endif
    enddo
    
end subroutine b4step1