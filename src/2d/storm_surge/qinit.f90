subroutine qinit(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    
    use multilayer_module, only: eta
    use geoclaw_module
    use hurricane_module

    implicit none
    
    ! Input
    integer, intent(in) :: maxmx,maxmy,meqn,mbc,mx,my,maux
    double precision, intent(in) :: xlower,ylower,dx,dy
    double precision, intent(inout) :: aux(maux,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc)

    ! Input/Output
    double precision, intent(inout) :: q(meqn,1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc)

    ! Locals
    integer :: i,j
    double precision :: x,y,wind_speed,tau,deta

    do i=1-mbc,mx+mbc
        x = xlower + (i-0.5d0)*dx
        do j=1-mbc,my+mbc
            y = ylower + (j-0.5d0)*dy
            q(1,i,j)=dmax1(0.d0,eta(1)-aux(1,i,j))
            q(2,i,j)=0.d0
            q(3,i,j)=0.d0
        enddo
    enddo
    
    ! Get and store hurricane wind field
    call hurricane_wind(maxmx,maxmy,maux,mbc,mx,my,xlower,ylower,dx,dy,-ramp_up_time,aux)
    
end subroutine qinit
