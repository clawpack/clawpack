subroutine qinit(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    
    use qinit_module
    use geoclaw_module
    
    implicit none
    
    ! Subroutine arguments
    integer, intent(in) :: maxmx,maxmy,meqn,mbc,mx,my,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy
    real(kind=8), intent(inout) :: q(meqn,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc)
    
    ! Locals
    integer :: i,j,m,layer_index,bottom_layer
    real(kind=8) :: xim,xip,xipc,ximc,xc,x,yim,yip,yjp,yjm,yjpc,yjmc,yc,y
    real(kind=8) :: dq
    
    ! Set flat state based on eta_init
    do i=1,mx
        do j=1,my
            ! This loop only executes if layers > 1, multiply by rho then
            do m=1,layers-1
                layer_index = 3*(m-1) + 1
                q(layer_index,i,j) = max(0.d0,eta_init(m)-eta_init(m+1)) * rho(m)
                q(layer_index+1,i,j) = 0.d0
                q(layer_index+2,i,j) = 0.d0
            enddo
            ! Take care of last layer
            bottom_layer = 3*layers-2
            q(bottom_layer,i,j) = max(0.d0,eta_init(m) - aux(1,i,j)) * rho(bottom_layer)
            q(bottom_layer+1,i,j) = 0.d0
            q(bottom_layer+2,i,j) = 0.d0
            
            ! Multiply by rho if this is a multi-layer run
            if (layers > 1) then
                q(bottom_layer,i,j) = q(bottom_layer,i,j) * rho(layers)
            endif
        enddo
    enddo
    
    ! Add perturbation to initial conditions
    call add_perturbation(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    
end subroutine qinit