subroutine qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    
    use qinit_module, only: qinit_type, add_perturbation
    
    implicit none
    
    ! Subroutine arguments
    integer, intent(in) :: meqn,mbc,mx,my,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    
    ! Locals
    integer :: i,j,m
    
    ! Set flat state based on eta_init
    q = 0.d0
    do i=1,mx
        do j=1,my
            ! Start with bottom layer and work up, set surface below for h
            eta_below = aux(1,i,j)
            do m=num_layers,1,-1
                layer_index = 3*(m-1) + 1
                q(layer_index,i,j) = max(0.d0,eta_init(m) - eta_below)
                eta_below = eta_init(m)
            enddo
        enddo
    enddo
    
    ! Add perturbation to initial conditions
    if (qinit_type > 0) then
        call add_perturbation(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    endif
    
end subroutine qinit