subroutine qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    
    use qinit_module, only: qinit_type, add_perturbation
    use geoclaw_module, only: rho
    use multilayer_module, only: num_layers, eta_init
    
    implicit none
    
    ! Subroutine arguments
    integer, intent(in) :: meqn,mbc,mx,my,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    
    ! Locals
    integer :: i,j,m,layer_index
    real(kind=8) :: eta_below
    
    ! Set flat state based on eta_init
    q = 0.d0
    do j = 1, my
        do i = 1, mx
            ! Start with bottom layer and work up, set surface below for h
            eta_below = aux(1,i,j)
            do m = num_layers, 1, -1
                layer_index = 3 * (m-1) + 1
                q(layer_index,i,j) = max(0.d0,eta_init(m) - eta_below)
                eta_below = q(layer_index,i,j) + eta_below
                q(layer_index,i,j) = q(layer_index,i,j) * rho(m)
            enddo
        enddo
    enddo
    
    ! Add perturbation to initial conditions
    if (qinit_type > 0) then
        call add_perturbation(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    endif

    if (.false.) then
        open(23, file='fort.aux', status='unknown', form='formatted')
        print *, 'Writing out aux arrays'
        print *, ' '
        do j = 1, my
            do i = 1, mx
                write(23, *) i, j, (q(m,i,j), m=1, meqn)
            enddo
        enddo
        close(23)
    endif
    
end subroutine qinit