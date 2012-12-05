
subroutine qinit(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    
    use qinit_module, only: qinit_type,add_perturbation
    use geoclaw_module, only: eta_init
    
    implicit none
    
    ! Subroutine arguments
    integer, intent(in) :: maxmx,maxmy,meqn,mbc,mx,my,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy
    real(kind=8), intent(inout) :: q(meqn,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc)
    
    ! Locals
    integer :: i,j,m,layer_index
    
    ! Set flat state based on eta_init
    q = 0.d0
    do j=1,my
        do i=1,mx
            q(1,i,j) = max(0.d0,eta_init - aux(1,i,j)) 
        enddo
    enddo
    
    ! Add perturbation to initial conditions
    if (qinit_type > 0) then
        call add_perturbation(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    endif

    if (.false.) then
        open(23, file='fort.aux',status='unknown',form='formatted')
        print *,'Writing out aux arrays'
        print *,' '
        do j=1,my
            do i=1,mx
                write(23,*) i,j,(q(m,i,j),m=1,meqn)
            enddo
        enddo
        close(23)
    endif
    
end subroutine qinit
