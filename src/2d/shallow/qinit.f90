
subroutine qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    
    use qinit_module, only: qinit_type,add_perturbation
    use geoclaw_module, only: sea_level
    
    implicit none
    
    ! Subroutine arguments
    integer, intent(in) :: meqn,mbc,mx,my,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    
    ! Locals
    integer :: i,j,m
    
    ! Set flat state based on sea_level
    q = 0.d0
    forall(i=1:mx, j=1:my)
        q(1,i,j) = max(0.d0, sea_level - aux(1,i,j))
    end forall
    
    ! Add perturbation to initial conditions
    if (qinit_type > 0) then
        call add_perturbation(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
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
